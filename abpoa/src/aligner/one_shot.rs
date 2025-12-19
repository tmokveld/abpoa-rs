use super::{Aligner, SequenceBatch, encode_sequences, to_i32};
use crate::encode::Alphabet;
use crate::result::{EncodedMsaResult, EncodedMsaView, MsaResult};
use crate::{Error, Result, sys};
use std::ptr;

impl Aligner {
    /// Run abPOA's one-shot MSA on a set of sequences and return the configured consensus/MSA output.
    ///
    /// This calls into abPOA `abpoa_msa`, so minimizer seeding, guide-tree
    /// partitioning, and progressive POA parameters will be used if enabled.
    pub fn msa(&mut self, batch: SequenceBatch<'_>) -> Result<MsaResult> {
        self.msa_one_shot_inner(batch, |abc, alphabet| unsafe {
            MsaResult::from_raw(abc, alphabet)
        })
    }

    /// Run abPOA's one-shot MSA on a set of sequences and return encoded consensus/MSA output.
    pub fn msa_encoded(&mut self, batch: SequenceBatch<'_>) -> Result<EncodedMsaResult> {
        self.msa_one_shot_inner(batch, |abc, _| unsafe { EncodedMsaResult::from_raw(abc) })
    }

    /// Run abPOA's one-shot MSA and return results as zero-copy encoded views.
    pub fn msa_view_encoded<'a>(
        &'a mut self,
        batch: SequenceBatch<'_>,
    ) -> Result<EncodedMsaView<'a>> {
        self.msa_one_shot_inner(batch, EncodedMsaView::new)
    }

    fn msa_one_shot_inner<T>(
        &mut self,
        batch: SequenceBatch<'_>,
        convert: impl FnOnce(*const sys::abpoa_cons_t, Alphabet) -> T,
    ) -> Result<T> {
        let seqs = batch.sequences();
        if seqs.is_empty() {
            let alphabet = self.alphabet();
            return Ok(convert(ptr::null(), alphabet));
        }
        let outputs = self.params.outputs();
        if outputs.is_empty() {
            return Err(Error::InvalidInput(
                "enable consensus and/or msa output to collect results".into(),
            ));
        }

        self.reset_cached_outputs()?;
        self.params
            .set_use_quality(batch.quality_weights().is_some());

        let alphabet = self.alphabet();
        let encoded = encode_sequences(seqs, alphabet);
        let mut seq_lens: Vec<i32> = encoded
            .iter()
            .map(|seq| to_i32(seq.len(), "sequence length exceeds i32"))
            .collect::<Result<_>>()?;

        // abpoa_msa only resets when abs->n_seq <= 0; clear it to avoid a redundant reset.
        // Safety: `self` owns the underlying aligner and its `abs` pointer is valid here.
        let abs_ptr = unsafe { (*self.as_mut_ptr()).abs };
        if abs_ptr.is_null() {
            return Err(Error::NullPointer("abpoa sequence container was null"));
        }
        self.graph_tracks_read_ids = true;
        self.read_id_count = 0;
        // Safety: `abs_ptr` is owned by this aligner; clearing n_seq forces abpoa_msa to reset.
        unsafe {
            (*abs_ptr).n_seq = 0;
        }

        let n_seq = to_i32(encoded.len(), "too many sequences for abpoa")?;
        let mut seq_ptrs: Vec<*mut u8> =
            encoded.iter().map(|seq| seq.as_ptr() as *mut u8).collect();

        let mut qual_ptrs: Vec<*mut i32> = batch
            .quality_weights()
            .map(|w| w.iter().map(|row| row.as_ptr() as *mut i32).collect())
            .unwrap_or_default();
        let qual_ptr = if qual_ptrs.is_empty() {
            ptr::null_mut()
        } else {
            qual_ptrs.as_mut_ptr()
        };

        let params_ptr = self.params.as_mut_ptr()?;
        let previous_out_gfa = unsafe { params_ptr.as_ref() }
            .ok_or(Error::NullPointer("abpoa parameters pointer was null"))?
            .out_gfa();
        // abpoa_msa routes all output through abpoa_output; if out_gfa is enabled it will skip
        // generating MSA/consensus unless it is cleared for the call.
        // Safety: `params_ptr` is uniquely owned and points to a live `abpoa_para_t`.
        unsafe {
            (*params_ptr).set_out_gfa(0);
        }
        // Safety: all pointers passed are valid for the duration of the call and match the
        // configured alphabet, abPOA will not keep them after returning.
        let status = unsafe {
            sys::abpoa_msa(
                self.as_mut_ptr(),
                params_ptr,
                n_seq,
                ptr::null_mut(),
                seq_lens.as_mut_ptr(),
                seq_ptrs.as_mut_ptr(),
                qual_ptr,
                ptr::null_mut(),
            )
        };
        // Safety: restore the caller-supplied value after the one-shot call completes.
        unsafe {
            (*params_ptr).set_out_gfa(previous_out_gfa);
        }
        if status != 0 {
            return Err(Error::Abpoa {
                func: "abpoa_msa",
                code: status,
            });
        }
        self.graph_tracks_read_ids = unsafe { params_ptr.as_ref() }
            .map(|raw| raw.use_read_ids() != 0)
            .unwrap_or(true);

        // Store the original sequences and optional names on the aligner for downstream graph
        // inspection helpers.
        self.store_batch_in_abs(&batch, 0, n_seq)?;

        let abc = unsafe { (*self.as_ptr()).abc };
        if abc.is_null() {
            return Err(Error::NullPointer(
                "abpoa returned a null consensus pointer",
            ));
        }

        Ok(convert(abc, alphabet))
    }
}
