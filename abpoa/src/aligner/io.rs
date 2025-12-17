use super::Aligner;
use crate::encode::{Alphabet, decode_aa, decode_dna};
use crate::output::pog::PogWriteOptions;
use crate::params::OutputMode;
use crate::result::MsaResult;
use crate::{Error, Result, sys};
use libc;
use std::{
    env,
    ffi::CString,
    fs::{self, OpenOptions},
    io::Write,
    os::raw::c_char,
    os::unix::io::AsRawFd,
    path::{Path, PathBuf},
    process, slice,
    time::{SystemTime, UNIX_EPOCH},
};

impl Aligner {
    /// Write consensus sequences in FASTA form to a Rust `Write`.
    pub fn write_consensus_fasta(&mut self, writer: &mut impl Write) -> Result<()> {
        if self.graph_is_empty() {
            return Ok(());
        }

        let outputs = self.params.outputs();
        if !outputs.contains(OutputMode::CONSENSUS) {
            return Err(Error::InvalidInput(
                "consensus output is disabled; enable OutputMode::CONSENSUS in Parameters".into(),
            ));
        }

        self.reset_cached_outputs()?;
        let alphabet = self.alphabet();
        let params_ptr = self.params.as_mut_ptr()?;
        let consensus_needs_msa_rank = unsafe { params_ptr.as_ref() }
            .ok_or(Error::NullPointer("abpoa parameters pointer was null"))
            .map(|raw| raw.max_n_cons > 1 || raw.cons_algrm == sys::ABPOA_MF as i32)?;
        if consensus_needs_msa_rank && !self.graph_tracks_read_ids {
            return Err(Error::InvalidInput(
                "cannot write clustered consensus output (multiple consensus sequences or MostFrequent consensus) from a graph built without read ids; rebuild the graph with read ids enabled before adding sequences"
                    .into(),
            ));
        }
        if consensus_needs_msa_rank {
            self.ensure_msa_rank_buffer()?;
        }

        unsafe { sys::abpoa_generate_consensus(self.as_mut_ptr(), params_ptr) };
        let abc = unsafe { (*self.as_ptr()).abc };
        if abc.is_null() {
            return Err(Error::NullPointer(
                "abpoa returned a null consensus pointer",
            ));
        }

        let result = unsafe { MsaResult::from_raw(abc, alphabet) };
        let batch_index = unsafe { params_ptr.as_ref() }
            .ok_or(Error::NullPointer("abpoa parameters pointer was null"))?
            .batch_index;

        (|| -> Result<()> {
            for (idx, cluster) in result.clusters.iter().enumerate() {
                write!(writer, ">Consensus_sequence")?;
                if batch_index > 0 {
                    write!(writer, "_{}", batch_index)?;
                }
                if result.clusters.len() > 1 {
                    write!(writer, "_{} ", idx + 1)?;
                    for (read_idx, read_id) in cluster.read_ids.iter().enumerate() {
                        if read_idx > 0 {
                            write!(writer, ",")?;
                        }
                        write!(writer, "{}", read_id)?;
                    }
                }
                writeln!(writer)?;
                writeln!(writer, "{}", cluster.consensus)?;
            }
            Ok(())
        })()
    }

    /// Write consensus sequences in FASTQ form to a Rust `Write`.
    pub fn write_consensus_fastq(&mut self, writer: &mut impl Write) -> Result<()> {
        if self.graph_is_empty() {
            return Ok(());
        }

        let outputs = self.params.outputs();
        if !outputs.contains(OutputMode::CONSENSUS) {
            return Err(Error::InvalidInput(
                "consensus output is disabled; enable OutputMode::CONSENSUS in Parameters".into(),
            ));
        }

        self.reset_cached_outputs()?;
        let alphabet = self.alphabet();
        let params_ptr = self.params.as_mut_ptr()?;
        let consensus_needs_msa_rank = unsafe { params_ptr.as_ref() }
            .ok_or(Error::NullPointer("abpoa parameters pointer was null"))
            .map(|raw| raw.max_n_cons > 1 || raw.cons_algrm == sys::ABPOA_MF as i32)?;
        if consensus_needs_msa_rank && !self.graph_tracks_read_ids {
            return Err(Error::InvalidInput(
                "cannot write clustered consensus output (multiple consensus sequences or MostFrequent consensus) from a graph built without read ids; rebuild the graph with read ids enabled before adding sequences"
                    .into(),
            ));
        }
        if consensus_needs_msa_rank {
            self.ensure_msa_rank_buffer()?;
        }

        unsafe { sys::abpoa_generate_consensus(self.as_mut_ptr(), params_ptr) };
        let abc = unsafe { (*self.as_ptr()).abc };
        if abc.is_null() {
            return Err(Error::NullPointer(
                "abpoa returned a null consensus pointer",
            ));
        }

        let result = unsafe { MsaResult::from_raw(abc, alphabet) };
        let batch_index = unsafe { params_ptr.as_ref() }
            .ok_or(Error::NullPointer("abpoa parameters pointer was null"))?
            .batch_index;

        (|| -> Result<()> {
            for (cluster_idx, cluster) in result.clusters.iter().enumerate() {
                let expected_len = cluster.consensus.len();
                if cluster.phred.len() != expected_len {
                    return Err(Error::InvalidInput(
                        "consensus phred length did not match consensus sequence length".into(),
                    ));
                }

                write!(writer, "@Consensus_sequence")?;
                if batch_index > 0 {
                    write!(writer, "_{}", batch_index)?;
                }
                if result.clusters.len() > 1 {
                    write!(writer, "_{} ", cluster_idx + 1)?;
                    for (read_idx, read_id) in cluster.read_ids.iter().enumerate() {
                        if read_idx > 0 {
                            write!(writer, ",")?;
                        }
                        write!(writer, "{}", read_id)?;
                    }
                }
                writeln!(writer)?;
                writeln!(writer, "{}", cluster.consensus)?;

                write!(writer, "+Consensus_sequence")?;
                if batch_index > 0 {
                    write!(writer, "_{}", batch_index)?;
                }
                if result.clusters.len() > 1 {
                    write!(writer, "_{} ", cluster_idx + 1)?;
                    for (read_idx, read_id) in cluster.read_ids.iter().enumerate() {
                        if read_idx > 0 {
                            write!(writer, ",")?;
                        }
                        write!(writer, "{}", read_id)?;
                    }
                }
                writeln!(writer)?;
                writer.write_all(&cluster.phred)?;
                writeln!(writer)?;
            }
            Ok(())
        })()
    }

    /// Write RC-MSA output (optionally including consensus rows) to a Rust `Write`.
    pub fn write_msa_fasta(&mut self, writer: &mut impl Write) -> Result<()> {
        if self.graph_is_empty() {
            return Ok(());
        }
        if !self.params.outputs().contains(OutputMode::MSA) {
            return Err(Error::InvalidInput(
                "MSA output is disabled; enable OutputMode::MSA in Parameters".into(),
            ));
        }
        if !self.graph_tracks_read_ids {
            return Err(Error::InvalidInput(
                "cannot generate MSA from a graph built without read ids; rebuild the graph with read ids enabled before adding sequences"
                    .into(),
            ));
        }

        self.reset_cached_outputs()?;
        self.ensure_msa_rank_buffer()?;
        let alphabet = self.alphabet();
        let decode_row: fn(&[u8]) -> String = match alphabet {
            Alphabet::Dna => decode_dna,
            Alphabet::AminoAcid => decode_aa,
        };

        let params_ptr = self.params.as_mut_ptr()?;
        unsafe { sys::abpoa_generate_rc_msa(self.as_mut_ptr(), params_ptr) };
        let abc_ptr = unsafe { (*self.as_ptr()).abc };
        let abc = unsafe { abc_ptr.as_ref() }.ok_or(Error::NullPointer(
            "abpoa returned a null consensus pointer",
        ))?;
        let abs_ptr = unsafe { (*self.as_ptr()).abs };
        let abs = unsafe { abs_ptr.as_ref() }
            .ok_or(Error::NullPointer("abpoa sequence container was null"))?;
        let n_seq = abc.n_seq.max(0) as usize;
        let msa_len = abc.msa_len.max(0) as usize;

        if abc.msa_len <= 0 {
            Ok(())
        } else if abc.msa_base.is_null() {
            Err(Error::NullPointer("abpoa returned a null msa buffer"))
        } else {
            (|| -> Result<()> {
                for idx in 0..n_seq {
                    let name = format_msa_name(abs, idx)?;
                    let row_ptr = unsafe { *abc.msa_base.add(idx) };
                    if row_ptr.is_null() {
                        return Err(Error::NullPointer("abpoa returned a null msa row"));
                    }
                    let row = unsafe { slice::from_raw_parts(row_ptr, msa_len) };
                    writeln!(writer, ">{}", name)?;
                    writeln!(writer, "{}", decode_row(row))?;
                }

                if abc.n_cons > 0 {
                    for cons_idx in 0..abc.n_cons.max(0) as usize {
                        let row_ptr = unsafe { *abc.msa_base.add(n_seq + cons_idx) };
                        if row_ptr.is_null() {
                            return Err(Error::NullPointer(
                                "abpoa returned a null consensus msa row",
                            ));
                        }
                        let ids = consensus_read_ids(abc, cons_idx);
                        let row = unsafe { slice::from_raw_parts(row_ptr, msa_len) };
                        write!(writer, ">Consensus_sequence")?;
                        if abc.n_cons > 1 {
                            write!(writer, "_{} ", cons_idx + 1)?;
                            for (idx, read_id) in ids.iter().enumerate() {
                                if idx > 0 {
                                    write!(writer, ",")?;
                                }
                                write!(writer, "{}", read_id)?;
                            }
                        }
                        writeln!(writer)?;
                        writeln!(writer, "{}", decode_row(row))?;
                    }
                }
                Ok(())
            })()
        }
    }

    /// Write GFA output to the given path.
    pub fn write_gfa_to_path<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        if self.graph_is_empty() {
            return Ok(());
        }
        if !self.graph_tracks_read_ids {
            return Err(Error::InvalidInput(
                "cannot generate GFA from a graph built without read ids; rebuild the graph with read ids enabled before adding sequences"
                    .into(),
            ));
        }

        let outputs = self.params.outputs();
        if !outputs.contains(OutputMode::MSA) {
            return Err(Error::InvalidInput(
                "GFA output requires MSA output; enable OutputMode::MSA in Parameters".into(),
            ));
        }

        self.reset_cached_outputs()?;
        let abs_ptr = unsafe { (*self.as_ptr()).abs };
        let has_consensus = has_consensus_sequence(abs_ptr)?;
        let wants_consensus = outputs.contains(OutputMode::CONSENSUS);
        let desired_outputs = if wants_consensus && !has_consensus {
            OutputMode::CONSENSUS | OutputMode::MSA
        } else {
            // Avoid regenerating consensus if the graph already carries one from a restore.
            OutputMode::MSA
        };
        if desired_outputs.contains(OutputMode::CONSENSUS) && self.consensus_needs_msa_rank()? {
            self.ensure_msa_rank_buffer()?;
        }
        let previous_outputs = outputs;
        self.params.set_outputs_for_call(desired_outputs);
        let write_result = (|| -> Result<()> {
            let params_ptr = self.params.as_mut_ptr()?;

            let file = OpenOptions::new()
                .read(true)
                .write(true)
                .create(true)
                .truncate(true)
                .open(path.as_ref())?;
            // Safety: `dup` creates an owned fd for `fdopen`/`fclose` without double-closing `file`
            let dup_fd = unsafe { libc::dup(file.as_raw_fd()) };
            if dup_fd < 0 {
                return Err(std::io::Error::last_os_error().into());
            }

            // Safety: `dup_fd` is a valid fd opened for read/write; `fdopen` takes ownership
            let fp = unsafe { libc::fdopen(dup_fd, c"w+".as_ptr() as *const c_char) };

            if fp.is_null() {
                unsafe { libc::close(dup_fd) };
                return Err(Error::NullPointer("failed to open FILE for GFA output"));
            }

            let write_result = {
                unsafe {
                    sys::abpoa_generate_gfa(self.as_mut_ptr(), params_ptr, fp as *mut sys::FILE);
                    libc::fflush(fp);
                }
                Ok(())
            };

            unsafe { libc::fclose(fp) };
            // abpoa_generate_gfa may allocate consensus buffers when out_cons is enabled
            unsafe { sys::abpoa_clean_msa_cons(self.as_mut_ptr()) };

            write_result
        })();
        self.params.set_outputs_for_call(previous_outputs);
        write_result
    }

    /// Write GFA output to a Rust `Write`.
    ///
    /// This is a convenience adapter around [`Self::write_gfa_to_path`] and uses a temporary file.
    pub fn write_gfa(&mut self, writer: &mut impl Write) -> Result<()> {
        if self.graph_is_empty() {
            return Ok(());
        }

        let path = temp_gfa_path()?;
        let write_result = self.write_gfa_to_path(&path);
        let read_result = match write_result {
            Ok(()) => {
                let contents = fs::read(&path)?;
                writer.write_all(&contents)?;
                Ok(())
            }
            Err(err) => Err(err),
        };

        let _ = fs::remove_file(&path);
        read_result
    }

    /// Write a GraphViz view of the current partial order graph (POG) to `path`.
    ///
    /// Supported extensions:
    /// - `.dot`: write GraphViz DOT
    /// - `.png` / `.pdf`: write an image using GraphViz `dot`
    ///
    /// When writing `.png`/`.pdf`, a DOT file is written next to `path` at `path + ".dot"`
    /// (mirroring upstream).
    ///
    /// `options` can be:
    /// - [`crate::PogDotOptions`] (or `&crate::PogDotOptions`) to use the Rust DOT emitter
    /// - [`PogWriteOptions::Upstream`] to call upstream `abpoa_dump_pog` (uses `system()` and may
    ///   `exit()` on failure).
    pub fn write_pog_to_path<P: AsRef<Path>>(
        &mut self,
        path: P,
        options: impl Into<PogWriteOptions>,
    ) -> Result<()> {
        if self.graph_is_empty() {
            return Ok(());
        }

        let path = path.as_ref();
        let ext = path
            .extension()
            .and_then(|ext| ext.to_str())
            .ok_or(Error::InvalidInput(
                "graph dump path must end with .dot, .png, or .pdf".into(),
            ))?;

        match options.into() {
            PogWriteOptions::Rust(dot_options) => {
                self.ensure_topological()?;
                let graph = self.graph()?;

                if ext == "dot" {
                    let mut file = OpenOptions::new()
                        .write(true)
                        .create(true)
                        .truncate(true)
                        .open(path)?;
                    return crate::output::dot::write_pog_dot(
                        &graph,
                        self.alphabet(),
                        &mut file,
                        &dot_options,
                    );
                }

                if ext != "png" && ext != "pdf" {
                    return Err(Error::InvalidInput(
                        "graph dump path must end with .dot, .png, or .pdf".into(),
                    ));
                }

                let mut dot_path = std::ffi::OsString::from(path.as_os_str());
                dot_path.push(".dot");
                let dot_path = PathBuf::from(dot_path);

                {
                    let mut file = OpenOptions::new()
                        .write(true)
                        .create(true)
                        .truncate(true)
                        .open(&dot_path)?;
                    crate::output::dot::write_pog_dot(
                        &graph,
                        self.alphabet(),
                        &mut file,
                        &dot_options,
                    )?;
                }

                let status = process::Command::new("dot")
                    .arg(format!("-T{ext}"))
                    .arg("-o")
                    .arg(path)
                    .arg(&dot_path)
                    .status()
                    .map_err(|err| {
                        std::io::Error::other(format!("failed to execute `dot`: {err}"))
                    })?;

                if !status.success() {
                    return Err(std::io::Error::other(format!(
                        "`dot` exited with status {status}"
                    ))
                    .into());
                }

                Ok(())
            }
            PogWriteOptions::Upstream => {
                if ext != "png" && ext != "pdf" {
                    return Err(Error::InvalidInput(
                        "upstream graph dump requires a .png or .pdf path".into(),
                    ));
                }

                let path_str = path.to_str().ok_or(Error::InvalidInput(
                    "graph dump path must be valid UTF-8".into(),
                ))?;
                let c_path = CString::new(path_str).map_err(|_| {
                    Error::InvalidInput("graph dump path cannot contain null bytes".into())
                })?;

                let params_ptr = self.params.as_mut_ptr()?;
                let previous = unsafe { (*params_ptr).out_pog };
                // Safety: `c_path` is a valid C string; `strdup` allocates an owned copy.
                let new = unsafe { libc::strdup(c_path.as_ptr()) };
                if new.is_null() {
                    return Err(Error::NullPointer("failed to store graph dump path"));
                }
                unsafe {
                    (*params_ptr).out_pog = new;
                }

                // Safety: `self`/`self.params` are valid for the duration of the call.
                unsafe { sys::abpoa_dump_pog(self.as_mut_ptr(), params_ptr) };

                unsafe {
                    libc::free((*params_ptr).out_pog as *mut libc::c_void);
                    (*params_ptr).out_pog = previous;
                }

                Ok(())
            }
        }
    }
}

fn has_consensus_sequence(abs: *const sys::abpoa_seq_t) -> Result<bool> {
    let Some(abs) = (unsafe { abs.as_ref() }) else {
        return Err(Error::NullPointer("abpoa sequence container was null"));
    };
    if abs.name.is_null() {
        return Ok(false);
    }
    let count = abs.n_seq.max(0) as usize;
    for idx in 0..count {
        let name = unsafe { abs.name.add(idx).as_ref() }
            .ok_or(Error::NullPointer("abpoa sequence name pointer was null"))?;
        if name.l > 0 && !name.s.is_null() {
            let bytes = unsafe { slice::from_raw_parts(name.s as *const u8, name.l as usize) };
            if bytes.starts_with(b"Consensus_sequence") {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn sequence_name(abs: &sys::abpoa_seq_t, idx: usize) -> Result<String> {
    if abs.name.is_null() {
        return Ok(format!("Seq_{}", idx + 1));
    }
    let kstr = unsafe { abs.name.add(idx).as_ref() }
        .ok_or(Error::NullPointer("abpoa sequence name pointer was null"))?;
    if kstr.l > 0 && !kstr.s.is_null() {
        // Safety: `s` points to a buffer at least `l` bytes long when `l > 0`.
        let bytes = unsafe { slice::from_raw_parts(kstr.s as *const u8, kstr.l as usize) };
        Ok(String::from_utf8_lossy(bytes).into_owned())
    } else {
        Ok(format!("Seq_{}", idx + 1))
    }
}

fn sequence_is_rc(abs: &sys::abpoa_seq_t, idx: usize) -> Result<bool> {
    if abs.is_rc.is_null() {
        return Ok(false);
    }
    // Safety: `is_rc` is sized to `n_seq` by abPOA; `idx` is derived from that bound.
    let flag = unsafe { *abs.is_rc.add(idx) };
    Ok(flag != 0)
}

fn format_msa_name(abs: &sys::abpoa_seq_t, idx: usize) -> Result<String> {
    let mut name = sequence_name(abs, idx)?;
    if sequence_is_rc(abs, idx)? {
        name.push_str("_reverse_complement");
    }
    Ok(name)
}

fn consensus_read_ids(abc: &sys::abpoa_cons_t, idx: usize) -> Vec<i32> {
    if abc.clu_n_seq.is_null() || abc.clu_read_ids.is_null() {
        return Vec::new();
    }
    // Safety: read count and pointers are filled by abPOA when consensus is generated.
    let count = unsafe { *abc.clu_n_seq.add(idx) }.max(0) as usize;
    if count == 0 {
        return Vec::new();
    }
    let ids_ptr = unsafe { *abc.clu_read_ids.add(idx) };
    if ids_ptr.is_null() {
        return Vec::new();
    }
    unsafe { slice::from_raw_parts(ids_ptr, count) }.to_vec()
}

fn temp_gfa_path() -> std::io::Result<PathBuf> {
    let mut path = env::temp_dir();
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    path.push(format!("abpoa_gfa_{}_{}.gfa", process::id(), timestamp));
    Ok(path)
}
