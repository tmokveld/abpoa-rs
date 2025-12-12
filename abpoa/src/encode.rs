//! Sequence encoding helpers mirroring abPOA's built-in tables

/// Supported alphabets for encoding/decoding
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Alphabet {
    Dna,
    AminoAcid,
}

/// Encode a DNA sequence (A/C/G/T/U/N) into abPOA's 0-4 integer alphabet
/// Unknown bases are coerced to `N` (`4`), matching `ab_nt4_table`
pub fn encode_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter().copied().map(encode_dna_base).collect()
}

/// Decode abPOA's DNA integer alphabet back into bases using `ab_nt256_table`
pub fn decode_dna(codes: &[u8]) -> String {
    codes.iter().map(|&c| decode_dna_code(c)).collect()
}

/// Encode an amino acid sequence into abPOA's 0-26 alphabet
/// Unknown residues are mapped to `*` (`26`), matching `ab_aa26_table`
pub fn encode_aa(seq: &[u8]) -> Vec<u8> {
    seq.iter().copied().map(encode_aa_residue).collect()
}

/// Decode abPOA's amino acid alphabet back into residue letters using `ab_aa256_table`
pub fn decode_aa(codes: &[u8]) -> String {
    codes.iter().map(|&c| decode_aa_code(c)).collect()
}

fn encode_dna_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        b'N' => 4,
        _ => 4,
    }
}

pub(crate) fn decode_dna_code(code: u8) -> char {
    match code {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        4 => 'N',
        // abPOA uses `m` (5) as the gap sentinel in MSA output; 27 also maps to '-'
        5 | 27 => '-',
        _ => 'N',
    }
}

fn encode_aa_residue(residue: u8) -> u8 {
    match residue.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        b'N' => 4,
        b'B' => 5,
        b'D' => 6,
        b'E' => 7,
        b'F' => 8,
        b'H' => 9,
        b'I' => 10,
        b'J' => 11,
        b'K' => 12,
        b'L' => 13,
        b'M' => 14,
        b'O' => 15,
        b'P' => 16,
        b'Q' => 17,
        b'R' => 18,
        b'S' => 19,
        b'U' => 20,
        b'V' => 21,
        b'W' => 22,
        b'X' => 23,
        b'Y' => 24,
        b'Z' => 25,
        b'*' => 26,
        _ => 26,
    }
}

pub(crate) fn decode_aa_code(code: u8) -> char {
    match code {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        4 => 'N',
        5 => 'B',
        6 => 'D',
        7 => 'E',
        8 => 'F',
        9 => 'H',
        10 => 'I',
        11 => 'J',
        12 => 'K',
        13 => 'L',
        14 => 'M',
        15 => 'O',
        16 => 'P',
        17 => 'Q',
        18 => 'R',
        19 => 'S',
        20 => 'U',
        21 => 'V',
        22 => 'W',
        23 => 'X',
        24 => 'Y',
        25 => 'Z',
        26 => '*',
        // abPOA uses `m` (27) as the gap sentinel in RC-MSA output
        27 => '-',
        _ => '*',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_round_trips_known_bases() {
        let seq = b"AcgTNu-?";
        let encoded = encode_dna(seq);
        assert_eq!(encoded, vec![0, 1, 2, 3, 4, 3, 4, 4]);
        let decoded = decode_dna(&encoded);
        assert_eq!(decoded, "ACGTNTNN");
    }

    #[test]
    fn dna_accepts_u_as_thymine() {
        let encoded = encode_dna(b"uUtT");
        assert_eq!(encoded, vec![3, 3, 3, 3]);
    }

    #[test]
    fn amino_acid_round_trips_known_residues() {
        let seq = b"ACGTNBDEFHIJKLMOPQRSUVWXYZ*";
        let encoded = encode_aa(seq);
        assert_eq!(
            encoded,
            (0u8..=26u8).collect::<Vec<u8>>(),
            "mapping should follow ab_aa26_table ordering"
        );
        let decoded = decode_aa(&encoded);
        assert_eq!(decoded, "ACGTNBDEFHIJKLMOPQRSUVWXYZ*");
    }

    #[test]
    fn amino_acid_unknown_defaults_to_stop() {
        let encoded = encode_aa(b"?");
        assert_eq!(encoded, vec![26]);
        let decoded = decode_aa(&[250]);
        assert_eq!(decoded, "*");
    }

    #[test]
    fn decodes_gap_sentinels() {
        assert_eq!(decode_dna(&[0, 1, 2, 3, 4, 5, 27]), "ACGTN--");
        assert_eq!(decode_aa(&[0, 26, 27, 250]), "A*-*");
    }
}
