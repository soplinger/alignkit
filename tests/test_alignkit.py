import random
from pathlib import Path

import pytest

import alignkit as ak

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def _recompute(al, sc):
    """Score an alignment column by column, charging affine gap runs."""
    total, prev = 0, None
    for x, y in zip(al.aligned_a, al.aligned_b):
        if x == "-" or y == "-":
            kind = "a" if x == "-" else "b"
            total += sc.gap_extend if kind == prev else sc.gap_open
            prev = kind
        else:
            total += sc.score(x, y)
            prev = None
    return total


def test_parse_fasta_normalises_and_keeps_line_numbers():
    records, issues = ak.parse_fasta(">one\nacgt\nTT\n>two\nGGCC\n")
    assert [r.name for r in records] == ["one", "two"]
    assert records[0].seq == "ACGTTT"
    assert records[0].line == 1 and records[1].line == 4
    assert [i.severity for i in issues] == ["warning"]  # lowercase note only


def test_parse_fasta_reports_each_error_with_line_numbers():
    text = ">seq1\nacgt\n>seq2\nAC-GT\n>seq1\nGGGG\nTTTT\n>seq3\n"
    records, issues = ak.parse_fasta(text)
    errors = {i.line: i.message for i in issues if i.severity == "error"}
    assert set(errors) == {4, 5, 8}
    assert "invalid character '-' at column 3" == errors[4]
    assert "duplicate record name 'seq1'" in errors[5]
    assert "no sequence" in errors[8]
    assert [r.name for r in records] == ["seq1", "seq2"]


def test_parse_fasta_flags_data_before_header_and_dna_symbols():
    _, issues = ak.parse_fasta("ACGT\n>x\nACGTRU\n", alphabet="dna")
    assert issues[0].line == 1 and issues[0].severity == "error"
    by_line = {(i.line, i.severity) for i in issues}
    assert (3, "warning") in by_line  # R is an ambiguity code
    assert (3, "error") in by_line    # U is not a DNA base


def test_identical_sequences_score_two_per_base():
    al = ak.align("ACGTAC", "ACGTAC")
    st = al.stats()
    assert al.score == 12
    assert st["identity"] == 1
    assert ak.cigar(al.aligned_a, al.aligned_b) == "6="


def test_local_alignment_finds_embedded_motif():
    al = ak.align("TTTTGATTACATTTT", "CCCGATTACACCC", mode="local")
    assert al.score == 14
    assert al.aligned_a == al.aligned_b == "GATTACA"
    assert (al.a_start, al.a_end, al.b_start, al.b_end) == (4, 11, 3, 10)


def test_local_alignment_with_no_positive_segment_is_empty():
    al = ak.align("AAAA", "TTTT", mode="local")
    assert al.score == 0
    assert al.aligned_a == "" and al.aligned_b == ""


def test_semiglobal_leading_overhang_is_free():
    al = ak.align("GATTACA", "TTTTGATTACA", mode="semiglobal")
    assert al.score == 14
    assert al.aligned_b == "TTTTGATTACA"
    assert al.aligned_a == "----GATTACA"
    assert (al.a_start, al.a_end, al.b_start, al.b_end) == (0, 7, 0, 11)


def test_affine_gap_prefers_one_long_gap_run():
    al = ak.align("ACGTACGTACGT", "ACGTACGT")
    assert al.score == 9  # 8 matches (+16) and one 4-long gap (-4 open, -3 extend)
    st = al.stats()
    assert st["gap_opens"] == 1
    assert st["gap_columns"] == 4
    assert st["matches"] == 8


@pytest.mark.parametrize("mode", ["global", "local"])
def test_random_alignments_rescore_to_reported_score(mode):
    rng = random.Random(1234)
    sc = ak.Scoring()
    for _ in range(60):
        a = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 12)))
        b = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 12)))
        al = ak.align(a, b, mode=mode, scoring=sc)
        assert _recompute(al, sc) == al.score
        assert al.aligned_a.replace("-", "") == a[al.a_start:al.a_end]
        assert al.aligned_b.replace("-", "") == b[al.b_start:al.b_end]


def test_levenshtein_known_distances_and_cigar():
    assert ak.levenshtein("kitten", "sitting")[0] == 3
    assert ak.levenshtein("flaw", "lawn")[0] == 2
    assert ak.levenshtein("", "abc")[0] == 3
    assert ak.cigar("AC-T", "ACGA") == "2=1I1X"


def test_reverse_complement_handles_iupac_and_rejects_rna():
    assert ak.reverse_complement("ACGTN") == "NACGT"
    assert ak.reverse_complement("acgtr") == "YACGT"
    with pytest.raises(ak.AlignError):
        ak.reverse_complement("ACGU")


def test_strand_both_picks_reverse_complement_match():
    a = "GGATCCAAGTTT"
    b = ak.reverse_complement(a)
    forward = ak.align_strands(a, b, "global", None, both=False)
    best = ak.align_strands(a, b, "global", None, both=True)
    assert best.strand == "-"
    assert best.score == 24
    assert forward.score < best.score


def test_scoring_and_cell_limits_are_enforced():
    with pytest.raises(ak.AlignError):
        ak.Scoring(gap_open=1).check()
    with pytest.raises(ak.AlignError):
        ak.align("A" * 100, "A" * 100, max_cells=100)
    with pytest.raises(ak.AlignError):
        ak.align("", "ACGT")


def test_render_shows_match_bars():
    al = ak.align("ACGT", "ACGT")
    text = ak.render_blocks(al, "a", "b")
    assert "||||" in text


def test_pairs_identity_gate_exit_codes(tmp_path):
    good = tmp_path / "good.fasta"
    good.write_text(">a\nACGT\n>b\nACGT\n", encoding="utf-8")
    bad = tmp_path / "bad.fasta"
    bad.write_text(">a\nACGTACGTAC\n>b\nTTTTGGGGCC\n", encoding="utf-8")
    assert ak.main(["pairs", str(good), "--min-identity", "0.9"]) == 0
    assert ak.main(["pairs", str(bad), "--min-identity", "0.9"]) == 1


def test_align_command_reports_bad_input_with_exit_code_two(tmp_path, capsys):
    broken = tmp_path / "broken.fasta"
    broken.write_text(">a\nAC GT\n", encoding="utf-8")
    code = ak.main(["align", str(broken), str(broken)])
    assert code == 2
    assert "invalid character" in capsys.readouterr().err


def test_lint_command_exit_codes():
    assert ak.main(["lint", str(EXAMPLES / "panel.fasta")]) == 0
    assert ak.main(["lint", str(EXAMPLES / "broken.fasta")]) == 1
