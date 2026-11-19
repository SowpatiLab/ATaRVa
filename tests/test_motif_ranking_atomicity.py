
import pytest
from ATARVA import decompose as mt
from ATARVA.decompose import (
    get_canonical_motif, get_most_frequent_motif, get_repeat_components,
    motif_decomposition, refine_decomposition, register_and_get_motif,
)


@pytest.fixture(autouse=True)
def clean_registry():
    mt.cyclic_motif_registry.clear()


def expand(decomp):
    out = []
    for part in decomp.split("-"):
        if part.startswith("("):
            motif, n = get_repeat_components(part)
            out.append(motif * n)
        else:
            out.append(part)
    return "".join(out)


def block_motifs(decomp):
    return [get_repeat_components(p)[0] for p in decomp.split("-") if p.startswith("(")]


def is_atomic(motif):
    """True if the motif is not a repeat of a shorter string (ACAC -> False)."""
    return motif not in (motif + motif)[1:-1]


@pytest.mark.parametrize("seq,size,expected", [
    ("AAAAAAAAAACCCGGGT",   1, ("A", "C")),    
    ("GGGGGAAAAATTTTTTTCC", 1, ("T", "G")),    
])
def test_motifs_are_ranked_by_count(seq, size, expected):
    result = get_most_frequent_motif(seq, size, "")
    print(f"\n{seq} (size {size}) -> {result}")
    assert result == expected


@pytest.mark.parametrize("seq,first,second", [
    ("ACACACGTGTGT", "AC", "GT"),
    ("GTGTGTACACAC", "GT", "AC"),
])
def test_tie_is_broken_by_first_appearance(seq, first, second):
    assert get_most_frequent_motif(seq, 2, "") == (first, second)


def test_tie_for_second_place_goes_to_the_motif_seen_first():
    assert get_most_frequent_motif("AAAAAAAAAACCCGGGT", 1, "")[1] == "C"
    assert get_most_frequent_motif("AAAAAAAAAAGGGCCCT", 1, "")[1] == "G"


def test_only_the_top_two_are_returned():
    result = get_most_frequent_motif("AAAAACCCCGGGT", 1, "")
    assert len(result) == 2
    assert result == ("A", "C")


def test_sequence_shorter_than_motif_size():
    assert get_most_frequent_motif("ACGT", 5, "") == ("ACGT", None)


def test_counting_uses_a_sliding_window():
    assert get_most_frequent_motif("ACACAC", 2, "") == ("AC", "CA")


@pytest.mark.parametrize("seq,motif,expected", [
    ("ACACACACACGTGT", "AC", ("AC", "GT")),
    ("ACACACACACGTGT", "CA", ("CA", "GT")),      
    ("CAGCAGCAGCAGCAGTTTTT", "CAG", ("CAG", "TTT")),
])
def test_given_motif_collapses_its_rotations(seq, motif, expected):
    assert get_most_frequent_motif(seq, len(motif), motif) == expected


@pytest.mark.parametrize("seq,size,expected", [
    ("ACACACACACGTGT", 2, ("AC", "CA")),
    ("CAGCAGCAGCAGCAGTTTTT", 3, ("CAG", "AGC")),
])
def test_without_a_given_motif_secondary_can_be_a_rotation_of_primary(seq, size, expected):
    """Documents CURRENT behaviour: rotations are only collapsed when a motif is given.
    If you want the second-ranked motif to always be a different motif, change this."""
    primary, secondary = get_most_frequent_motif(seq, size, "")
    assert (primary, secondary) == expected
    assert get_canonical_motif(primary) == get_canonical_motif(secondary)


def test_primary_motif_is_registered():
    primary, _ = get_most_frequent_motif("CAGCAGCAGCAGCAGTTTTT", 3, "CAG")
    assert mt.cyclic_motif_registry[get_canonical_motif("CAG")] == primary


def test_registry_keeps_first_rotation_seen():
    assert [register_and_get_motif(m) for m in ("CAG", "AGC", "GCA", "CAG")] == ["CAG"] * 4
    mt.cyclic_motif_registry.clear()
    assert [register_and_get_motif(m) for m in ("GCA", "CAG", "AGC")] == ["GCA"] * 3


def test_registry_ignores_empty_motif():
    assert register_and_get_motif("") == ""
    assert mt.cyclic_motif_registry == {}


@pytest.mark.parametrize("block,size,expected", [
    ("(AA)5",        1, "(A)10"),          # single-base unit: counts multiply
    ("(AAA)2",       1, "(A)6"),
    ("(ATAT)3",      2, "(AT)6"),          # even length, halves equal
    ("(CAGCAG)4",    3, "(CAG)8"),
    ("(AGGAGG)3",    3, "(AGG)6"),
    ("(CAGCAGCAG)2", 3, "(CAG)6"),         # odd length, divisible by motif_size
])
def test_block_is_reduced_to_its_atomic_unit(block, size, expected):
    out, _ = refine_decomposition([block, "T"], size, 40)
    print(f"\n{block} (size {size}) -> {out}")
    assert out == [expected, "T"]
    assert expand("-".join(out)) == expand(block + "-T")


@pytest.mark.parametrize("block,size", [
    ("(CAGT)3", 4),                        
    ("(CAG)5",  3),
    ("(AT)8",   2),
])
def test_atomic_block_is_left_alone(block, size):
    out, _ = refine_decomposition([block, "T"], size, 40)
    assert out == [block, "T"]


def test_block_longer_than_motif_size_is_not_reduced_when_size_too_large():
    """Documents CURRENT behaviour: a motif no longer than motif_size is kept as is."""
    out, _ = refine_decomposition(["(CAGCAGCAG)2", "T"], 10, 40)
    assert out == ["(CAGCAGCAG)2", "T"]


@pytest.mark.xfail(strict=True, reason="even-length motif whose halves differ "
                   "((ACACAC), (ATATAT)) is never reduced to its atomic unit")
@pytest.mark.parametrize("block,expected", [
    ("(ACACAC)2", "(AC)6"),
    ("(ATATAT)2", "(AT)6"),
])
def test_even_length_block_with_unequal_halves_is_reduced(block, expected):
    out, _ = refine_decomposition([block, "T"], 2, 40)
    assert out == [expected, "T"]


@pytest.mark.xfail(strict=True, reason="atomic reduction of an odd-length block "
                   "only happens when its length is divisible by motif_size")
def test_reduction_does_not_depend_on_motif_size():
    out, _ = refine_decomposition(["(CAGCAGCAG)2", "T"], 2, 40)
    assert out == ["(CAG)6", "T"]


@pytest.mark.parametrize("seq,size,expected", [
    ("AT" * 8,        4, "(AT)8"),         
    ("AT" * 8,        6, "(AT)8"),
    ("AC" * 8,        6, "(AC)8"),
    ("CAG" * 6,       6, "(CAG)6"),
    ("CAG" * 6,       9, "(CAG)6"),
    ("A" * 12,        1, "(A)12"),
    ("ACACACACACACGGG", 2, "(AC)6-GGG"),
    ("ATATATATATATCCC", 2, "(AT)6-CCC"),
])
def test_decomposition_reports_atomic_motif(seq, size, expected):
    decomp, _ = motif_decomposition(seq, size)
    print(f"\n{seq} (size {size}) -> {decomp}")
    assert decomp == expected
    assert expand(decomp) == seq


@pytest.mark.parametrize("seq,size", [
    ("AT" * 8, 4), ("AC" * 8, 6), ("CAG" * 6, 6), ("CAG" * 6, 9), ("A" * 12, 1),
    ("ATATATATATATTT", 4), ("CAGCAGCAGCAGCAGCAGTT", 6),
])
def test_every_repeat_block_motif_is_atomic(seq, size):
    decomp, _ = motif_decomposition(seq, size)
    assert all(is_atomic(m) for m in block_motifs(decomp)), decomp


@pytest.mark.xfail(strict=True, reason="a sequence that is ONE repeat block skips "
                   "refine_decomposition (len(fseq)==1), so (AAAAAA)2 is not "
                   "reduced to (A)12")
@pytest.mark.parametrize("size", [2, 3, 4, 6])
def test_single_block_homopolymer_is_atomic(size):
    decomp, _ = motif_decomposition("A" * 12, size)
    assert decomp == "(A)12"