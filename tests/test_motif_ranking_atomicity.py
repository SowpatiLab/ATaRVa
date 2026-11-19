import itertools

import pytest
from ATARVA import decompose as mt
from ATARVA.decompose import (
    get_atomic_motif, get_canonical_motif, get_most_frequent_motif,
    get_repeat_components, is_valid_repeat_block, motif_decomposition,
    reduce_block, refine_decomposition, register_and_get_motif,
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


def test_reduction_does_not_depend_on_a_large_motif_size():
    out, _ = refine_decomposition(["(CAGCAGCAG)2", "T"], 10, 40)
    assert out == ["(CAG)6", "T"]


@pytest.mark.parametrize("block,expected", [
    ("(ACACAC)2", "(AC)6"),
    ("(ATATAT)2", "(AT)6"),
])
def test_even_length_block_with_unequal_halves_is_reduced(block, expected):
    out, _ = refine_decomposition([block, "T"], 2, 40)
    assert out == [expected, "T"]


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


@pytest.mark.parametrize("size", [2, 3, 4, 6])
def test_single_block_homopolymer_is_atomic(size):
    decomp, _ = motif_decomposition("A" * 12, size)
    assert decomp == "(A)12"


def reference_atomic(motif):
    """Independent, brute-force definition: shortest prefix that tiles the motif."""
    for d in range(1, len(motif) + 1):
        if len(motif) % d == 0 and motif[:d] * (len(motif) // d) == motif:
            return motif[:d]
    return motif


def all_motifs(max_len, alphabet="ACGT"):
    for n in range(1, max_len + 1):
        for tup in itertools.product(alphabet, repeat=n):
            yield "".join(tup)


@pytest.mark.parametrize("motif,expected", [
    # homopolymers
    ("A", "A"), ("AA", "A"), ("AAAAAA", "A"), ("TTTTTTTTTT", "T"),
    # dinucleotides repeated
    ("ATAT", "AT"), ("ACACAC", "AC"), ("GTGTGTGT", "GT"),
    # trinucleotides repeated
    ("CAGCAG", "CAG"), ("CAGCAGCAG", "CAG"), ("CGGCGGCGGCGG", "CGG"),
    # longer units repeated
    ("GCGAGCAGCGAGCA", "GCGAGCA"), ("TATCTGATATCTGA", "TATCTGA"),
    # already atomic: must come back unchanged
    ("AC", "AC"), ("CAG", "CAG"), ("CAGT", "CAGT"), ("ACG", "ACG"),
    ("ACGTACGA", "ACGTACGA"), ("AAC", "AAC"), ("ACA", "ACA"),
    # N bases are ordinary characters
    ("NN", "N"), ("ACNACN", "ACN"),
    # empty
    ("", ""),
])
def test_get_atomic_motif_known_cases(motif, expected):
    assert get_atomic_motif(motif) == expected


def test_atomic_motif_of_a_rotation_has_the_same_length():
    # rotating a motif does not change how many bases the atomic unit has
    for motif in ("CAGCAG", "ACACAC", "AAAA", "CAGT", "GAGAGA"):
        n = len(get_atomic_motif(motif))
        for i in range(len(motif)):
            assert len(get_atomic_motif(motif[i:] + motif[:i])) == n, (motif, i)


def test_matches_brute_force_for_every_motif_up_to_length_7():
    bad = [m for m in all_motifs(7) if get_atomic_motif(m) != reference_atomic(m)]
    assert not bad, bad[:5]


def test_atomic_motif_properties_for_every_motif_up_to_length_6():
    for m in all_motifs(6):
        a = get_atomic_motif(m)
        assert m.startswith(a), m                       # is a prefix
        assert len(m) % len(a) == 0, m                  # divides the length
        assert a * (len(m) // len(a)) == m, m           # rebuilds the motif
        assert get_atomic_motif(a) == a, m              # idempotent: a is atomic


def test_atomic_motif_is_never_longer_than_the_motif():
    assert all(len(get_atomic_motif(m)) <= len(m) for m in all_motifs(6))


def test_binary_alphabet_exhaustive_to_length_12():
    # many periodic strings with a tiny alphabet; checks longer lengths cheaply
    for m in all_motifs(12, "AC"):
        assert get_atomic_motif(m) == reference_atomic(m), m


@pytest.mark.parametrize("block,expected", [
    ("(ACAC)3", "(AC)6"),
    ("(AAAAAA)2", "(A)12"),
    ("(AA)5", "(A)10"),
    ("(CAGCAG)4", "(CAG)8"),
    ("(CAGCAGCAG)2", "(CAG)6"),
    ("(GCGAGCAGCGAGCA)2", "(GCGAGCA)4"),
    ("(ACACAC)2", "(AC)6"),
    ("(ATATAT)2", "(AT)6"),
    ("(NN)3", "(N)6"),
])
def test_reduce_block_reduces_and_multiplies_the_count(block, expected):
    assert reduce_block(block) == expected


@pytest.mark.parametrize("block", ["(A)7", "(AC)8", "(CAG)5", "(CAGT)3", "(GCGAGCA)9"])
def test_reduce_block_leaves_atomic_blocks_alone(block):
    assert reduce_block(block) == block


@pytest.mark.parametrize("block,expected", [
    ("(ACAC)1", "(AC)2"),               # count of 1 still becomes 2 copies of AC
    ("(AA)100", "(A)200"),              # large counts
    ("(CAGCAG)1000", "(CAG)2000"),
])
def test_reduce_block_count_edge_cases(block, expected):
    assert reduce_block(block) == expected


def test_reduce_block_empty_motif_is_returned_unchanged():
    assert reduce_block("()5") == "()5"


@pytest.mark.parametrize("block", ["(AC)", "AC3", "(AC)x", "(AC", "", ")", "(AC)3.5"])
def test_reduce_block_rejects_malformed_blocks(block):
    with pytest.raises((ValueError, IndexError)):
        reduce_block(block)


def test_reduce_block_returns_a_valid_block_string():
    out = reduce_block("(ACAC)3")
    assert isinstance(out, str)
    assert is_valid_repeat_block(out)


def test_reduce_block_keeps_the_sequence_for_every_small_block():
    for motif in all_motifs(6):
        for count in (1, 2, 3, 7):
            block = f"({motif}){count}"
            out = reduce_block(block)
            m, n = get_repeat_components(out)
            assert m * n == motif * count, block                    # same sequence
            assert get_atomic_motif(m) == m, block                  # motif is atomic
            assert reduce_block(out) == out, block                  # idempotent
            assert len(m) * n == len(motif) * count, block          # same length


@pytest.fixture
def spy(monkeypatch):
    calls = []
    real = mt.reduce_block

    def wrapper(block):
        calls.append(block)
        return real(block)

    monkeypatch.setattr(mt, "reduce_block", wrapper)
    return calls


def test_refine_decomposition_reduces_every_block(spy):
    out, _ = refine_decomposition(["(ACAC)3", "T", "(CAGCAG)2"], 2, 30)
    assert "(ACAC)3" in spy and "(CAGCAG)2" in spy
    assert out == ["(AC)6", "T", "(CAG)4"]


def test_refine_decomposition_reduces_blocks_it_builds_from_plain_segments(spy):
    out, _ = refine_decomposition(["AAAAAA", "T"], 2, 7)
    assert "(AA)3" in spy                        # built inside refine, then reduced
    assert out == ["(A)6", "T"]


def test_motif_decomposition_reduces_a_single_block_result(spy):
    decomp, non_rep = motif_decomposition("A" * 12, 4)
    assert spy, "reduce_block was never called"
    assert decomp == "(A)12"
    assert non_rep == 0


def test_motif_decomposition_does_not_reduce_when_there_is_no_block(spy):
    decomp, non_rep = motif_decomposition("ACGT", 2)
    assert spy == []
    assert non_rep == 1


def test_malformed_block_inside_refine_is_kept_not_raised():
    out, _ = refine_decomposition(["(AC)", "GG"], 2, 10)    # no count after ')'
    assert out[0] == "(AC)"


@pytest.mark.parametrize("seq,size", [
    ("CTCCATCTCCATCTCCATCTCCATCTCCATCTCCATTTGGGGGCTA", 6),
    ("CAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCTCCCCCCCCCGA", 7),
    ("ACCACCACCACCACCGTTTTTTTTTTGC", 7),
    ("ATCTGAATCTGAATCTGAATCTGAATCTGAATCTGAATCTGACGAGAGAGAGAGAGACCTTACTTA", 5),
])
def test_roundtrip_when_a_block_is_followed_by_a_segment_starting_with_its_motif(seq, size):
    decomp, _ = motif_decomposition(seq, size)
    print(f"\n{seq} (size {size}) -> {decomp}")
    assert expand(decomp) == seq


def test_refine_keeps_the_remainder_of_a_partly_absorbed_segment():
    out, _ = refine_decomposition(["(G)4", "GCTA"], 1, 8)
    assert expand("-".join(out)) == "GGGG" + "GCTA"


def _random_cases(seed=7, n=400):
    import random
    rnd = random.Random(seed)

    def bases(k):
        return "".join(rnd.choice("ACGT") for _ in range(k))

    cases = []
    for _ in range(n):
        parts = []
        for _ in range(rnd.randint(1, 3)):
            parts.append(bases(rnd.randint(1, 6)) * rnd.randint(2, 9))
            if rnd.random() < .7:
                parts.append(bases(rnd.randint(1, 4)))
        cases.append(("".join(parts), rnd.randint(1, 8)))
    return cases


RANDOM_CASES = _random_cases()


def test_random_sequences_roundtrip_and_are_atomic():
    failures = []
    for seq, size in RANDOM_CASES:
        mt.cyclic_motif_registry.clear()
        decomp, non_rep = motif_decomposition(seq, size)
        ok = (expand(decomp) == seq and 0 <= non_rep <= 1
              and all(is_atomic(m) for m in block_motifs(decomp)))
        if not ok:
            failures.append((seq, size, decomp))
    assert not failures, f"{len(failures)} of {len(RANDOM_CASES)} failed, e.g. {failures[:3]}"