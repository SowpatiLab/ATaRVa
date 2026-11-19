
import pytest
from ATARVA import decompose as mt
from ATARVA.decompose import (
    get_canonical_motif, kmp_search_non_overlapping,
    is_valid_repeat_block, get_repeat_components,
    convert_to_bitset, motif_decomposition,
)


@pytest.fixture(autouse=True)
def clean_registry():
    mt.cyclic_motif_registry.clear()

def test_canonical_motif():
    assert get_canonical_motif("CAG") == "AGC"
    assert get_canonical_motif("GCA") == "AGC"

def test_kmp_non_overlapping():
    assert kmp_search_non_overlapping("ATATAT", "AT") == [0, 2, 4]
    assert kmp_search_non_overlapping("ACATATATACAC", "AT") == [2, 4, 6]

def test_repeat_block_parsing():
    assert is_valid_repeat_block("(CAG)5")
    assert not is_valid_repeat_block("(CAG)")
    assert get_repeat_components("(CAG)5") == ("CAG", 5)

# def test_bitset():
#     l, r = convert_to_bitset("ACGT")
#     assert l.to01() == "0011"
#     assert r.to01() == "0101"

def expand(decomp):
    out = []
    for part in decomp.split("-"):
        if part.startswith("("):
            motif, n = get_repeat_components(part)
            out.append(motif * n)
        else:
            out.append(part)
    return "".join(out)

@pytest.mark.parametrize("seq,size", [
    ("CAGCAGCAGCAGCAG", 3),
    ("ATATATATAT", 2),
    ("AAAAAAAAAA", 1),
    ("CAGCAGCAGTTCAGCAG", 3),
    ("ACGTTGCA", 2),
])
def test_roundtrip(seq, size):
    decomp, non_rep = motif_decomposition(seq, size)
    print(f"\n{seq} (size {size}) -> {decomp}   non_repeat={non_rep}")   # <- added
    assert expand(decomp) == seq
    assert 0 <= non_rep <= 1