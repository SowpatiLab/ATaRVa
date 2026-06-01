import re
from dataclasses import dataclass

CIGAR_OPS = {
    'M': 0, 'I': 1, 'D': 2, 'N': 3,
    'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8
}

# ops that consume reference
REF_CONSUMING   = {0, 2, 3, 7, 8}
# ops that consume query
QUERY_CONSUMING = {0, 1, 4, 7, 8}


@dataclass
class MDValidationResult:
    is_valid:         bool
    errors:           list[str]
    cigar_ref_len:    int     # ref bases consumed by CIGAR
    md_ref_len:       int     # ref bases described by MD tag
    md_matches:       int     # matched bases in MD
    md_mismatches:    int     # mismatched bases in MD
    md_deletions:     int     # deleted bases in MD

    def __str__(self) -> str:
        status = '✅ VALID' if self.is_valid else '❌ INVALID'
        lines  = [
            f'{status}',
            f'  CIGAR ref len : {self.cigar_ref_len}',
            f'  MD    ref len : {self.md_ref_len}',
            f'  MD matches    : {self.md_matches}',
            f'  MD mismatches : {self.md_mismatches}',
            f'  MD deletions  : {self.md_deletions}',
        ]
        if self.errors:
            lines.append('  Errors:')
            for err in self.errors:
                lines.append(f'    ✗ {err}')
        return '\n'.join(lines)


def _parse_cigar(cigar_str: str) -> list[tuple[int, int]]:
    """Parse CIGAR string into (op_int, length) tuples."""
    return [
        (CIGAR_OPS[op], int(length))
        for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_str)
    ]


def _parse_md_tag(md_tag: str) -> list[tuple[str, str | int]]:
    """
    Parse MD tag into list of (type, value) tokens.
    Types:
        'match'    → int  number of matching bases
        'mismatch' → str  mismatched reference base
        'deletion' → str  deleted reference bases (^ACGT)
    """
    tokens = []
    for match in re.finditer(r'(\d+)|(\^[A-Za-z]+)|([A-Za-z])', md_tag):
        if match.group(1) is not None:
            tokens.append(('match',    int(match.group(1))))
        elif match.group(2) is not None:
            tokens.append(('deletion', match.group(2)[1:]))  # strip ^
        else:
            tokens.append(('mismatch', match.group(3)))
    return tokens


def _cigar_ref_len(cigartuples: list[tuple]) -> int:
    """Total reference bases consumed by CIGAR."""
    return sum(length for op, length in cigartuples if op in REF_CONSUMING)


def _md_ref_len(md_tokens: list[tuple]) -> tuple[int, int, int, int]:
    """
    Total ref bases described by MD tag.
    Returns (total, matches, mismatches, deletions)
    """
    matches    = sum(v         for t, v in md_tokens if t == 'match')
    mismatches = sum(1         for t, v in md_tokens if t == 'mismatch')
    deletions  = sum(len(v)    for t, v in md_tokens if t == 'deletion')
    total      = matches + mismatches + deletions
    return total, matches, mismatches, deletions


def validate_md_tag(md_tag:    str,
                    cigar_str: str) -> MDValidationResult:
    """
    Validate MD tag consistency against CIGAR string.

    Checks:
    - MD tag format is valid
    - Reference bases in MD match reference bases in CIGAR
    - Deletions in MD are consistent with CIGAR D operations
    - No mismatches in CIGAR = operations

    :param md_tag:    MD tag string e.g. '5A3^ACC2'
    :param cigar_str: CIGAR string e.g. '5=1X3=3D2='
    :return:          MDValidationResult
    """
    errors = []

    # ── parse CIGAR ───────────────────────────────────────────────────
    cigartuples = _parse_cigar(cigar_str)
    if not cigartuples:
        errors.append(f'Invalid CIGAR string: {cigar_str}')
        return MDValidationResult(False, errors, 0, 0, 0, 0, 0)

    # ── parse MD tag ──────────────────────────────────────────────────
    md_tokens = _parse_md_tag(md_tag)
    if not md_tokens:
        errors.append(f'Invalid or empty MD tag: {md_tag}')
        return MDValidationResult(False, errors,
                                  _cigar_ref_len(cigartuples),
                                  0, 0, 0, 0)

    # ── compute lengths ───────────────────────────────────────────────
    cigar_ref  = _cigar_ref_len(cigartuples)
    md_ref, matches, mismatches, deletions = _md_ref_len(md_tokens)

    # ── check 1 — ref length match ────────────────────────────────────
    if cigar_ref != md_ref:
        errors.append(
            f'Reference length mismatch — '
            f'CIGAR={cigar_ref}  MD={md_ref}'
        )

    # ── check 2 — deletion consistency ───────────────────────────────
    cigar_del_bases = sum(
        length for op, length in cigartuples if op == 2
    )
    if cigar_del_bases != deletions:
        errors.append(
            f'Deletion length mismatch — '
            f'CIGAR deletions={cigar_del_bases}  '
            f'MD deletions={deletions}'
        )

    # ── check 3 — no mismatches in = operations ───────────────────────
    cigar_match_bases = sum(
        length for op, length in cigartuples if op == 7   # = ops
    )
    if cigar_match_bases > 0 and mismatches > 0:
        # count mismatches that fall within = regions
        cigar_mismatch_bases = sum(
            length for op, length in cigartuples if op == 8  # X ops
        )
        if mismatches != cigar_mismatch_bases:
            errors.append(
                f'Mismatch count inconsistency — '
                f'CIGAR X ops={cigar_mismatch_bases}  '
                f'MD mismatches={mismatches}'
            )

    # ── check 4 — MD tag format validity ─────────────────────────────
    if not re.fullmatch(r'[\d\^A-Za-z]+', md_tag):
        errors.append(f'MD tag contains invalid characters: {md_tag}')

    # ── check 5 — MD must not start with mismatch or deletion ─────────
    # if md_tokens and md_tokens[0][0] != 'match':
    #     errors.append(
    #         f'MD tag should start with match count — '
    #         f'got {md_tokens[0][0]}'
    #     )

    return MDValidationResult(
        is_valid      = len(errors) == 0,
        errors        = errors,
        cigar_ref_len = cigar_ref,
        md_ref_len    = md_ref,
        md_matches    = matches,
        md_mismatches = mismatches,
        md_deletions  = deletions
    )

CIGAR_OPS = {
    'M': 0, 'I': 1, 'D': 2, 'N': 3,
    'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8
}

REF_CONSUMING   = {0, 2, 3, 7, 8}
QUERY_CONSUMING = {0, 1, 4, 7, 8}


@dataclass
class CSCigarValidationResult:
    is_valid:        bool
    errors:          list[str]
    cigar_ref_len:   int
    cigar_query_len: int
    cs_ref_len:      int
    cs_query_len:    int
    cs_matches:      int
    cs_mismatches:   int
    cs_insertions:   int
    cs_deletions:    int

    def __str__(self) -> str:
        status = '✅ VALID' if self.is_valid else '❌ INVALID'
        lines  = [
            f'{status}',
            f'  CIGAR  ref={self.cigar_ref_len}   query={self.cigar_query_len}',
            f'  CS     ref={self.cs_ref_len}       query={self.cs_query_len}',
            f'  CS ops  matches={self.cs_matches}  '
            f'mismatches={self.cs_mismatches}  '
            f'ins={self.cs_insertions}  '
            f'del={self.cs_deletions}',
        ]
        if self.errors:
            lines.append('  Errors:')
            for err in self.errors:
                lines.append(f'    ✗ {err}')
        return '\n'.join(lines)


def _parse_cigar(cigar_str: str) -> list[tuple[int, int]]:
    """Parse CIGAR string into (op_int, length) tuples."""
    return [
        (CIGAR_OPS[op], int(length))
        for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_str)
    ]


def _parse_cs_tag(cs_tag: str) -> list[tuple[str, str | int]]:
    """
    Parse CS tag into (op, value) tokens.
    Handles both short and long format:
        Short: =ATCG  *ag  +GG  -ACC
        Long:  :4     *ag  +GG  -ACC
    """
    tokens = []
    for match in re.finditer(r'([=\*\+\-])([A-Za-z]+)|:([0-9]+)', cs_tag):
        if match.group(3) is not None:
            tokens.append((':',  int(match.group(3))))  # long form match
        else:
            tokens.append((match.group(1), match.group(2)))
    return tokens


def _cigar_lengths(cigartuples: list[tuple]) -> tuple[int, int, int, int]:
    """
    Compute ref and query lengths from CIGAR.
    Returns (ref_len, query_len, del_bases, ins_bases)
    """
    ref_len   = sum(l for op, l in cigartuples if op in REF_CONSUMING)
    query_len = sum(l for op, l in cigartuples if op in QUERY_CONSUMING)
    del_bases = sum(l for op, l in cigartuples if op == 2)
    ins_bases = sum(l for op, l in cigartuples if op == 1)
    return ref_len, query_len, del_bases, ins_bases


def _cs_lengths(cs_ops: list[tuple]) -> tuple[int, int, int, int, int, int]:
    """
    Compute ref/query lengths and op counts from CS tag.
    Handles both short (=ATCG) and long (:4) match formats.
    """
    ref_len = query_len = matches = mismatches = insertions = deletions = 0

    for op, val in cs_ops:
        if op == ':':                    # long form match count
            n           = int(val)
            ref_len    += n
            query_len  += n
            matches    += n
        elif op == '=':                  # short form match
            n           = len(val)
            ref_len    += n
            query_len  += n
            matches    += n
        elif op == '*':                  # mismatch
            ref_len    += 1
            query_len  += 1
            mismatches += 1
        elif op == '+':                  # insertion
            query_len  += len(val)
            insertions += len(val)
        elif op == '-':                  # deletion
            ref_len    += len(val)
            deletions  += len(val)

    return ref_len, query_len, matches, mismatches, insertions, deletions


def validate_cs_cigar(cs_tag:    str,
                       cigar_str: str) -> CSCigarValidationResult:
    """
    Validate CS tag consistency against CIGAR string.

    Checks:
    - CS tag format is valid
    - Reference lengths match
    - Query lengths match (excluding soft clips)
    - Insertion counts match
    - Deletion counts match
    - Mismatch tokens are exactly 2 bases
    - CS mismatches consistent with CIGAR X ops

    :param cs_tag:    CS tag string e.g. '=ATC*ag+GG-ACC=AT'
    :param cigar_str: CIGAR string e.g. '3=1X2I3D2='
    :return:          CSCigarValidationResult
    """
    errors = []

    # ── parse CIGAR ───────────────────────────────────────────────────
    cigartuples = _parse_cigar(cigar_str)
    if not cigartuples:
        errors.append(f'Invalid CIGAR string: {cigar_str}')
        return CSCigarValidationResult(False, errors, 0, 0, 0, 0, 0, 0, 0, 0)

    # ── parse CS tag ──────────────────────────────────────────────────
    cs_ops = _parse_cs_tag(cs_tag)
    if not cs_ops:
        errors.append(f'Invalid or empty CS tag: {cs_tag}')
        cigar_ref, cigar_query, _, _ = _cigar_lengths(cigartuples)
        return CSCigarValidationResult(False, errors,
                                       cigar_ref, cigar_query,
                                       0, 0, 0, 0, 0, 0)

    # ── compute lengths ───────────────────────────────────────────────
    cigar_ref, cigar_query, cigar_del, cigar_ins = _cigar_lengths(cigartuples)
    cs_ref, cs_query, matches, mismatches, insertions, deletions = \
        _cs_lengths(cs_ops)

    # soft clips consume query but are absent from CS tag
    soft_clip_bases     = sum(l for op, l in cigartuples if op == 4)
    adj_cigar_query_len = cigar_query - soft_clip_bases

    # ── check 1 — ref length ──────────────────────────────────────────
    if cigar_ref != cs_ref:
        errors.append(
            f'Reference length mismatch — '
            f'CIGAR={cigar_ref}  CS={cs_ref}'
        )

    # ── check 2 — query length ────────────────────────────────────────
    if adj_cigar_query_len != cs_query:
        errors.append(
            f'Query length mismatch — '
            f'CIGAR={adj_cigar_query_len} (excl. soft clips)  '
            f'CS={cs_query}'
        )

    # ── check 3 — insertion count ─────────────────────────────────────
    if cigar_ins != insertions:
        errors.append(
            f'Insertion length mismatch — '
            f'CIGAR={cigar_ins}  CS={insertions}'
        )

    # ── check 4 — deletion count ──────────────────────────────────────
    if cigar_del != deletions:
        errors.append(
            f'Deletion length mismatch — '
            f'CIGAR={cigar_del}  CS={deletions}'
        )

    # ── check 5 — mismatch token format ──────────────────────────────
    for op, val in cs_ops:
        if op == '*' and len(val) != 2:
            errors.append(
                f'Mismatch token *{val!r} must be exactly 2 bases '
                f'(ref+query) — got {len(val)}'
            )

    # ── check 6 — X ops vs CS mismatches ─────────────────────────────
    cigar_x_bases = sum(l for op, l in cigartuples if op == 8)
    if cigar_x_bases > 0 and cigar_x_bases != mismatches:
        errors.append(
            f'Mismatch count inconsistency — '
            f'CIGAR X ops={cigar_x_bases}  '
            f'CS mismatches={mismatches}'
        )

    # ── check 7 — M ops should not coexist with = or X in CS ─────────
    cigar_m_bases = sum(l for op, l in cigartuples if op == 0)
    has_eq_or_x   = any(op in ('=', '*') for op, _ in cs_ops)
    if cigar_m_bases > 0 and has_eq_or_x:
        errors.append(
            f'CIGAR uses M ops but CS uses = and * — '
            f'CIGAR should use = and X for CS consistency'
        )

    # ── check 8 — CS tag format validity ─────────────────────────────
    if not re.fullmatch(r'(?:[=\*\+\-][A-Za-z]+|:[0-9]+)+', cs_tag):
        errors.append(f'CS tag contains invalid characters: {cs_tag}')

    return CSCigarValidationResult(
        is_valid        = len(errors) == 0,
        errors          = errors,
        cigar_ref_len   = cigar_ref,
        cigar_query_len = cigar_query,
        cs_ref_len      = cs_ref,
        cs_query_len    = cs_query,
        cs_matches      = matches,
        cs_mismatches   = mismatches,
        cs_insertions   = insertions,
        cs_deletions    = deletions
    )