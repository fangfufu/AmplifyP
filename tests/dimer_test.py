from amplifyp.dimer import analyze_group, calculate_dimer
from amplifyp.dna import Primer
from tests.examples.amplify4_examples import (
    primer_11bp,
    primer_1701,
    primer_10289,
    primer_10290,
)


def test_calculate_dimer_strong() -> None:
    # AAAA vs TTTT -> 4 A-T pairs = 4 * 20 = 80. Serious > 60.
    p1 = Primer("AAAA")
    p2 = Primer("TTTT")

    res = calculate_dimer(p1, p2)
    assert res.quality == 80.0
    assert res.overlap == 4
    assert res.serious is True
    # p1 is shorter/equal? Equal. logic uses p2 as long.
    # calculate_dimer logic: if len(p1) < len(p2): p1 short. else: p2 short.
    # So p1=AAAA, p2=TTTT -> len same -> short=TTTT(p2), long=AAAA(p1).
    # Wait, let's check my implementation logic again.
    # if len(p1) < len(p2): short=p1
    # else: short=p2, long=p1.
    # So passing (AAAA, TTTT) -> short=TTTT, long=AAAA.
    # Result object has p1=short, p2=long.
    # So res.p1 should be TTTT, res.p2 should be AAAA.
    assert res.p1.seq == "TTTT"
    assert res.p2.seq == "AAAA"


def test_calculate_dimer_self_palindromic() -> None:
    # 5' GGCC 3'
    # 3' CCGG 5'
    # Perfect self-dimer. 4 G-C pairs = 4 * 30 = 120.
    p = Primer("GGCC")
    res = calculate_dimer(p, p)

    assert res.quality == 120.0
    assert res.overlap == 4
    assert res.serious is True
    assert res.p1 == p
    assert res.p2 == p


def test_calculate_dimer_none() -> None:
    # AAAA vs AAAA.
    # A-A pairing is -20.
    # 4 overlaps -> -80.
    p = Primer("AAAA")
    res = calculate_dimer(p, p)

    assert res.quality < 0
    assert res.serious is False


def test_calculate_dimer_group() -> None:
    p1 = Primer("AAAA", name="P1")
    p2 = Primer("TTTT", name="P2")
    p3 = Primer("CCCC", name="P3")  # Self-dimer? C-C is -20.

    # P1-P2 is serious (80).
    # Others are not.

    results = analyze_group([p1, p2, p3])

    assert len(results) >= 1
    # Check that we found (P1, P2) or (P2, P1)
    # The result stores pointers to the primer objects.
    names = {(r.p1.name, r.p2.name) for r in results}
    # Logic: short=TTTT, long=AAAA. So P2, P1.
    # But analyze_group calls combinations(P1, P2). (A, B).
    # If passed in that order, calc(A, B) -> A is short if len <.
    # Here len equal. So short=B(P2), long=A(P1).
    # So result should be (P2, P1).

    # Wait, is there any other interaction?
    # P3-P3 (CCCC-CCCC) -> C-C -20.
    # P1-P3 (AAAA-CCCC) -> A-C -20.
    # P2-P3 (TTTT-CCCC) -> T-C -20.

    assert ("P2", "P1") in names or ("P1", "P2") in names

    # Ensure score is 80
    serious_results = [r for r in results if r.quality == 80]
    assert len(serious_results) == 1


def test_calculate_dimer_offset_alignment() -> None:
    # P1: 5' AT 3'
    # P2: 5' GGATCC 3'
    # P1 should bind to the AT in the middle?
    # Indexing:
    # P2: G G A T C C
    #     0 1 2 3 4 5
    # P1 (3'->5'): T A
    # If aligned at pos 2 (A):
    # p1(T) vs p2(A) (pos 2)
    # p1(A) vs p2(T) (pos 3)
    # match!

    p1 = Primer("AT")
    p2 = Primer("GGATCC")

    res = calculate_dimer(p1, p2)

    # Expected: 2 overlaps.
    # A-T (20) + T-A (20) = 40.
    # Serious? Cutoff 60. So False.
    assert res.quality == 40.0
    assert res.p1_pos == 2
    assert res.serious is False

    # Let's lower cutoff to verify logic works
    from dataclasses import replace

    from amplifyp.settings import DEFAULT_PRIMER_DIMER_SETTINGS

    custom_settings = replace(
        DEFAULT_PRIMER_DIMER_SETTINGS, threshold=30.0, overlap=2
    )
    res_low = calculate_dimer(p1, p2, settings=custom_settings)
    assert res_low.serious is True


def test_calculate_dimer_real_examples() -> None:
    # Real examples from Amplify4
    result_10289_10289 = calculate_dimer(primer_10289, primer_10289)
    assert result_10289_10289.quality == 160
    assert result_10289_10289.overlap == 10

    result_10290_10290 = calculate_dimer(primer_10290, primer_10290)
    assert result_10290_10290.quality == 200
    assert result_10290_10290.overlap == 12

    result_10290_11bp = calculate_dimer(primer_10290, primer_11bp)
    assert result_10290_11bp.quality == 70
    assert result_10290_11bp.overlap == 6

    result_1701_1701 = calculate_dimer(primer_1701, primer_1701)
    assert result_1701_1701.quality == 120
    assert result_1701_1701.overlap == 8

    result_dimer_group = analyze_group(
        [primer_10289, primer_10290, primer_11bp, primer_1701]
    )
    assert len(result_dimer_group) == 4
