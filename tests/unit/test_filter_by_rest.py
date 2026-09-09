# -*- coding: utf-8 -*-

"""Tests for filtering by the loci `filter_top_3` holds aside.

The point of these is byte equivalence: `finalize_results` replaces the pass
`change_output_by_extra_gl` used to make over `don.pmug`, and for every donor
whose correct answer was already inside the top `number_of_results` it has to
produce exactly the same four output files as that pass did.
"""

import io

import pytest

from grim.filter_by_rest import (
    aggregate_muug_results_from_haps,
    aggregate_pop_results_from_haps,
    create_haps,
    filter_results,
    finalize_results,
    truncate_haps,
    write_filter,
    write_umug,
    write_umug_pops,
)
from grim.imputation.impute import (
    write_best_hap_race_pairs,
    write_best_prob,
    write_best_prob_genotype,
)

NUMBER_OF_RESULTS = 20
NUMBER_OF_POP_RESULTS = 100

CONFIG = {
    "output_MUUG": True,
    "output_haplotypes": True,
    "number_of_results": NUMBER_OF_RESULTS,
    "number_of_pop_results": NUMBER_OF_POP_RESULTS,
}


def hap(a, b, c, dqb1, drb1):
    return "~".join(
        [
            "A*" + a,
            "B*" + b,
            "C*" + c,
            "DQB1*" + dqb1,
            "DRB1*" + drb1,
        ]
    )


def build_res_haps(consistent_rank):
    """A candidate set where only the pair at `consistent_rank` carries
    DQB1*06:02 on both chromosomes - the phase the worked example needs.

    `consistent_rank` of None means no pair is consistent.
    """
    haps, probs, pops = [], [], []
    for rank in range(40):
        dqb1_2 = "06:02" if rank == consistent_rank else "05:01"
        haps.append(
            [
                hap("01:01", "08:01", "04:01", "06:02", "03:01"),
                hap("33:03", "44:03", "07:01", dqb1_2, "07:01"),
            ]
        )
        pops.append(["CAU", "CAU"])
        probs.append(1.0 / (rank + 1) ** 2)
    return {"MaxProb": probs[0], "Haps": haps, "Probs": probs, "Pops": pops}


# The two genotypes `build_res_haps` can produce, written out rather than
# aggregated: `build_res_muugs` is the fixture for tests of the aggregation, so
# it must not be built by the function under test.
MUUG_MISMATCHED = (
    "A*01:01+A*33:03^B*08:01+B*44:03^C*04:01+C*07:01"
    "^DQB1*05:01+DQB1*06:02^DRB1*03:01+DRB1*07:01"
)
MUUG_CONSISTENT = (
    "A*01:01+A*33:03^B*08:01+B*44:03^C*04:01+C*07:01"
    "^DQB1*06:02+DQB1*06:02^DRB1*03:01+DRB1*07:01"
)


def build_res_muugs(res_haps, consistent_rank=None):
    """What the MUUG pass produces for the same subject.

    Keys go in in order of first appearance, the way the imputation inserts
    them, because the writers break ties on insertion order.
    """
    muugs = {}
    for rank, prob in enumerate(res_haps["Probs"]):
        muug = MUUG_CONSISTENT if rank == consistent_rank else MUUG_MISMATCHED
        muugs[muug] = muugs.get(muug, 0) + prob
    return {
        "MaxProb": res_haps["MaxProb"],
        "Haps": muugs,
        "Pops": {"CAU,CAU": sum(res_haps["Probs"])},
    }


def build_phase_swaps(n_pairs):
    """`n_pairs` different phasings of one genotype, all of them consistent
    with `EXTRA_GL`.

    Every pair aggregates to the same MUUG, so how far down the ranking the
    aggregation reaches is the only thing that can change the probability
    written for it - which is exactly what `umug_from_full_results` decides.
    """
    swappable = [
        ("01:01", "33:03"),
        ("08:01", "44:03"),
        ("04:01", "07:01"),
        ("03:01", "07:01"),
    ]
    haps, probs, pops = [], [], []
    for rank in range(n_pairs):
        picked = [pair[(rank >> bit) & 1] for bit, pair in enumerate(swappable)]
        other = [pair[1 - ((rank >> bit) & 1)] for bit, pair in enumerate(swappable)]
        haps.append(
            [
                hap(picked[0], picked[1], picked[2], "06:02", picked[3]),
                hap(other[0], other[1], other[2], "06:02", other[3]),
            ]
        )
        pops.append(["CAU", "CAU"])
        probs.append(1.0 / (rank + 1) ** 2)
    return {"MaxProb": probs[0], "Haps": haps, "Probs": probs, "Pops": pops}


EXTRA_GL = "C*04:01+C*07:01^DQB1*06:02+DQB1*06:02"


def run_old_path(subject_id, res_muugs, res_haps, extra_gl, tmp_path):
    """Impute writes the four files, then `change_output_by_extra_gl` rewrites
    them off the top `number_of_results` rows of `don.pmug`.
    """
    pmug_path = str(tmp_path / "don.pmug")
    with open(pmug_path, "w") as fout:
        write_best_hap_race_pairs(
            subject_id,
            res_haps["Haps"],
            res_haps["Pops"],
            res_haps["Probs"],
            NUMBER_OF_RESULTS,
            fout,
        )

    # The body of `change_output_by_extra_gl`, for one subject.
    read_back = create_haps(pmug_path)["res_haps"][0]
    if extra_gl:
        read_back = filter_results(read_back, extra_gl)

    files = {name: io.StringIO() for name in ("pmug", "pmug.pops", "umug", "umug.pops")}
    missed = len(read_back["Haps"]) == 0
    if not missed:
        write_filter(
            subject_id,
            read_back,
            files["pmug"],
            files["pmug.pops"],
            files["umug"],
            files["umug.pops"],
            NUMBER_OF_RESULTS,
            NUMBER_OF_POP_RESULTS,
            True,
            True,
        )
    return {name: buf.getvalue() for name, buf in files.items()}, missed


def run_new_path(
    subject_id, res_muugs, res_haps, extra_gl, umug_from_full_results=False
):
    """`finalize_results` in the worker, then the writers in `impute.py`."""
    res_muugs, res_haps, missed = finalize_results(
        res_muugs,
        res_haps,
        extra_gl,
        NUMBER_OF_RESULTS,
        umug_from_full_results,
    )

    files = {name: io.StringIO() for name in ("pmug", "pmug.pops", "umug", "umug.pops")}
    write_best_hap_race_pairs(
        subject_id,
        res_haps["Haps"],
        res_haps["Pops"],
        res_haps["Probs"],
        NUMBER_OF_RESULTS,
        files["pmug"],
    )
    write_best_prob(
        subject_id, res_haps["Pops"], res_haps["Probs"], 1, files["pmug.pops"]
    )
    write_best_prob_genotype(
        subject_id, res_muugs["Haps"], NUMBER_OF_RESULTS, files["umug"]
    )
    write_best_prob_genotype(
        subject_id, res_muugs["Pops"], NUMBER_OF_POP_RESULTS, files["umug.pops"]
    )
    return {name: buf.getvalue() for name, buf in files.items()}, missed


@pytest.mark.parametrize(
    "extra_gl, consistent_rank",
    [
        ("", None),  # nothing held aside - the filter must not touch this donor
        (EXTRA_GL, 0),  # correct answer at the top
        (EXTRA_GL, NUMBER_OF_RESULTS - 1),  # correct answer on the cut
    ],
)
def test_new_path_matches_old_output(tmp_path, extra_gl, consistent_rank):
    res_haps = build_res_haps(consistent_rank)
    old, old_missed = run_old_path(
        "25014660",
        build_res_muugs(res_haps, consistent_rank),
        res_haps,
        extra_gl,
        tmp_path,
    )

    res_haps = build_res_haps(consistent_rank)
    new, new_missed = run_new_path(
        "25014660", build_res_muugs(res_haps, consistent_rank), res_haps, extra_gl
    )

    assert new_missed == old_missed is False
    assert new == old


def test_answer_below_the_cut_is_recovered(tmp_path):
    """The bug itself: a consistent pair ranked past `number_of_results`."""
    below_the_cut = NUMBER_OF_RESULTS + 5

    res_haps = build_res_haps(below_the_cut)
    old, old_missed = run_old_path(
        "25014660",
        build_res_muugs(res_haps, below_the_cut),
        res_haps,
        EXTRA_GL,
        tmp_path,
    )
    assert old_missed
    assert old["pmug"] == ""

    res_haps = build_res_haps(below_the_cut)
    new, new_missed = run_new_path(
        "25014660", build_res_muugs(res_haps, below_the_cut), res_haps, EXTRA_GL
    )
    assert not new_missed
    assert new["pmug"].count("\n") == 1
    assert "DQB1*06:02~DRB1*07:01" in new["pmug"]


def test_no_consistent_pair_anywhere_is_still_a_miss(tmp_path):
    res_haps = build_res_haps(None)
    _, old_missed = run_old_path(
        "25014660", build_res_muugs(res_haps), res_haps, EXTRA_GL, tmp_path
    )
    assert old_missed

    res_haps = build_res_haps(None)
    new, new_missed = run_new_path(
        "25014660", build_res_muugs(res_haps), res_haps, EXTRA_GL
    )
    assert new_missed
    assert new == {"pmug": "", "pmug.pops": "", "umug": "", "umug.pops": ""}


def test_aggregations_match_the_writers_they_were_lifted_from():
    """`write_umug`/`write_umug_pops` aggregate inline; the new path aggregates
    first and writes with `write_best_prob_genotype`. Same bytes either way.
    """
    res_haps = truncate_haps(build_res_haps(3), NUMBER_OF_RESULTS)

    old_umug, old_pops = io.StringIO(), io.StringIO()
    write_umug("25014660", res_haps, old_umug, NUMBER_OF_RESULTS)
    write_umug_pops("25014660", res_haps, old_pops, NUMBER_OF_POP_RESULTS)

    new_umug, new_pops = io.StringIO(), io.StringIO()
    write_best_prob_genotype(
        "25014660",
        aggregate_muug_results_from_haps(res_haps),
        NUMBER_OF_RESULTS,
        new_umug,
    )
    write_best_prob_genotype(
        "25014660",
        aggregate_pop_results_from_haps(res_haps),
        NUMBER_OF_POP_RESULTS,
        new_pops,
    )

    assert new_umug.getvalue() == old_umug.getvalue()
    assert new_pops.getvalue() == old_pops.getvalue()


def test_truncate_sorts_even_when_there_is_nothing_to_cut():
    res_haps = build_res_haps(None)
    res_haps["Haps"].reverse()
    res_haps["Probs"].reverse()
    res_haps["Pops"].reverse()

    returned = truncate_haps(res_haps, 1000)

    assert returned is res_haps
    assert res_haps["Probs"] == sorted(res_haps["Probs"], reverse=True)
    assert len(res_haps["Haps"]) == 40


def test_which_of_the_two_filters_rewrites_its_argument():
    """`truncate_haps` rewrites the dict it is handed and returns that same
    dict, so a caller may use either. `filter_results` beside it does not: when
    nothing survives it returns a *different* dict, which is why
    `finalize_results` assigns the return value of both.
    """
    res_haps = build_res_haps(0)
    assert truncate_haps(res_haps, 5) is res_haps

    survives = build_res_haps(0)
    assert filter_results(survives, EXTRA_GL) is survives

    none_survive = build_res_haps(None)
    assert filter_results(none_survive, EXTRA_GL) is not none_survive


def test_umug_from_full_results_reaches_past_the_cut(tmp_path):
    """The second switch: where `don.umug` gets its probabilities from.

    Every pair here is a phasing of the same genotype and consistent with
    `EXTRA_GL`, so the MUUG is identical in both modes and only the probability
    written for it can move - False stops at `number_of_results`, True sums the
    whole filtered set.
    """
    n_pairs = NUMBER_OF_RESULTS + 20

    def muugs_for(res_haps):
        total = sum(res_haps["Probs"])
        return {
            "MaxProb": res_haps["MaxProb"],
            "Haps": {MUUG_CONSISTENT: total},
            "Pops": {"CAU,CAU": total},
        }

    res_haps = build_phase_swaps(n_pairs)
    probs = list(res_haps["Probs"])
    old, _ = run_old_path("25014660", muugs_for(res_haps), res_haps, EXTRA_GL, tmp_path)

    res_haps = build_phase_swaps(n_pairs)
    top_n, _ = run_new_path("25014660", muugs_for(res_haps), res_haps, EXTRA_GL)

    res_haps = build_phase_swaps(n_pairs)
    full, _ = run_new_path(
        "25014660",
        muugs_for(res_haps),
        res_haps,
        EXTRA_GL,
        umug_from_full_results=True,
    )

    # Off, the default: `don.umug` is what the old post-imputation pass wrote.
    assert top_n == old
    assert top_n["umug"] == "25014660,{muug},{prob},0\n".format(
        muug=MUUG_CONSISTENT, prob=sum(probs[:NUMBER_OF_RESULTS])
    )

    # On: the phased files are untouched, the MUUG files carry the full sum.
    assert full["pmug"] == old["pmug"]
    assert full["pmug.pops"] == old["pmug.pops"]
    assert full["umug"] == "25014660,{muug},{prob},0\n".format(
        muug=MUUG_CONSISTENT, prob=sum(probs)
    )
    assert full["umug.pops"] == "25014660,CAU,CAU,{prob},0\n".format(prob=sum(probs))
    assert full["umug"] != top_n["umug"]


def test_results_without_haps_are_left_alone():
    """`haps_output` off, or an imputation that produced no phased pairs: there
    is nothing to filter or aggregate, so nothing may be overwritten.
    """
    res_muugs = {
        "MaxProb": 0.5,
        "Haps": {"A*01:01+A*01:01": 0.5},
        "Pops": {"CAU,CAU": 0.5},
    }
    res_haps = {"Haps": "Nan", "Probs": 0, "Pops": {}}

    out_muugs, out_haps, missed = finalize_results(
        res_muugs, res_haps, EXTRA_GL, NUMBER_OF_RESULTS
    )

    assert not missed
    assert out_muugs["Haps"] == {"A*01:01+A*01:01": 0.5}
    assert out_haps["Haps"] == "Nan"
