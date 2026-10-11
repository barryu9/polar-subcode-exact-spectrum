"""One-command, zero-start reproduction of every code and metric in Table I."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
from statistics import mean
from time import perf_counter

from deterministic_design import construct, search_exclusions
from polar_subcodes import (canonical_checks, complete_directory, encode, extend,
                           parity, polar_transform, right_data, validate_code)

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def expression_list(code):
    expressions = [0] * code["N"]
    for row, sources in code["frozen_expressions"].items():
        expressions[int(row)] = sum(1 << i for i in sources)
    return expressions


def count_enumeration(stage):
    """Recompute S_Q in its original counting scope, including ALL kernel lifts."""
    counts = Counter()
    n, information = stage["N"], stage["information"]
    if stage["kind"] == "V":
        mother = sorted(set(information) | set(stage["rows"]))
        words = complete_directory(n, canonical_checks(n, mother, [0] * n), stage["cutoff"],
                                   counts=counts, shortcuts=False)
    else:
        left_words, seed = stage["count_fixture"]
        released = {j for j in range(n // 2, n) if j not in information
                    and j.bit_count() >= stage["threshold"]}
        # Count the domain-induced COMMON cover, not just the zero-T member.
        checks = tuple(row for row in canonical_checks(n, information, seed)
                       if row.bit_length() - 1 not in released)
        offsets, columns = right_data(n // 2, [row for row in checks if row.bit_length() > n // 2])
        words = []
        for u in left_words:
            words.extend(extend(n // 2, u, offsets, columns, stage["cutoff"], counts))
    assert len(words) == len(set(words))
    counts["final_directory_outputs"] = len(words)
    counts["output_sha256"] = hashlib.sha256("\n".join(map(str, sorted(words))).encode()).hexdigest()
    counts["SQ_exact_fraction"] = str(Fraction(counts["direct_candidates"], counts["quotient_lifts"]))
    return dict(counts)


def one_replay(profile, archived_code):
    """Directories -> exact exclusions -> solve -> recover -> verify -> witness."""
    n, information, target = profile["N"], profile["information_set"], profile["d_target"]
    threshold = min(i.bit_count() for i in information)
    stages = []
    tick = perf_counter()
    expressions = construct(n, information, threshold, target, stages)
    generators = validate_code(n, information, expressions)
    witnesses = []
    for k, word in enumerate(generators):
        if word.bit_count() == target:
            u, _ = encode(1 << k, n, information, expressions)
            witnesses.append({"source": information[k], "input_hex": hex(u),
                              "word_hex": hex(word), "weight": target})
    assert witnesses, "A target-weight witness is needed for equality, not just a lower bound"
    elapsed = perf_counter() - tick
    # This comparison happens AFTER zero-start construction; the saved design
    # is an oracle for regression, never a seed or a substitute for solving.
    assert expressions == expression_list(archived_code), "Design differs from the archived experiment"
    witness = archived_code["target_weight_witness"]
    u, word = int(witness["input_hex"], 16), int(witness["word_hex"], 16)
    assert polar_transform(u, n) == word and word.bit_count() == target
    assert all(not parity(row & u) for row in canonical_checks(n, information, expressions))
    return expressions, stages, elapsed, witnesses[0]


def reproduce(case, repeats, expected, output):
    profile = read(ROOT / "profiles" / f"{case}.json")
    archived_code = read(ROOT / "codes" / f"{case}.json")
    totals, final_stages = [], None
    for run in range(1, repeats + 1):
        print(f"{case}: zero-start replay {run}/{repeats} ...", flush=True)
        expressions, stages, elapsed, witness = one_replay(profile, archived_code)
        public_stages = [{k: v for k, v in stage.items() if k not in ("solver_fixture", "count_fixture")}
                         for stage in stages]
        for actual, saved in zip(public_stages, expected["stages"]):
            for key, value in saved.items():
                assert actual[key] == value, (case, key, actual[key], value)
        assert len(public_stages) == len(expected["stages"])
        totals.append(elapsed)
        final_stages = stages
        write(output / case / f"replay_{run}.json", {"N": profile["N"], "K": profile["K"],
              "frozen_expressions": {str(j): [i for i in profile["information_set"] if v >> i & 1]
                                     for j, v in enumerate(expressions) if v},
              "stages": public_stages, "exact_distance": profile["d_target"],
              "target_weight_witness": witness, "teaching_total_seconds": elapsed})
    last = final_stages[-1]
    audit = count_enumeration(last)
    for key, value in expected["enumeration"].items():
        assert audit[key] == value, (case, key, audit[key], value)
    prefixes = None
    if last["kind"] == "T":
        conditions, width = last["solver_fixture"]
        zero = search_exclusions(conditions, width, ordering="zero")
        assert zero["status"] == "sat" and zero["nodes"] == expected["zero_first_nodes"]
        prefixes = {"zero_first": zero["nodes"], "conditional_mean": last["nodes"]}
    archived_totals = expected["original_total_seconds"]
    assert len(archived_totals) == 5
    assert abs(mean(archived_totals) - expected["original_mean_total_seconds"]) < 1e-12
    result = {"case": case, "distance_verified": profile["d_target"], "stages": public_stages,
              "enumeration": audit, "visited_prefixes_last_stage": prefixes,
              "original_five_replay_seconds": archived_totals,
              "paper_mean_total_seconds": mean(archived_totals),
              "fresh_teaching_replay_seconds": totals, "fresh_teaching_mean_seconds": mean(totals)}
    write(output / case / "summary.json", result)
    print(f"{case}: PASS; d_min={profile['d_target']}; S_Q={audit['SQ_exact_fraction']}; "
          f"fresh mean={mean(totals):.3f}s; archived paper mean={mean(archived_totals):.3f}s", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", help="One case, e.g. n256_k64; default: all six cases")
    parser.add_argument("--replays", type=int, default=5, help="Fresh full constructions per code (default: 5)")
    parser.add_argument("--output", type=Path, default=ROOT / "results", help="Output directory")
    args = parser.parse_args()
    if args.replays < 1:
        parser.error("--replays must be positive")
    evidence = read(ROOT / "data" / "published_table_i.json")
    cases = {row["case"]: row for row in evidence["rows"]}
    if args.case and args.case not in cases:
        parser.error("Unknown case; choose " + ", ".join(cases))
    selected = [args.case] if args.case else list(cases)
    results = [reproduce(case, args.replays, cases[case], args.output) for case in selected]
    write(args.output / "table_i.json", {"python": platform.python_version(), "platform": platform.platform(),
          "fresh_replays": args.replays, "timing_note": "Pure-Python teaching times differ from archived paper times.",
          "rows": results})
    print(f"All {len(selected)} selected cases passed. Results: {args.output.resolve()}", flush=True)


if __name__ == "__main__":
    main()
