"""Paper Sections III--IV: exact exclusions, shared compression and design.

Start reading at construct(), then follow either solve_direct() (V) or
solve_cross_half() (T). No published reference expression is used as a seed.
"""
from fractions import Fraction
from functools import lru_cache
from math import prod

from polar_subcodes import (canonical_checks, complete_directory,
                           polar_transform, quotient_exclusion, right_data,
                           set_bits, solve_binary, validate_code)


class Inconclusive(RuntimeError):
    """A search budget has ended; this is NOT an infeasibility certificate."""


def suffix_coordinates(columns):
    """Rightmost independent columns give coordinates for each suffix image.

    A suffix of unassigned variables spans coordinate bits 0,...,rank-1.
    This makes the image-containment test and exact conditional mean explicit.
    """
    basis, indices = {}, []
    for index in reversed(range(len(columns))):
        value, coordinates = columns[index], 0
        while value:
            p = value.bit_length() - 1
            if p not in basis:
                coordinates ^= 1 << len(indices)
                basis[p] = value, coordinates
                indices.append(index)
                break
            row, label = basis[p]
            value, coordinates = value ^ row, coordinates ^ label
    def encode_label(value):
        result = 0
        while value:
            p = value.bit_length() - 1
            if p not in basis:
                raise ValueError("Label is outside the observation image")
            row, label = basis[p]
            value, result = value ^ row, result ^ label
        return result
    coordinates = tuple(encode_label(v) for v in columns)
    return coordinates, {index: k for k, index in enumerate(indices)}, encode_label


def prepare_conditions(maps, forbidden, offsets=None):
    """Affine exclusions c_g + A_g gamma NOT IN B_g (Section IV)."""
    if offsets is None:
        offsets = [0] * len(maps)
    conditions = []
    for columns, bad, c in zip(maps, forbidden, offsets):
        coordinates, pivots, encode_label = suffix_coordinates(columns)
        targets = set()
        for label in bad:
            try:
                targets.add(encode_label(label ^ c))
            except ValueError:
                pass  # This forbidden label can never occur in the allowed image.
        conditions.append((coordinates, pivots, tuple(sorted(targets))))
    return conditions


def conditional_mean(conditions, prefix=0, length=0):
    """Exact Phi(pi): expected number of violated groups under uniform completion."""
    total = Fraction(0)
    for columns, pivots, targets in conditions:
        offset = 0
        for p in range(length):
            if prefix >> p & 1:
                offset ^= columns[p]
        remaining_rank = sum(p >= length for p in pivots)
        reachable = sum((label ^ offset) >> remaining_rank == 0 for label in targets)
        total += Fraction(reachable, 1 << remaining_rank)
    return total


def search_exclusions(conditions, width, *, ordering="mean", node_limit=None):
    """Algorithm 2: finite DFS; ONLY image containment prunes a branch.

    Both children are retained. Phi >= 1 never prunes; Phi < 1 is a sufficient
    completion guarantee (Lemma 1). Return unsat only after full exhaustion.
    """
    if ordering not in ("mean", "zero"):
        raise ValueError("Unknown branch ordering")
    max_rank = max((len(pivots) for _, pivots, _ in conditions), default=0)
    root = [(k, 0, len(pivots), targets) for k, (_, pivots, targets) in enumerate(conditions) if targets]
    def potential(state):
        return sum(len(bad) << (max_rank - rank) for _, _, rank, bad in state)
    initial_mean = Fraction(potential(root), 1 << max_rank)
    frontier, nodes = [(0, 0, root)], 0
    while frontier:
        if node_limit is not None and nodes >= node_limit:
            raise Inconclusive("Node budget ended with unexplored branches")
        p, value, state = frontier.pop()
        nodes += 1
        # All labels of a residual image are forbidden: sound containment prune.
        if any(len(bad) == (1 << rank) for _, _, rank, bad in state):
            continue
        if p == width:
            assert not state
            return {"status": "sat", "value": value, "nodes": nodes, "initial_mean": str(initial_mean)}
        children = []
        for bit in (0, 1):
            child = []
            for k, offset, rank, bad in state:
                columns, pivots, _ = conditions[k]
                next_offset = offset ^ (columns[p] if bit else 0)
                next_rank = rank - int(p in pivots)
                reachable = tuple(t for t in bad if (t ^ next_offset) >> next_rank == 0)
                if reachable:
                    child.append((k, next_offset, next_rank, reachable))
            children.append((bit, child))
        if ordering == "mean":
            children.sort(key=lambda item: (potential(item[1]), item[0]))
        for bit, child in reversed(children):
            frontier.append((p + 1, value | (bit << p), child))
    return {"status": "unsat", "nodes": nodes, "initial_mean": str(initial_mean)}


def joint_observation(maps, width):
    """Theorem 3: rank(L_D), representative columns and ker(L_D).

    Columns at row-echelon pivot indices span the complete observation image.
    Setting other coefficients to zero chooses a representative, not a smaller
    set of attainable spectra. The kernel is the uniform spectral freedom.
    """
    basis = {}
    for columns in maps:
        for b in range(max((v.bit_length() for v in columns), default=0)):
            row = sum(((v >> b) & 1) << p for p, v in enumerate(columns))
            while row:
                pivot = row.bit_length() - 1
                if pivot not in basis:
                    basis[pivot] = row
                    break
                row ^= basis[pivot]
    return sorted(basis), list(basis.values())


def uniform_spectral_kernel(maps, width):
    """Theorem 3: every returned direction preserves all observed cutoff spectra.

    Constituents must stay fixed. This is uniform additive equivalence over
    the entire outer coefficient space, not a statement about SCL decoding.
    """
    _, rows = joint_observation(maps, width)
    return solve_binary([(row, 0) for row in rows], width)[1]


class SourcePatterns:
    """Theorem 1 exclusions grouped by source scope for the small V stages.

    A variable is the pattern of one source across all released frozen rows.
    This is the experiment's exact grouped zero-first/sparse/product ordering;
    the binary conditional-mean ordering is implemented separately above.
    """
    def __init__(self, information, rows, words):
        self.width = len(rows)
        self.sources = sorted((i for i in information if any(i < j for j in rows)),
                              key=lambda i: (-i.bit_count(), -i))
        domains = []
        for source in self.sources:
            allowed = sum((source < j) << q for q, j in enumerate(rows))
            domains.append(sum(1 << v for v in range(1 << self.width) if not v & ~allowed))
        exclusions, groups = {}, set()
        for u in words:
            scope = sum(((u >> source) & 1) << k for k, source in enumerate(self.sources))
            defect = sum(((u >> j) & 1) << q for q, j in enumerate(rows))
            groups.add((scope, defect))
            exclusions[scope] = exclusions.get(scope, 0) | (1 << defect)
        self.groups = len(groups)
        self.root_conflict = bool(exclusions.get(0, 0) & 1)
        exclusions.pop(0, None)
        for scope, forbidden in exclusions.items():
            if not scope & (scope - 1):
                k = scope.bit_length() - 1
                domains[k] &= ~forbidden
        self.root_conflict |= any(not domain for domain in domains)
        self.initial = ([-1] * len(self.sources), domains, exclusions)
        full = (1 << (1 << self.width)) - 1
        self.masks = [full // ((1 << (1 << q)) + 1) for q in range(self.width)]

    @lru_cache(maxsize=32768)
    def translate(self, forbidden, value):
        """Bitset XOR permutation: new forbidden labels are old labels XOR value."""
        while value:
            bit = value & -value
            mask = self.masks[bit.bit_length() - 1]
            forbidden = ((forbidden & mask) << bit) | ((forbidden >> bit) & mask)
            value ^= bit
        return forbidden

    def assign(self, state, k, value):
        old_values, old_domains, old_exclusions = state
        values, domains = old_values.copy(), old_domains.copy()
        values[k], domains[k] = value, 1 << value
        exclusions = {}
        for scope, forbidden in old_exclusions.items():
            if scope >> k & 1:
                scope ^= 1 << k
                forbidden = self.translate(forbidden, value)
            if not scope:
                if forbidden & 1:
                    return None
                continue
            # Union identical residual scopes; this is exact conjunction.
            exclusions[scope] = exclusions.get(scope, 0) | forbidden
        for scope, forbidden in exclusions.items():
            if not scope & (scope - 1):
                j = scope.bit_length() - 1
                domains[j] &= ~forbidden
                if not domains[j]:
                    return None
        return values, domains, exclusions

    def frame(self, state):
        values, domains, _ = state
        k = values.index(-1)
        options = sorted(set_bits(domains[k]), key=lambda v: (v.bit_count(), v))
        if not (domains[k] & 1 and self.assign(state, k, 0) is not None):
            scored = []
            for v in options:
                child = self.assign(state, k, v)
                if child is not None:
                    sizes = [domain.bit_count() for j, domain in enumerate(child[1]) if child[0][j] == -1]
                    scored.append((-prod(sizes), v.bit_count(), v))
            options = [v for _, _, v in sorted(scored)]
        return [state, k, options, 0]

    def solve(self, node_limit=None):
        if self.root_conflict:
            return None, 0
        if not self.sources:
            return [], 0
        frontier, nodes = [self.frame(self.initial)], 0
        while frontier:
            if node_limit is not None and nodes >= node_limit:
                raise Inconclusive("Grouped search budget ended with unexplored branches")
            state, k, options, cursor = frontier[-1]
            if cursor == len(options):
                frontier.pop()
                continue
            frontier[-1][3] += 1
            nodes += 1
            child = self.assign(state, k, options[cursor])
            if child is None:
                continue
            if -1 not in child[0]:
                return child[0], nodes
            frontier.append(self.frame(child))
        return None, nodes


def exclusion_cutoff(n, threshold, target):
    """Optional RM mother: use its gap only in the stated admissible range."""
    dp = 1 << threshold
    order = n.bit_length() - 1 - threshold
    gap = 2 <= order <= n.bit_length() - 3 and dp < target <= 3 * dp // 2
    return (dp if gap else target - 1), gap


def solve_direct(n, information, rows, cutoff):
    """Release selected checks, enumerate the common supercode, then exclude."""
    mother = sorted(set(information) | set(rows))
    words = complete_directory(n, canonical_checks(n, mother, [0] * n), cutoff)
    # A surviving canonical input must begin at an information index.
    words = [u for u in words if u and (u & -u).bit_length() - 1 in information]
    problem = SourcePatterns(information, rows, words)
    values, nodes = problem.solve()
    if values is None:
        raise RuntimeError("The complete declared V domain is infeasible")
    expressions = [0] * n
    for source, value in zip(problem.sources, values):
        for q in set_bits(value):
            expressions[rows[q]] |= 1 << source
    problem.translate.cache_clear()
    return expressions, {"q": sum(i < j for j in rows for i in information),
                         "groups": problem.groups, "nodes": nodes}


def solve_cross_half(n, information, seed, threshold, cutoff):
    """Fix both constituents; solve the shared cross-half T observation system."""
    half = n // 2
    sources = [i for i in information if i < half]
    right_info = [i - half for i in information if i >= half]
    right_checks = canonical_checks(half, right_info, [v >> half for v in seed[half:]])
    _, columns = right_data(half, [row << half for row in right_checks])
    frozen = [j for j in range(half) if j not in right_info]
    rows = [half + j for j in frozen if (half + j).bit_count() >= threshold]
    selected = [frozen.index(j - half) for j in rows]
    words = complete_directory(half, canonical_checks(half, sources, seed[:half]), cutoff)
    words = [u for u in words if u]
    pure_right = complete_directory(half, right_checks, cutoff // 2)
    assert pure_right == [0], "The fixed right constituent must satisfy the pure-right bound"
    maps = []
    for u in words:
        # Table I's optional RM gap gives radius zero: B_a(D) = {0}.
        assert polar_transform(u, half).bit_count() == cutoff
        images, bad, _ = quotient_exclusion(u, half, columns, len(frozen), cutoff)
        assert bad == {0}
        row_images = [images[k] for k in selected]
        maps.append(tuple(image if u >> source & 1 else 0 for source in sources for image in row_images))
    q = len(sources) * len(rows)
    positions, joint_rows = joint_observation(maps, q)
    conditions = prepare_conditions([[m[p] for p in positions] for m in maps], [{0}] * len(maps))
    answer = search_exclusions(conditions, len(positions))
    if answer["status"] != "sat":
        raise RuntimeError("The complete declared T domain is infeasible")
    expressions = seed.copy()
    for b, p in enumerate(positions):
        if answer["value"] >> b & 1:
            source, row = divmod(p, len(rows))
            expressions[rows[row]] ^= 1 << sources[source]
    # Re-evaluate every original exclusion, independently of DFS propagation.
    for m in maps:
        label = 0
        for b, p in enumerate(positions):
            if answer["value"] >> b & 1:
                label ^= m[p]
        assert label != 0
    return expressions, {"q": q, "r": len(positions), "groups": len(words),
                         "nodes": answer["nodes"], "Phi": answer["initial_mean"],
                         "kernel_dimension": q - len(positions)}, (conditions, len(positions)), (words, seed)


def construct(n, information, threshold, target, stages):
    """Section V's zero-start recursion; the short tutorial's main design flow."""
    dp = min((1 << i.bit_count() for i in information), default=n + 1)
    if dp >= target:
        return [0] * n  # The unprecoded constituent already meets its target.
    rows = [j for j in range(n) if j not in information and j.bit_count() >= threshold
            and any(i < j for i in information)]
    cutoff, gap = exclusion_cutoff(n, threshold, target)
    if len(rows) <= 12:
        expressions, metrics = solve_direct(n, information, rows, cutoff)
        stage = {"kind": "V", "N": n, "K": len(information), "rows": rows,
                 "information": information, "cutoff": cutoff, "target": target,
                 "RM_gap_used": gap, **metrics}
    else:
        half = n // 2
        right_info = [i - half for i in information if i >= half]
        right = construct(half, right_info, threshold - 1, (target + 1) // 2, stages)
        seed = [0] * half + [v << half for v in right]
        left_dp = min((1 << i.bit_count() for i in information if i < half), default=n + 1)
        if left_dp >= target:
            # Static Plotkin lift: min(left_dp, 2*right_target) >= target.
            expressions = seed
            stage = None
        else:
            assert gap, "The Table I cross-half tutorial uses radius-zero RM exclusions"
            expressions, metrics, solver_fixture, count_fixture = solve_cross_half(
                n, information, seed, threshold, cutoff)
            stage = {"kind": "T", "N": n, "K": len(information), "cutoff": cutoff,
                     "target": target, "RM_gap_used": gap, "information": information,
                     "threshold": threshold, **metrics,
                     "solver_fixture": solver_fixture, "count_fixture": count_fixture}
    validate_code(n, information, expressions)
    assert all(i.bit_count() >= threshold for i in information)
    assert all(not v or j.bit_count() >= threshold for j, v in enumerate(expressions))
    # Complete lower-layer exclusion plus the optional RM gap proves the bound.
    remaining = complete_directory(n, canonical_checks(n, information, expressions), cutoff)
    assert remaining == [0], "A nonzero word below the target has survived"
    if stage is not None:
        stages.append(stage)
    return expressions
