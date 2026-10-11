"""Paper Sections II--III: canonical checks and complete quotient enumeration.

Bit i of an integer is coordinate i. Indices start at zero; there is no bit
reversal. This module uses only the Python standard library.
"""
from collections import Counter
from itertools import combinations
from math import comb


def set_bits(value):
    """Visit the nonzero coordinates of a binary vector, in increasing order."""
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def parity(value):
    return value.bit_count() & 1


def polar_transform(value, n):
    """u -> u F^(tensor log2(n)); the binary polar transform is its own inverse."""
    if n < 1 or n & (n - 1) or value < 0 or value >> n:
        raise ValueError("Expected a binary vector of power-of-two length")
    full, step = (1 << n) - 1, 1
    while step < n:
        # Alternating blocks of 'step' ones and zeros select the left halves.
        left_halves = full // ((1 << step) + 1)
        value ^= (value >> step) & left_halves
        step *= 2
    return value


def encode(message, n, information, expressions):
    """Canonical frozen expressions: u_j = sum_i s_ji u_i over earlier info bits."""
    if message < 0 or message >> len(information):
        raise ValueError("Message exceeds the information-set dimension")
    u = sum(((message >> k) & 1) << i for k, i in enumerate(information))
    for j, sources in enumerate(expressions):
        u |= parity(u & sources) << j
    return u, polar_transform(u, n)


def canonical_checks(n, information, expressions):
    """Section II: h_j = e_j + s_j for each frozen index j."""
    info = set(information)
    return tuple((1 << j) ^ expressions[j] for j in range(n) if j not in info)


def validate_code(n, information, expressions):
    """Check dimension and causality before any distance claim is made."""
    if n < 1 or n & (n - 1) or len(expressions) != n:
        raise ValueError("Invalid polar code length")
    if sorted(set(information)) != list(information) or any(i < 0 or i >= n for i in information):
        raise ValueError("Information indices must be sorted, distinct and in range")
    info_mask = sum(1 << i for i in information)
    for j, sources in enumerate(expressions):
        if sources < 0 or sources & ~info_mask or sources >> j or (j in information and sources):
            raise ValueError(f"Noncanonical or noncausal expression at index {j}")
    # Canonical checks have distinct pivots, so the code dimension is |I|.
    generators = [encode(1 << k, n, information, expressions)[1] for k in range(len(information))]
    if ImageFactor(generators).rank != len(information):
        raise AssertionError("Generator dimension does not match |I|")
    return generators


def gray_span(basis):
    """Enumerate EVERY vector of a span, including zero (no kernel lift is lost)."""
    value = 0
    yield value
    for index in range(1, 1 << len(basis)):
        value ^= basis[(index & -index).bit_length() - 1]
        yield value


def embed(value, positions):
    return sum(((value >> k) & 1) << j for k, j in enumerate(positions))


def solve_binary(equations, width):
    """Solve row*x = rhs over F_2; return one origin and a nullspace basis."""
    rows = {}
    for value, rhs in equations:
        for pivot in sorted(rows):
            if value >> pivot & 1:
                row, bit = rows[pivot]
                value, rhs = value ^ row, rhs ^ bit
        if not value:
            if rhs:
                return None
            continue
        pivot = (value & -value).bit_length() - 1
        for old, (row, bit) in list(rows.items()):
            if row >> pivot & 1:
                rows[old] = row ^ value, bit ^ rhs
        rows[pivot] = value, rhs
    origin = sum(rhs << p for p, (_, rhs) in rows.items())
    basis = [(1 << f) | sum(((row >> f) & 1) << p for p, (row, _) in rows.items())
             for f in range(width) if f not in rows]
    return origin, basis


class ImageFactor:
    """Factor a column map once; retain both image representatives and its kernel."""
    def __init__(self, columns):
        self.pivots, self.kernel = {}, []
        for index, column in enumerate(columns):
            residual, preimage = self.reduce(column)
            label = preimage ^ (1 << index)
            if residual:
                self.pivots[residual.bit_length() - 1] = residual, label
            else:
                self.kernel.append(label)
        self.rank = len(self.pivots)

    def reduce(self, value):
        """Return the quotient residual and an inside-space preimage."""
        label = 0
        for p in sorted(self.pivots, reverse=True):
            if value >> p & 1:
                row, preimage = self.pivots[p]
                value, label = value ^ row, label ^ preimage
        return value, label


def right_data(half, coupled_checks):
    """Canonical input checks -> physical right-half check columns P[:,t]."""
    offsets = tuple(row & ((1 << half) - 1) for row in coupled_checks)
    columns = tuple(sum(parity((row >> half) & polar_transform(1 << t, half)) << q
                        for q, row in enumerate(coupled_checks)) for t in range(half))
    return offsets, columns


def extend(half, left, offsets, columns, cutoff, counts=None):
    """Algorithm 1's exact extension: outside supports, one preimage, ALL lifts.

    Plotkin: x = (a+b, b), wt(x) = wt(a) + 2 wt(b outside supp(a)).
    P_S xi = rhs + P_O z is consistent iff its quotient residual vanishes.
    """
    a = polar_transform(left, half)
    inside = tuple(set_bits(a))
    outside = tuple(t for t in range(half) if not a >> t & 1)
    radius = (cutoff - len(inside)) // 2
    if radius < 0:
        return []
    rhs = sum(parity(row & left) << q for q, row in enumerate(offsets))
    factor = ImageFactor([columns[t] for t in inside])
    kernel = [embed(v, inside) for v in factor.kernel]
    if counts is not None:
        # Table I S_Q counts the specified direct loop without executing it.
        outside_count = sum(comb(len(outside), t) for t in range(min(radius, len(outside)) + 1))
        counts["direct_candidates"] += outside_count * (1 << len(inside))
        counts["outside_tests"] += outside_count
        counts["factorizations"] += 1
    outputs = []
    for weight in range(min(radius, len(outside)) + 1):
        for support in combinations(outside, weight):
            syndrome, z = rhs, 0
            for t in support:
                syndrome ^= columns[t]
                z |= 1 << t
            residual, xi = factor.reduce(syndrome)
            if residual:
                continue
            origin = embed(xi, inside) ^ z
            for direction in gray_span(kernel):
                beta = polar_transform(origin ^ direction, half)
                outputs.append(left | (beta << half))
    if counts is not None:
        counts["quotient_lifts"] += len(outputs)
    return outputs


def complete_directory(n, checks, cutoff, *, counts=None, shortcuts=True):
    """Algorithm 1: complete input list with wt(u G_N) <= cutoff, including zero.

    Exact shortcuts avoid a large zero-left support ball during construction.
    Set shortcuts=False for the paper's prescribed recursive S_Q counting scope.
    """
    if cutoff < 0:
        return []
    if n == 1:
        return [u for u in (0, 1) if u <= cutoff and all(not parity(row & u) for row in checks)]
    if shortcuts:
        # Every possible nonzero input row belongs to this polar mother code.
        pivots = {row.bit_length() - 1 for row in checks}
        active = set(range(n)) - pivots
        active.update(j for row in checks for j in set_bits(row ^ (1 << (row.bit_length() - 1))))
        if not active or cutoff < min(1 << j.bit_count() for j in active):
            return [0]
    half = n // 2
    left_checks = tuple(row for row in checks if row.bit_length() <= half)
    coupled = tuple(row for row in checks if row.bit_length() > half)
    left_words = complete_directory(half, left_checks, cutoff, counts=counts, shortcuts=shortcuts)
    offsets, columns = right_data(half, coupled)
    outputs = []
    for left in left_words:
        if shortcuts and left == 0:
            # a=0 gives x=(b,b); the pure-right cutoff is floor(D/2).
            right = complete_directory(half, tuple(row >> half for row in coupled),
                                       cutoff // 2, counts=counts, shortcuts=True)
            outputs.extend(beta << half for beta in right)
        else:
            outputs.extend(extend(half, left, offsets, columns, cutoff, counts))
    return outputs


def relative_quotient(left, half, check_columns, redundancy):
    """Section III: Q_a annihilates im(P_S); keep the exact lift multiplicity."""
    a = polar_transform(left, half)
    inside = tuple(set_bits(a))
    _, annihilator = solve_binary([(check_columns[t], 0) for t in inside], redundancy)
    def project(label):
        return sum(parity(label & row) << q for q, row in enumerate(annihilator))
    row_images = tuple(project(1 << q) for q in range(redundancy))
    outside_images = tuple(project(check_columns[t]) for t in range(half) if not a >> t & 1)
    lift_dimension = len(inside) - redundancy + len(annihilator)
    return row_images, outside_images, lift_dimension


def relative_spectrum(left, half, offsets, columns, cutoff):
    """Theorem 2: exact low-weight counts for one nonzero left input (a+b,b)."""
    weight = polar_transform(left, half).bit_count()
    radius = (cutoff - weight) // 2
    if radius < 0:
        return {}
    row_images, errors, lift_dimension = relative_quotient(left, half, columns, len(offsets))
    rhs = sum(parity(row & left) << q for q, row in enumerate(offsets))
    syndrome = 0
    for q in set_bits(rhs):
        syndrome ^= row_images[q]
    result = Counter()
    for t in range(min(radius, len(errors)) + 1):
        for support in combinations(range(len(errors)), t):
            label = 0
            for j in support:
                label ^= errors[j]
            if label == syndrome:
                result[weight + 2 * t] += 1 << lift_dimension
    return dict(result)


def quotient_exclusion(left, half, columns, redundancy, cutoff):
    """Theorem 2: build Q_a e_j and the EXACT distinct forbidden set B_a(D)."""
    weight = polar_transform(left, half).bit_count()
    radius = (cutoff - weight) // 2
    row_images, errors, lift_dimension = relative_quotient(left, half, columns, redundancy)
    forbidden = set()
    for t in range(min(radius, len(errors)) + 1):
        for support in combinations(range(len(errors)), t):
            label = 0
            for j in support:
                label ^= errors[j]
            forbidden.add(label)
    return row_images, forbidden, lift_dimension


def affine_region_safe(offsets, maps, forbidden, origin, directions):
    """Theorem 3: no forbidden label intersects any affine observation image.

    maps[g][p] is column p of L_g. This checks an affine domain without
    enumerating its points, and also handles nonzero affine offsets.
    """
    for c, columns, bad in zip(offsets, maps, forbidden):
        def image(value):
            result = 0
            for p in set_bits(value):
                result ^= columns[p]
            return result
        factor = ImageFactor([image(v) for v in directions])
        center = c ^ image(origin)
        if any(factor.reduce(label ^ center)[0] == 0 for label in bad):
            return False
    return True
