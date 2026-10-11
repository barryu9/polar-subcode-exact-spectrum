"""Independent small-code checks of the paper's enumeration and certificates."""
from fractions import Fraction
import random
import unittest

from deterministic_design import (Inconclusive, conditional_mean, prepare_conditions,
                                  search_exclusions, uniform_spectral_kernel)
from polar_subcodes import (ImageFactor, affine_region_safe, canonical_checks,
                            complete_directory, encode, gray_span, parity,
                            polar_transform, quotient_exclusion, relative_spectrum, right_data,
                            solve_binary, validate_code)


def slow_transform(value, n):
    bits = [(value >> i) & 1 for i in range(n)]
    step = 1
    while step < n:
        for start in range(0, n, 2 * step):
            for j in range(start, start + step):
                bits[j] ^= bits[j + step]
        step *= 2
    return sum(bit << i for i, bit in enumerate(bits))


def image(columns, value):
    result = 0
    for p, column in enumerate(columns):
        if value >> p & 1:
            result ^= column
    return result


class PaperAlgorithmTests(unittest.TestCase):
    def test_transform_and_index_convention(self):
        for n in (1, 2, 4, 8):
            for u in range(1 << n):
                self.assertEqual(polar_transform(u, n), slow_transform(u, n))
                self.assertEqual(polar_transform(polar_transform(u, n), n), u)

    def test_complete_directories_against_all_messages(self):
        rng = random.Random(2027)
        for n in (2, 4, 8):
            for _ in range(12):
                info = [i for i in range(n) if rng.randrange(2)]
                expressions = [sum((rng.randrange(2) << i) for i in info if i < j)
                               if j not in info else 0 for j in range(n)]
                validate_code(n, info, expressions)
                checks = canonical_checks(n, info, expressions)
                all_words = [encode(m, n, info, expressions) for m in range(1 << len(info))]
                for cutoff in range(n + 1):
                    expected = {u for u, word in all_words if word.bit_count() <= cutoff}
                    for shortcuts in (False, True):
                        actual = complete_directory(n, checks, cutoff, shortcuts=shortcuts)
                        self.assertEqual(set(actual), expected)
                        self.assertEqual(len(actual), len(expected))

    def test_canonical_causality_rejection(self):
        with self.assertRaises(ValueError):
            validate_code(4, [1, 3], [2, 0, 8, 0])

    def test_linear_factor_keeps_all_kernel_lifts(self):
        columns = [3, 5, 6, 3]
        factor = ImageFactor(columns)
        for label in range(8):
            residual, origin = factor.reduce(label)
            exact = {v for v in range(16) if image(columns, v) == label}
            lifted = set() if residual else {origin ^ v for v in gray_span(factor.kernel)}
            self.assertEqual(lifted, exact)
        solution = solve_binary([(0b1011, 1), (0b0110, 0)], 4)
        origin, kernel = solution
        expected = {x for x in range(16) if parity(x & 0b1011) == 1 and parity(x & 0b0110) == 0}
        self.assertEqual({origin ^ v for v in gray_span(kernel)}, expected)

    def test_relative_spectrum_with_nonzero_offset(self):
        n, info = 8, [1, 3, 5, 7]
        expressions = [0, 0, 2, 0, 2, 0, 0b101000, 0]
        checks = canonical_checks(n, info, expressions)
        offsets, columns = right_data(4, [row for row in checks if row.bit_length() > 4])
        for left in (0b0110, 0b1000, 0b1110):
            for cutoff in range(9):
                expected = {}
                for right in range(16):
                    u = left | (right << 4)
                    if all(not parity(row & u) for row in checks):
                        weight = slow_transform(u, n).bit_count()
                        if weight <= cutoff:
                            expected[weight] = expected.get(weight, 0) + 1
                self.assertEqual(relative_spectrum(left, 4, offsets, columns, cutoff), expected)
                images, forbidden, _ = quotient_exclusion(left, 4, columns, len(offsets), cutoff)
                rhs = sum(parity(row & left) << q for q, row in enumerate(offsets))
                self.assertEqual(image(images, rhs) in forbidden, bool(expected))

    def test_exact_mean_and_complete_solver_against_all_assignments(self):
        rng = random.Random(6)
        for _ in range(30):
            width = 4
            maps = [[rng.randrange(8) for _ in range(width)] for _ in range(3)]
            offsets = [rng.randrange(8) for _ in maps]
            forbidden = [set(rng.sample(range(8), rng.randrange(5))) for _ in maps]
            conditions = prepare_conditions(maps, forbidden, offsets)
            solutions = [v for v in range(1 << width)
                         if all((c ^ image(m, v)) not in bad for c, m, bad in zip(offsets, maps, forbidden))]
            for length in range(width + 1):
                for prefix in range(1 << length):
                    violations = []
                    for tail in range(1 << (width - length)):
                        v = prefix | (tail << length)
                        violations.append(sum((c ^ image(m, v)) in bad for c, m, bad in zip(offsets, maps, forbidden)))
                    self.assertEqual(conditional_mean(conditions, prefix, length),
                                     Fraction(sum(violations), len(violations)))
            for ordering in ("mean", "zero"):
                answer = search_exclusions(conditions, width, ordering=ordering)
                self.assertEqual(answer["status"] == "sat", bool(solutions))
                if solutions:
                    self.assertIn(answer["value"], solutions)

    def test_mean_above_one_never_prunes_and_interrupt_is_inconclusive(self):
        conditions = prepare_conditions([[1]] * 3, [{0}] * 3)
        self.assertEqual(conditional_mean(conditions), Fraction(3, 2))
        self.assertEqual(search_exclusions(conditions, 1)["value"], 1)
        with self.assertRaises(Inconclusive):
            search_exclusions(conditions, 1, node_limit=1)

    def test_affine_region_and_uniform_kernel(self):
        maps, offsets, bad = [[1, 2, 3], [2, 1, 3]], [1, 2], [{0, 3}, {1}]
        for origin in range(8):
            for directions in ([], [1], [3, 4], [1, 2, 4]):
                exact = all(all((c ^ image(m, origin ^ v)) not in forbidden
                                for c, m, forbidden in zip(offsets, maps, bad))
                            for v in gray_span(directions))
                self.assertEqual(affine_region_safe(offsets, maps, bad, origin, directions), exact)
        kernel = uniform_spectral_kernel(maps, 3)
        self.assertEqual(set(gray_span(kernel)),
                         {v for v in range(8) if all(image(m, v) == 0 for m in maps)})


if __name__ == "__main__":
    unittest.main()
