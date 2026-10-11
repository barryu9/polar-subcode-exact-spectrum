# A guided route from equations to a regenerated code

## Start with `[256,64]`

Run `python reproduce_table_i.py --case n256_k64 --replays 1` from the repository
root. Open `profiles/n256_k64.json`, then read `construct()` in
`deterministic_design.py`. The profile supplies only information indices and
the target 48. Every frozen expression initially equals zero.

1. The smallest information-row weight is `d_P = 32`. The optional RM mother
   keeps rows of smaller weight frozen to zero. Its gap means that excluding
   weight 32 suffices for target 48.
2. The left length-128 constituent is unprecoded. The right profile has
   dimension 52 and target 24. Its seven released rows give 226 causal binary
   coefficients. `solve_direct()` builds a **complete common-supercode
   directory**, groups its 20,816 source/defect exclusions and performs the
   deterministic grouped search. This is the `V[128,52]` stage.
3. Fix these two constituents. `solve_cross_half()` lists all 28 critical
   left inputs and checks the pure-right distance bound. For each input,
   `relative_quotient()` computes an annihilator of the inside check image.
   With radius zero, the exact forbidden quotient set is `{0}`.
4. The same cross-half coefficients enter many exclusions. `joint_observation()`
   compresses their shared map from 144 coefficients to rank 87. The omitted
   57 dimensions are the uniform spectral kernel for these fixed constituents.
   Keeping representative columns preserves every attainable observation.
5. `prepare_conditions()` and `search_exclusions()` implement the affine
   exclusion system in Section IV. The initial full-system mean is `7/32`,
   so Lemma 1 ensures completion without backtracking. The deterministic DFS
   visits 88 prefixes, counting the root and the complete assignment.
6. Recover the cross-half `T` coefficients and verify canonical causality,
   dimension, RM containment, complete lower-layer exclusion and a weight-48
   witness. Only then compare the regenerated design with `codes/n256_k64.json`.

The other five cases follow the same `construct()` flow. The `[1024,512]`
case finishes its active search at `T[512,369]`; the final length-1024 lift is
static because its constituent distance bounds already imply target 24.

## Complete enumeration and the counting boundary

`complete_directory()` returns inputs, including zero. To obtain a spectrum,
count `polar_transform(u, N).bit_count()` over the returned list. For a fixed
nonzero left input, `relative_spectrum()` implements Theorem 2 directly:
each consistent outside support contributes `2^(|S_a|-rank(P_Sa))` inside
lifts. Repeated projected syndrome labels keep their physical-support
multiplicities. `extend()` constructs all these lifts explicitly.

For normal construction, two exact shortcuts avoid unnecessary work: a polar
mother lower bound can establish that only zero remains, and the zero-left
Plotkin branch is enumerated recursively with half the weight cutoff. The
paper's `S_Q` counting scope is audited **separately** with shortcuts disabled:

- For a final `V` stage, count every recursive extension of the complete common
  mother directory, including intermediate lifts and zero-left extensions.
- For a final `T` stage, release the editable cross checks to obtain the
  domain-induced common cover; count all critical-left extensions, including
  every lift. Count this cover, not just the fixed seed code.
- Exclude scalar base cases, factorization work and outside-syndrome tests
  from the ratio. They are still performed as needed by the algorithm.

The exact direct loop size per extension is
`sum(comb(|O_a|,t), t=0..radius) * 2^|S_a|`.
This is the loop of the paper's compared direct enumerator, not a lower bound
on every conceivable direct method. The quotient denominator is a freshly
generated count. Both the directory cardinality and its sorted-input SHA-256
are checked against the independent original experiment audit.

## Solving, guarantees and scope

The general binary `search_exclusions()` accepts affine offsets and arbitrary
finite forbidden sets, not only the `{0}` sets of the Table I T stages.
`conditional_mean()` implements the expected count of violated constraints
under a uniform remaining binary assignment. It does not assume independence
of violation events. `affine_region_safe()` checks Theorem 3 by testing whether
any forbidden label belongs to a residual affine observation image.

The small V stages group binary coefficients into patterns for one information
source. They use the local experiment's order: decreasing source row weight,
then decreasing source index; legal zero first, then sparse siblings; otherwise
exact remaining-domain-product lookahead. This is an exact implementation of
the common-supercode exclusions, with merged residual scopes and unit-scope
propagation. It does **not** use the T stage's binary mean ordering. Both
searches retain all alternatives and only report infeasibility on exhaustion.

The RM gap is used only when its stated parameter conditions hold. For the
`V[64,44]` constituent of `[128,64]`, target 8 exceeds the gap above minimum
4, so the algorithm excludes **all weights through 7**. The underlying
canonical directory, relative quotient, affine safety and binary solver
helpers themselves do not require RM containment. The Table I cross-half
constructor is deliberately limited to the paper's radius-zero experimental
domains; it rejects an unsupported nonzero-radius stage rather than silently
replacing its forbidden set by `{0}`.

Search budgets are optional in the standalone solver. Hitting a budget raises
`Inconclusive`, never returns `unsat`. A mean at least one never prunes, and a
mean below one is sufficient, not necessary: the `[256,128]` T stage has mean
`115/8` but still completes without backtracking. Distance certificates do not
establish decoding-performance gains.

## Inspect or use a saved design

No matrix package is needed. For example, from the repository root:

```python
import json
from pathlib import Path
from polar_subcodes import encode
from reproduce_table_i import expression_list

code = json.loads(Path("codes/n256_k64.json").read_text())
expressions = expression_list(code)
u, x = encode(1, code["N"], code["information_set"], expressions)
print("Input:", hex(u), "Codeword:", hex(x), "Weight:", x.bit_count())
```

An integer stores a binary vector with coordinate `i` in bit `i`. The JSON
expression lists are intended for reading; the integers are just a compact
standard-library representation of the same equations.
