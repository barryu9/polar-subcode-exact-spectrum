# Frozen expressions for polar subcodes: a reproducible tutorial

**Paper:** *Design Frozen Expressions for Polar Subcodes via Exact Low-weight Spectrum Characterization*

**Authors:** Junhua You, Shaohua Wu, Aolin Liu, Qinyu Zhang

**Status:** Submitted to the IEEE Information Theory Workshop (ITW 2027).

This repository is a small, standalone teaching implementation of the paper.
It connects three objects: polar subcodes with causal frozen expressions,
polar quotients that eliminate inside-support variables without losing kernel
lifts, and deterministic solving of the resulting exact exclusions. It
regenerates all six Table I designs from **zero expressions**, then verifies
their distances and structural metrics. The archived designs are checked only
after construction; they are never used to guide the search.

## Run all of Table I

Install Python **3.10 or later**. There are **no external packages**, compiled
extensions, downloads of reference designs, or platform-specific libraries.
After cloning this repository, run:

```sh
python reproduce_table_i.py
```

The default runs **five fresh full construction replays for each of the six
codes**. It also audits enumeration counts and the last T stage's zero-first
control. Results are written to `results/table_i.json`, per-code summaries and
per-replay designs. Run time depends on the machine; the full run takes about
one to two minutes on the authors' development machine. Every run checks:

1. Canonical causality, generator rank and the optional RM mother restriction.
2. Complete exclusion below the target, with the RM gap used only where valid.
3. A nonzero target-weight witness, proving equality of minimum distance.
4. Exact equality with the archived final expressions, enumeration-directory
   hashes, stage sizes, shared ranks, conditional means and search counts.

For a shorter first pass or one code:

```sh
python reproduce_table_i.py --replays 1
python reproduce_table_i.py --case n256_k64 --replays 1
python -m unittest discover -s tests -v
```

Keep Python's assertions enabled (do not use `python -O`). A failed check stops
the command with a nonzero exit status. The tests compare tiny codes and affine
systems with exhaustive independent enumeration, including nonzero offsets,
all kernel lifts, exact conditional means, safe affine regions and cases with
an initial mean above one.

## Follow the paper through the code

Read the files in this order; the complete design pipeline uses three Python
modules and does not require an experiment framework.

| Paper | Where to read | What to look for |
|---|---|---|
| II: canonical dynamic frozen expressions | `polar_subcodes.py`: `encode`, `canonical_checks`, `validate_code` | Zero-based indices, earlier information sources, `h_j = e_j + s_j`, dimension |
| II: common-supercode exclusion, Theorem 1 | `deterministic_design.py`: `solve_direct`, `SourcePatterns` | Release checks; list every low-weight input; exclude each surviving source/defect pattern |
| III: polar quotient, Theorem 2 and Algorithm 1 | `polar_subcodes.py`: `relative_quotient`, `quotient_exclusion`, `relative_spectrum`, `extend`, `complete_directory` | Exact forbidden labels, outside supports, quotient consistency, one preimage and **every** inside-kernel lift |
| III: shared observation and spectral freedom, Theorem 3 | `deterministic_design.py`: `joint_observation`, `uniform_spectral_kernel`; `polar_subcodes.py`: `affine_region_safe` | Observation rank, kernel directions for fixed constituents, exact affine safety |
| IV: deterministic solving, Algorithm 2 and Lemma 1 | `deterministic_design.py`: `prepare_conditions`, `conditional_mean`, `search_exclusions` | Exact affine exclusions; image-containment pruning; retain both children; mean ordering |
| V: Table I construction and verification | `deterministic_design.py`: `construct`; `reproduce_table_i.py`: `one_replay`, `count_enumeration` | Zero start, recursive right design, fixed left constituent, final certificates and exact metrics |

For a guided worked example and precise counting/timing boundaries, see
[the paper-to-code notes](docs/paper_to_code.md).

### Table I values

The timing column below is the paper's **archived** arithmetic mean of five
replays, in seconds. It is not the run time of this teaching implementation.

| `[N,K]` | `d_P -> d_target` | Active solve sequence | `q -> r` | `S_Q` | Initial `Phi` | `n_0 -> n_Phi` | Paper mean total (s) |
|---|---|---|---|---|---|---|---|
| `[64,32]` | `8 -> 12` | `V[64,32]` | N/A | `79991/3626 = 22.06` | N/A | N/A | 0.058 |
| `[128,64]` | `8 -> 16` | `V[64,44]` | N/A | `171886/84693 = 2.03` | N/A | N/A | 5.326 |
| `[256,64]` | `32 -> 48` | `V[128,52] -> T[256,64]` | `144 -> 87` | `2^26 = 67108864` | `7/32` | `92 -> 88` | 1.399 |
| `[256,128]` | `16 -> 24` | `V[128,89] -> T[256,128]` | `390 -> 302` | `2048` | `115/8` | `485 -> 303` | 5.768 |
| `[256,171]` | `8 -> 12` | `V[128,108] -> T[256,171]` | `756 -> 422` | `16` | `19/32` | `436 -> 423` | 1.552 |
| `[1024,512]` | `16 -> 24` | `T[256,223] -> T[512,369]` | `3504 -> 1772` | `16` | `19/65536` | `1802 -> 1773` | 2.166 |

`V` denotes the canonical expression **increment matrix** solved directly
from a complete common-supercode directory. `T` denotes the **cross-half
coupling matrix** with both constituents fixed. These are the paper's matrices,
not newly introduced method names. `[N,K]` after either symbol is the current
constituent's block length and dimension.

The **stage metrics** describe the **last active solve** in the displayed
sequence: the number of allowed binary coefficients `q`, observation rank `r`,
enumeration ratio `S_Q`, initial full-system mean `Phi`, and visited search
prefixes under zero-first and conditional-mean ordering (`n_0`, `n_Phi`). For
`[1024,512]`, this is the `[512,369]` solve; the length-1024 outer lift needs no
new solve. Total construction time includes every construction stage.
`Phi` counts expected violated exclusion groups, not codewords or decoding
errors. N/A denotes inapplicable T metrics.

`S_Q` is **direct enumeration candidates / explicitly generated quotient
lifts**, in the paper's specified counting scope. It is not a timing speedup.
For `[256,64]`, 28 critical left words each have 32 inside coordinates. The
direct loop checks `28 * 2^32 = 120259084288` candidates; the inside map has
rank 26, so quotient enumeration generates `28 * 2^6 = 1792` lifts. Their ratio
is exactly `2^26`. The script computes the direct count analytically and
actually generates the quotient lifts; it does not execute the 120-billion
candidate loop.

### Timing and evidence

`data/published_table_i.json` preserves all five original total times per code,
the exact count audits, stage metrics and hashes of the original local records.
The original timing platform was a single-threaded Intel Core i7-12700K with
32 GiB RAM. Its implementation used native-assisted enumeration and a larger
experimental dependency graph. This release instead favors transparent
standard-library Python. Its measured speed can be either faster or slower
for a given case and must **not** be presented as the paper's original timing.

The script reports archived paper means and **fresh teaching-implementation
means separately**. Fresh totals cover directories, shared observation
compression, solving, recovery and verification. File I/O, regression
comparisons, the separate `S_Q` audit and the zero-first control are outside
those totals. Five fresh replays rebuild every code independently.

### Information sets and saved codes

`profiles/` contains only the information sets used by the paper, target
distances and source-file metadata. The authors thank **Dr. Vera Miloslavskaya**
for making reference designs available on her
[code-specification website](https://sites.google.com/site/veradmiloslavskaya/specifications-of-error-correcting-codes).
The referenced published results establish that the selected distance targets
are attainable. Reference expression coefficients are discarded; this
repository does not redistribute those reference design files.

`codes/` contains **our independently constructed designs**. Every entry
`"j": [i_1, i_2, ...]` means `u_j = u_i1 XOR u_i2 XOR ...`; absent frozen
rows are zero. All indices are zero-based, and the polar transform has no bit
reversal. Each file also includes the information set and an exact-distance
witness. The reproduction command regenerates these expressions and checks
them against the archive.

## Scope and AI Statement

The public code is limited to the algorithms, six designs and experiments
reported in this paper. This is a tutorial, not a general-purpose high-speed
decoder or a large-scale search framework. The complete-directory and affine
exclusion helpers do not require an RM mother. The Table I construction uses
that **additional restriction** to simplify the lower-weight exclusion where
its gap applies; `[128,64]` still requires a complete cutoff-7 constituent
directory. General search and enumeration may remain exponential. A verified
minimum distance is not a claim of improved SCL error rate.

The authors used OpenAI's GPT-6 model to audit the paper's mathematical
rigor and organize an open-source repository for readers
to study and reproduce its algorithms and experiments
(see the repository footnote in Section V).
The authors independently developed the technical content and take
responsibility for all content in this paper.

**OpenAI's GPT-6 created this repository from the authors' local experimental
materials** to help readers study and reproduce the paper's algorithms and
experiments and help reviewers assess its accuracy.

The implementation is released under the [MIT License](LICENSE). Please cite
the paper when using its results; the submission has no proceedings DOI yet.
