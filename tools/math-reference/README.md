# `tools/math-reference` — what `pantometry_core::math` is held to

`generate.py` writes four binary tables, and all four are **committed**:

| what | where | what reads it |
| --- | --- | --- |
| `exp-reference.bin` | `crates/pantometry-core/tests/data/` | `exp_and_ln_against_mpmath.rs` |
| `ln-reference.bin` | `crates/pantometry-core/tests/data/` | the same test |
| `sin-reference.bin` | `crates/pantometry-core/tests/data/` | `sin_and_cos_against_mpmath.rs` |
| `cos-reference.bin` | `crates/pantometry-core/tests/data/` | the same test, at the same points |

The two tests share their recipe, table reader and comparison in `tests/reference/mod.rs`.

```sh
python tools/math-reference/generate.py                        # rewrites the four tables
python tools/math-reference/generate.py --scale 10 --out DIR   # ten times the random points, elsewhere
python tools/math-reference/generate.py --exp-table            # prints math.rs's EXP_TABLE
python tools/math-reference/generate.py --ln-table             # prints math.rs's LN_TABLE
python tools/math-reference/generate.py --constants            # prints math.rs's other constants
```

It needs mpmath (tested with 1.3.0) and nothing else. Every value is computed at 50 significant
digits and rounded to a double by integer arithmetic. mpmath's own conversion goes through
`math.ldexp`, which would round a subnormal twice. The script's docstring has the table layout and
the recipe for the points.

## This is the third exception to "nothing generated is committed"

The first two are `tools/presets` and `tools/parts`. Like them, these files are **inputs**: the
test cannot run without them, and making them needs Python with mpmath, which neither CI nor a
contributor's `cargo test` has.

**The points are not stored, only the answers.** The test regenerates the inputs from the same
integer recipe and checks an FNV-1a digest of them against the one each table starts with. If the
recipe drifts on one side, the test fails and names the digest, rather than comparing against the
wrong points.

**A row's `frac` is 0 only where the value is exact.** The test reads from its sign which side of
the rounded double the exact value is on, and just below a power of two the spacing there is half
the rounded value's. Until the sin and cos tables were added, a `frac` nearer 0 than `2⁻¹⁶` ulp
was stored as 0, and a value one whole ulp off measured as half an ulp. 2048 of sin's rows and
10 473 of cos's are that close to `±1` from inside, where `1 − 2⁻⁵³` in place of 1 would have
passed. Such a `frac` is now stored as `±2⁻¹⁵`, and below
`|x| = 2⁻¹⁰` mpmath works at `2 log₂(1/|x|)` more bits, since at 50 digits `exp(−10⁻²⁰¹)` is
exactly 1 and the side is lost before the rounding. Regenerating moved 2905 rows of the exp table
and 169 of ln's, every one in `frac` alone. Each test checks, by `table_marks_exact_values`, that
the rows with `frac` 0 are exactly the points where the function is exact.

**They are kept out of the published kernel.** The tables come to 2.4 MB, and every user of
`pantometry-core` would download them. `Cargo.toml` excludes `tests/data/`, the two tests that
read it and the module they share, so a published crate holds no test that cannot compile. CI and
both gates run the tests from the repository, on Linux, macOS, Windows, release and
`wasm32-wasip1`, where the pinned digests of `exp`, `ln` and `sin_cos_turns` over fixed inputs are
what make "bit-identical on every platform" a check rather than a claim.

## When to run it

Run it when `math.rs`'s algorithm, its constant tables or the test's point recipe change. Changing
the recipe changes the digest, and the test then refuses the old tables until they are regenerated.
The `--exp-table`, `--ln-table` and `--constants` output must equal the constants in `math.rs` word
for word. That is how they were checked, and it is the first thing to check after editing either.
They are not checked against this script in CI: `math.rs`'s own tests hold every table entry and
constant to a closed form summed there in double-double — `2^(j/128)` by its Taylor series and by
exact products, `ln c` by the `atanh` series, `ln 2` by `Σ 1/(n 2ⁿ)`, `π/2` by Machin's formula
and each series coefficient from the last — and hold the `ln` table to the two properties the
method needs of it, cell by cell.

## The design the constants serve

`math.rs`'s `exp`, `ln` and `sin_cos_turns` are original to this workspace, written from the
published table-driven method (P. T. P. Tang, *ACM TOMS* 15(2), 1989, and 16(4), 1990),
Cody–Waite reduction, Dekker's exact sum and product, and Taylor's series; the module
documentation has the references and the design, and no library's code was read to write them. What the script computes for it:

| constant | what it is |
| --- | --- |
| `EXP_TABLE` | `2^(j/128)` for `j` in `0..128`, as the nearest double and the nearest double to the rest |
| `LN_TABLE` | per `F` from 0.703125 to 1.40625 (257 of them), `c`, `1/F` rounded to nine bits — 1 in the two cells about 1 — and `−ln c` as a multiple of `2⁻³⁵` and the nearest double to the rest |
| `LN2_HI`, `LN2_LO` | `ln 2` rounded to a multiple of `2⁻³⁵`, and the nearest double to the rest |
| `INV_STEP` | the nearest double to `128/ln 2` |
| `S1`, `S3`, `S5` (`_HI`, `_LO`) | the coefficients of `u`, `u³` and `u⁵` in `sin((π/2) u)`, `(−1)ᵏ (π/2)²ᵏ⁺¹/(2k + 1)!`, as the nearest double and the nearest double to the rest |
| `C2`, `C4` (`_HI`, `_LO`) | the coefficients of `u²` and `u⁴` in `cos((π/2) u)`, `(−1)ᵏ (π/2)²ᵏ/(2k)!`, likewise |
| `SIN_TAIL`, `COS_TAIL` | the coefficients of `u⁷` to `u¹⁹` and of `u⁶` to `u²⁰`, each the nearest double |

`sin_cos_turns` takes its argument in turns, so that its reduction, `4|t| = n + u`, is exact and
needs no table; it is in the module documentation with the rest. Its tables are mpmath's
`sinpi(2t)` and `cospi(2t)`, which reduce exactly too: `sin(2πt)` at 50 digits would round `2πt`
first, and give `−1.0 · 10⁻⁵¹` rather than `0` at `t = ½`.

**The point recipe was extended for this design** with its own seams — every boundary of `exp`'s
table index `(n + ½) ln 2/128` for `|n| ≤ 512`, `|x| = 708`, every edge of `ln`'s 256 cells at
`k = −1, 0, 1`, and `0.703125 · 2ᵉ` in every binade — appended after the points that were there,
whose rows in the regenerated tables are byte for byte the ones before. The thresholds and cell
edges of the first implementation are kept as points, so that nothing that was tested stops
being.
