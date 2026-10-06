# `tools/erfc-reference` — what `pantometry_forcefield::ewald::erfc` is made of and held to

`generate.py` writes no file. It prints, and what it prints is committed by hand in two places:

| command | what it prints | where it is committed |
| --- | --- | --- |
| `python tools/erfc-reference/generate.py --tables` | `ERF_SERIES`, `ERFCX_TERMS`, `ERFCX_TABLE`, `ERFCX_FAR` | `crates/pantometry-forcefield/src/ewald.rs`, as decimal literals |
| `python tools/erfc-reference/generate.py --points` | `REFERENCE`: 438 points `(x, hi, lo)`, the exact erfc as a double-double | `crates/pantometry-forcefield/tests/the_ewald_sum_against_closed_forms.rs` |
| `python tools/erfc-reference/generate.py --measure FILE` | the worst ulp error, by region, of `x_bits value_bits` lines | nowhere: a measurement, run once |

It needs mpmath (tested with 1.3.0) and nothing else; every value is computed at 50 digits.

**The tables in `ewald.rs` are `--tables`'s output word for word**, doc comments aside. They are
decimal literals rather than `f64::from_bits` because that is not `const` at the crate's MSRV,
1.78. Each literal is the shortest decimal that reads back as the double, so it is the same double.
Check that first after editing either side.

**How it was fitted, and what was not read.** erf's Taylor series below 0.4375 (Abramowitz and
Stegun 7.1.5). Above that, Chebyshev coefficients of `erfcx(x) = exp(x²) erfc(x)` on nine
intervals up to 6, by the discrete cosine sum at 64 nodes, cut below `2⁻⁶²` of the first. Beyond 6,
the same for `x erfcx(x)` in `u = 1/x²`. No library's code or table was read; the script's docstring
has the method.

**The measurement.** `cargo test --release -p pantometry-forcefield --test
the_ewald_sum_against_closed_forms erfc_dump_for_mpmath -- --ignored --nocapture` prints a million
points, and `--measure` on them gave a worst error of 3.21 ulp over the 992 240 with a normal result.
The test holds its 438 committed points to 3.5 ulp (`ERFC_ULP_BOUND`).

Run it when the method, the intervals or the point recipe change.
