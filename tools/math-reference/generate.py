"""The reference values `pantometry_core::math::{exp, ln}` are held to.

    python tools/math-reference/generate.py            # writes the two committed tables
    python tools/math-reference/generate.py --scale 10 --out DIR  # ten times the random points
    python tools/math-reference/generate.py --exp-table   # prints math.rs's EXP_TABLE
    python tools/math-reference/generate.py --ln-table    # prints math.rs's LN_TABLE
    python tools/math-reference/generate.py --constants   # prints math.rs's LN2_HI, LN2_LO, INV_STEP

The constants in `crates/pantometry-core/src/math.rs` are what the last three print, word for
word; the reference tables are what the test reads.

Needs mpmath (tested with 1.3.0) and nothing else. Every value is computed at 50 significant
digits, about 166 bits, and rounded to a double here, by integer arithmetic, rather than by
mpmath's own conversion -- which goes through `math.ldexp` and would round a subnormal twice.

**The points are not stored.** They are drawn from a recipe of integer operations on bit patterns
(splitmix64, ranges of bit patterns, a fixed list of bases with their neighbours) which
`crates/pantometry-core/tests/exp_and_ln_against_mpmath.rs` repeats. Each table starts with an
FNV-1a digest of the inputs it was made for, and the test recomputes it: a recipe that drifted on
one side is a failure that names the digest, not a comparison of the wrong points.

Per point the table holds the correctly rounded value's bits (u64) and where the exact value lies
from it, `(exact - rounded) / ulp(exact)`, in units of 2^-15 (i16). The ulp is the spacing of
doubles at the exact value, `2^(max(floor(log2 |t|), -1022) - 52)`; the test rebuilds it from the
rounded value and the sign of that fraction.

Layout, little-endian: b"PMREF1\\0\\0", u64 input digest, u32 count, then count x (u64, i16).
"""

import argparse
import math
import os
import struct
import sys

import mpmath

mpmath.mp.dps = 50

MASK = (1 << 64) - 1


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def double(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


def splitmix(seed, i):
    z = (seed + (i + 1) * 0x9E3779B97F4A7C15) & MASK
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
    return z ^ (z >> 31)


def uniform(lo, hi, seed, n, sign=0):
    """n bit patterns uniform in [lo, hi], with the sign bit `sign`; sign 2 alternates."""
    out = []
    for i in range(n):
        b = lo + splitmix(seed, i) % (hi - lo + 1)
        s = (i & 1) if sign == 2 else sign
        out.append(b | (s << 63))
    return out


def around(base, reach, signs=(0,)):
    return [((base + j) & MASK) | (s << 63) for s in signs for j in range(-reach, reach + 1)]


SIGN = 1 << 63
LN_OFF = 0x3FE6000000000000  # 0.6875, where the first implementation's cells began
LN_2 = 0.6931471805599453  # std::f64::consts::LN_2, bits 0x3fe62e42fefa39ef


def exp_inputs(scale):
    xs = []
    # The whole domain, both signs, log-uniform in |x| from 2^-30 to past each threshold.
    xs += uniform(bits(2.0**-30), bits(746.0), 1, 16384 * scale, sign=2)
    # Densely where the result is neither 1 nor 0 nor infinite: |x| from 1 to the thresholds.
    xs += uniform(bits(1.0), bits(709.8), 2, 8192 * scale, sign=0)
    xs += uniform(bits(1.0), bits(745.2), 3, 8192 * scale, sign=1)
    # Every bit pattern that is not a NaN: almost all overflow, underflow or return 1.
    xs += [b for b in uniform(0, MASK, 4, 4096 * scale) if not math.isnan(double(b))]
    # Every reduction boundary: k changes at (k + 1/2) ln 2. Four doubles either side.
    for k in range(-1075, 1025):
        xs += around(bits((k + 0.5) * LN_2), 4)
    # Branch thresholds, each side, both signs: |x| = 1/2 ln 2, 3/2 ln 2 and 2^-28, which were the
    # first implementation's and are kept so that the points stay the same; 708.39 (where the
    # result turns subnormal); and the overflow and underflow limits.
    for base in [
        0x3FD62E4300000000,
        0x3FF0A2B200000000,
        0x3E30000100000000,
        0x4086232B00000000,
        bits(708.3964185322641),
        0x40862E42FEFA39EF,
        0x40874910D52D3051,
    ]:
        xs += around(base, 64, signs=(0, 1))
    # The subnormal results.
    xs += uniform(bits(708.3964185322641), bits(745.1332191019412), 5, 2048 * scale, sign=1)
    # The seams of math.rs's method that the lines above do not reach: every boundary of the
    # table index j, (n + 1/2) ln2/128, for |n| <= 512, two doubles either side; and |x| = 708,
    # where the main path ends, both signs.
    for n in range(-512, 512):
        xs += around(bits((n + 0.5) * LN_2 / 128), 2)
    xs += around(bits(708.0), 64, signs=(0, 1))
    return xs


def ln_inputs(scale):
    xs = []
    # Every positive finite double, log-uniform.
    xs += uniform(1, 0x7FEFFFFFFFFFFFFF, 11, 16384 * scale)
    # Every subnormal.
    xs += uniform(1, 0x000FFFFFFFFFFFFF, 12, 2048 * scale)
    # About 1, where ln x is small and its relative accuracy is the hard part: 4096 doubles each
    # side, and a wider band uniform in bits from 15/16 to 17/16.
    xs += around(bits(1.0), 4096)
    xs += uniform(bits(0.9375), bits(1.0625), 13, 8192 * scale)
    # Every binade's reduction boundary: k steps where the mantissa crosses sqrt(2), high word
    # 0x6a09e. Three doubles either side, in every normal binade.
    for e in range(1, 2047):
        xs += around((e << 52) | (0x6A09E << 32), 3)
    # The descreening integral's range, ln(U/L) with U/L from 1 to 64.
    xs += uniform(bits(1.0), bits(64.0), 14, 8192 * scale)
    # Extremes: the smallest normal, the largest finite, the smallest subnormals, and every power
    # of two (f = 0, ln x = k ln 2).
    xs += around(0x0010000000000000, 64)
    xs += [0x7FEFFFFFFFFFFFFF - j for j in range(129)]
    xs += list(range(1, 65))
    xs += [bits(math.ldexp(1.0, k)) for k in range(-1074, 1024)]
    # Every edge of the first implementation's 128 cells of [0.6875, 1.375), at k = -1, 0 and 1,
    # eight doubles either side: kept, so that the points stay the same, and dense about 1.
    for k in (-1, 0, 1):
        for i in range(129):
            xs += around(((LN_OFF + (i << 45)) + (k << 52)) & MASK, 8)
    # The seams of math.rs's method: every edge of its 256 cells, F's bits plus half a step, at
    # k = -1, 0 and 1, two doubles either side; and where k steps, 0.703125 2^e, in every binade,
    # one either side.
    for k in (-1, 0, 1):
        for i in range(256):
            xs += around(((LN_F_FIRST + (i << LN_F_SHIFT) + (1 << (LN_F_SHIFT - 1))) + (k << 52)) & MASK, 2)
    for e in range(1, 2047):
        xs += around((e << 52) | 0x0006800000000000, 1)
    return xs


def fnv1a(words):
    h = 0xCBF29CE484222325
    for w in words:
        for byte in struct.pack("<Q", w):
            h ^= byte
            h = (h * 0x100000001B3) & MASK
    return h


def rounded(t):
    """(the double nearest t, (t - it)/ulp(t) as a float). t an mpf."""
    if t == 0:
        return 0.0, 0.0
    big = mpmath.mpf(2) ** 1024 - mpmath.mpf(2) ** 970
    if abs(t) >= big:
        return math.copysign(math.inf, t), 0.0
    _, e = mpmath.frexp(t)  # |t| = m 2^e, 1/2 <= m < 1
    q = max(int(e) - 1, -1022) - 52
    unit = mpmath.ldexp(mpmath.mpf(1), q)
    n = int(mpmath.nint(t / unit))
    hi = math.ldexp(float(n), q)  # n < 2^54, so exact
    assert mpmath.mpf(hi) == n * unit
    frac = float((t - n * unit) / unit)
    assert abs(frac) <= 0.5
    return hi, frac


def table(name, xs, function):
    digest = fnv1a(xs)
    body = []
    for b in xs:
        hi, frac = rounded(function(mpmath.mpf(double(b))))
        body.append(struct.pack("<Qh", bits(hi), int(round(frac * 32768.0))))
    print(f"{name}: {len(xs)} points, input digest {digest:016x}", file=sys.stderr)
    return b"PMREF1\0\0" + struct.pack("<QI", digest, len(xs)) + b"".join(body)


def rust_words(words, per_line):
    """u64 words as Rust hex literals, `per_line` to a line."""
    lines = []
    for i in range(0, len(words), per_line):
        lines.append(
            "    "
            + " ".join(
                "0x%04x_%04x_%04x_%04x," % (w >> 48, (w >> 32) & 0xFFFF, (w >> 16) & 0xFFFF, w & 0xFFFF)
                for w in words[i : i + per_line]
            )
        )
    return "\n".join(lines)


# The design these tables serve is stated in `math.rs`'s module documentation; in short:
#
# exp: x = (128 k + j) ln2/128 + r, |r| <= ln2/256, exp x = 2^k 2^(j/128) e^r (Tang, ACM TOMS
#      15(2), 1989). The table is 2^(j/128) as a double-double, T_hi + T_lo.
# ln:  x = 2^k m, m in [0.703125, 1.40625); F = m rounded to 8 bits after the leading one, so
#      F = 1 when m is near 1; c is a short approximation of 1/F, and
#      ln x = k ln2 - ln c + ln(1 + r), r = m c - 1 = (F c - 1) + (m - F) c, which is exact
#      (Tang, ACM TOMS 16(4), 1990, for the table-driven form). The table is, per F, c and -ln c
#      as a double-double whose high part is a multiple of 2^-35, so that k ln2_hi + (-ln c)_hi
#      is exact.

LN_F_FIRST = 0x3FE6800000000000  # 0.703125, the first F
LN_F_SHIFT = 44  # F has 8 fraction bits: consecutive F are 1 << 44 apart in their bits
LN_F_COUNT = 257  # F from 0.703125 to 1.40625: 152 in [1/2, 1) at 2^-9, 105 in [1, 2) at 2^-8
HI_QUANTUM = 35  # ln2_hi and every (-ln c)_hi are multiples of 2^-35
# The cell just below 1, m in [1 - 3 2^-10, 1 - 2^-10): c = 1 there rather than the nearest
# 1 + 2^-8, so that r = m - 1 and nothing cancels where ln x is as small as 2^-10.
LN_C_IS_ONE = (0x3FEFF00000000000, 0x3FF0000000000000)  # 511/512 and 1


def to_quantum(t, q):
    """The multiple of 2^-q nearest t, as a double (exact: |t| < 2^18)."""
    n = int(mpmath.nint(t * mpmath.mpf(2) ** q))
    v = math.ldexp(float(n), -q)
    assert mpmath.mpf(v) == n * mpmath.mpf(2) ** -q
    return v


def exp_table():
    """`math::EXP_TABLE`: per j in 0..128, the bits of T_hi and of T_lo, 2^(j/128) = T_hi + T_lo."""
    words = []
    for j in range(128):
        t = mpmath.power(2, mpmath.mpf(j) / 128)
        hi = rounded(t)[0]
        lo = rounded(t - mpmath.mpf(hi))[0]
        words += [bits(hi), bits(lo)]
    return words


def ln_c(fb):
    """c for the F whose bits are fb: 1/F rounded to a multiple of 2^-9 below 1 and of 2^-8 at or
    above it -- nine significant bits either way -- or 1 in the cells LN_C_IS_ONE names."""
    if fb in LN_C_IS_ONE:
        return 1.0
    inv = 1 / mpmath.mpf(double(fb))
    q = 9 if inv < 1 else 8
    return to_quantum(inv, q)


def ln_table():
    """`math::LN_TABLE`: per F, the bits of c, of (-ln c)_hi and of (-ln c)_lo.

    F's bits are LN_F_FIRST + (i << LN_F_SHIFT) for i in 0..LN_F_COUNT.
    """
    words = []
    for i in range(LN_F_COUNT):
        c = ln_c(LN_F_FIRST + (i << LN_F_SHIFT))
        g = -mpmath.log(mpmath.mpf(c))
        ghi = to_quantum(g, HI_QUANTUM)
        glo = rounded(g - mpmath.mpf(ghi))[0]
        words += [bits(c), bits(ghi), bits(glo)]
    return words


def constants():
    """`math`'s scalar constants, as `const NAME: u64 = bits;` lines."""
    ln2 = mpmath.log(2)
    hi = to_quantum(ln2, HI_QUANTUM)
    out = [
        ("LN2_HI", bits(hi)),
        ("LN2_LO", bits(rounded(ln2 - mpmath.mpf(hi))[0])),
        ("INV_STEP", bits(rounded(128 / ln2)[0])),
    ]
    return "\n".join(
        "const %s: u64 = 0x%04x_%04x_%04x_%04x;"
        % (name, w >> 48, (w >> 32) & 0xFFFF, (w >> 16) & 0xFFFF, w & 0xFFFF)
        for name, w in out
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scale", type=int, default=1)
    parser.add_argument("--out", default=None)
    parser.add_argument("--exp-table", action="store_true", help="print math::EXP_TABLE and stop")
    parser.add_argument("--ln-table", action="store_true", help="print math::LN_TABLE and stop")
    parser.add_argument("--constants", action="store_true", help="print math's constants and stop")
    args = parser.parse_args()
    if args.exp_table:
        print(rust_words(exp_table(), 2))
        return
    if args.ln_table:
        print(rust_words(ln_table(), 3))
        return
    if args.constants:
        print(constants())
        return
    here = os.path.dirname(os.path.abspath(__file__))
    out = args.out or os.path.join(here, "..", "..", "crates", "pantometry-core", "tests", "data")
    os.makedirs(out, exist_ok=True)

    def exp(t):
        return mpmath.exp(t)

    def ln(t):
        if t == 0:
            raise ValueError("ln 0 is not in a table")
        return mpmath.log(t)

    for name, xs, f in [("exp", exp_inputs(args.scale), exp), ("ln", ln_inputs(args.scale), ln)]:
        with open(os.path.join(out, f"{name}-reference.bin"), "wb") as fh:
            fh.write(table(name, xs, f))


if __name__ == "__main__":
    main()
