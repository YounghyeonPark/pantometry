"""The constants of `pantometry_forcefield::ewald::erfc`, and the values it is held to.

    python tools/erfc-reference/generate.py --tables          # prints ewald.rs's ERF_SERIES, ERFCX_TABLE, ERFCX_FAR
    python tools/erfc-reference/generate.py --points          # prints the test's reference points
    python tools/erfc-reference/generate.py --measure FILE    # reads `x_bits value_bits` lines, reports the worst ulp

Needs mpmath (tested with 1.3.0) and nothing else. Every value is computed at 50 significant
digits and rounded to a double by mpmath's `float`, which is correct for every normal double --
and no constant here is subnormal.

**What is fitted, and how.** Nothing is read from a library or a paper's table:

- `ERF_SERIES` is Taylor's series of erf, `erf x = (2/sqrt(pi)) sum (-1)^n x^(2n+1) / (n! (2n+1))`
  (Abramowitz and Stegun 7.1.5), fourteen terms, for |x| < 0.4375.
- `ERFCX_TABLE` holds, for each interval `[a, b)` of `EDGES`, the Chebyshev coefficients of the
  scaled complement `erfcx(x) = exp(x^2) erfc(x)` on it, by the discrete cosine sum at 64
  Chebyshev nodes, cut where a coefficient falls below `2^-62` of the first.
- `ERFCX_FAR` is the same for `g(u) = x erfcx(x)` as a function of `u = 1/x^2` on `[0, 1/36]`,
  i.e. `x >= 6`, where `g(0) = 1/sqrt(pi)` -- the asymptotic expansion's leading term
  (A&S 7.1.23) -- and the function is smooth in `u`.

`--points` prints, for x drawn by a fixed recipe, the exact erfc as a double-double `(hi, lo)`, so
that the test can measure an error in ulps with no reference arithmetic of its own beyond one
exact subtraction.
"""

import argparse
import struct
import sys

import mpmath

mpmath.mp.dps = 50

EDGES = [0.4375, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6]
SMALL = 0.4375
FAR = 6.0


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def double(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


def erfcx(x):
    x = mpmath.mpf(x)
    return mpmath.erfc(x) * mpmath.exp(x * x)


def chebyshev(f, a, b, nodes=64):
    a = mpmath.mpf(a)
    b = mpmath.mpf(b)
    theta = [mpmath.pi * (k + mpmath.mpf(1) / 2) / nodes for k in range(nodes)]
    fx = [f((a + b) / 2 + (b - a) / 2 * mpmath.cos(t)) for t in theta]
    c = [2 * sum(fx[k] * mpmath.cos(j * theta[k]) for k in range(nodes)) / nodes for j in range(nodes)]
    c[0] /= 2
    tol = abs(c[0]) * mpmath.mpf(2) ** -62
    n = max(j for j in range(len(c)) if abs(c[j]) > tol)
    return [float(v) for v in c[: n + 1]]


def far(u):
    if u == 0:
        return 1 / mpmath.sqrt(mpmath.pi)
    x = 1 / mpmath.sqrt(u)
    return x * erfcx(x)


def tables():
    series = [
        float(2 / mpmath.sqrt(mpmath.pi) * (-1) ** n / (mpmath.factorial(n) * (2 * n + 1)))
        for n in range(14)
    ]
    rows = [chebyshev(erfcx, EDGES[i], EDGES[i + 1]) for i in range(len(EDGES) - 1)]
    tail = chebyshev(far, 0, mpmath.mpf(1) / 36)
    return series, rows, tail


def rust(x):
    """The shortest decimal that reads back as `x`: `f64::from_bits` is not `const` at the
    crate's MSRV, 1.78, and a correctly rounded literal is the same double."""
    r = repr(x)
    assert float(r) == x
    return r if ("." in r or "e" in r) else r + ".0"


def print_tables():
    series, rows, tail = tables()
    # The first term is 2/sqrt(pi), which clippy knows by name; it is printed as computed.
    print("#[allow(clippy::approx_constant)]")
    print("const ERF_SERIES: [f64; %d] = [" % len(series))
    for v in series:
        print("    %s," % rust(v))
    print("];")
    width = max(len(r) for r in rows)
    print("const ERFCX_TERMS: [usize; %d] = [%s];" % (len(rows), ", ".join(str(len(r)) for r in rows)))
    print("const ERFCX_TABLE: [[f64; %d]; %d] = [" % (width, len(rows)))
    for r in rows:
        print("    [")
        for v in r + [0.0] * (width - len(r)):
            print("        %s," % rust(v))
        print("    ],")
    print("];")
    print("const ERFCX_FAR: [f64; %d] = [" % len(tail))
    for v in tail:
        print("    %s," % rust(v))
    print("];")


def splitmix(seed, i):
    mask = (1 << 64) - 1
    z = (seed + (i + 1) * 0x9E3779B97F4A7C15) & mask
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & mask
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & mask
    return z ^ (z >> 31)


def recipe():
    """The points: every seam of the method and its neighbours, and a uniform spread."""
    xs = []
    for e in EDGES + [SMALL, FAR, 0.0, 26.5, 27.0]:
        for s in (1, -1):
            b = bits(s * e)
            xs += [double(b), double(b + 1)] + ([double(b - 1)] if e != 0 else [])
    for i in range(240):
        # uniform in [-3, 27)
        xs.append(-3.0 + 30.0 * (splitmix(0xE7FC, i) >> 11) / float(1 << 53))
    for i in range(120):
        # denser where Ewald lives, [0, 5)
        xs.append(5.0 * (splitmix(0xE7FD, i) >> 11) / float(1 << 53))
    out = []
    for x in xs:
        t = mpmath.erfc(mpmath.mpf(x))
        if t == 0 or abs(t) < mpmath.mpf(2) ** -1022:
            continue
        hi = float(t)
        lo = float(t - mpmath.mpf(hi))
        out.append((x, hi, lo))
    return out


def print_points():
    pts = recipe()
    print("const REFERENCE: [(u64, u64, u64); %d] = [" % len(pts))
    for x, hi, lo in pts:
        print("    (0x%016x, 0x%016x, 0x%016x)," % (bits(x), bits(hi), bits(lo)))
    print("];")


def measure(path):
    worst = {}
    count = 0
    for line in open(path):
        parts = line.split()
        if len(parts) != 2:
            continue
        x = double(int(parts[0], 16))
        v = double(int(parts[1], 16))
        t = mpmath.erfc(mpmath.mpf(x))
        if t == 0 or abs(t) < mpmath.mpf(2) ** -1022:
            continue
        count += 1
        e = int(mpmath.floor(mpmath.log(abs(t), 2)))
        ulp = mpmath.mpf(2) ** (e - 52)
        err = float(abs(mpmath.mpf(v) - t) / ulp)
        band = (
            "x < -0.4375" if x <= -SMALL else
            "|x| < 0.4375" if x < SMALL else
            "0.4375 <= x < 6" if x < FAR else
            "x >= 6"
        )
        if err > worst.get(band, (0.0, 0.0))[0]:
            worst[band] = (err, x)
    print("points:", count)
    for band, (err, x) in sorted(worst.items()):
        print("%-18s worst %.4f ulp at x = %r" % (band, err, x))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tables", action="store_true")
    p.add_argument("--points", action="store_true")
    p.add_argument("--measure")
    a = p.parse_args()
    if a.tables:
        print_tables()
    elif a.points:
        print_points()
    elif a.measure:
        measure(a.measure)
    else:
        p.print_help()
        sys.exit(2)


if __name__ == "__main__":
    main()
