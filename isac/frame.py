"""Transmit frame for the OFDM ISAC simulation: bits -> Gray-coded QAM -> N x M grid."""
from pathlib import Path

import numpy as np

from isac.params import OFDMRadarParams


def _gray(n):
    """Binary-reflected Gray code of an integer or integer array."""
    return n ^ (n >> 1)


def _check_order(order):
    k = int(np.log2(order))
    m = int(round(np.sqrt(order)))
    if 2**k != order or m * m != order:
        raise ValueError("order must be 4, 16, 64 or 256 (square QAM)")
    return k, m


def qam_constellation(order):
    """Gray-coded square QAM with unit average power. Index = bit label, value = point."""
    k, m = _check_order(order)
    half = k // 2
    levels = 2 * np.arange(m) - (m - 1)          # 16-QAM: -3, -1, 1, 3
    scale = np.sqrt(2 * (order - 1) / 3)         # average energy of the raw grid
    const = np.empty(order, dtype=complex)
    for i in range(m):
        for q in range(m):
            label = (_gray(i) << half) | _gray(q)
            const[label] = (levels[i] + 1j * levels[q]) / scale
    return const


def bits_to_symbols(bits, order):
    """Group bits into k-bit labels (MSB first) and look up their constellation points."""
    k, _ = _check_order(order)
    bits = np.asarray(bits).reshape(-1, k)
    weights = 1 << np.arange(k - 1, -1, -1)
    labels = bits @ weights
    return qam_constellation(order)[labels]


def symbols_to_bits(symbols, order):
    """Hard-decision demapping: slice I and Q to the nearest level, undo Gray code, output bits."""
    k, m = _check_order(order)
    half = k // 2
    scale = np.sqrt(2 * (order - 1) / 3)
    symbols = np.ravel(symbols)

    def level_index(x):
        return np.clip(np.rint((x * scale + (m - 1)) / 2), 0, m - 1).astype(int)

    labels = (_gray(level_index(symbols.real)) << half) | _gray(level_index(symbols.imag))
    shifts = np.arange(k - 1, -1, -1)
    return ((labels[:, None] >> shifts) & 1).astype(np.uint8).ravel()


def generate_frame(params, order=4, rng=None):
    """Random bits mapped to an N x M grid. Column m = OFDM symbol m, row n = subcarrier n."""
    rng = np.random.default_rng() if rng is None else rng
    k, _ = _check_order(order)
    bits = rng.integers(0, 2, params.N * params.M * k, dtype=np.uint8)
    symbols = bits_to_symbols(bits, order)
    X = symbols.reshape(params.N, params.M, order="F")   # fill column by column
    return bits, X


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    p = OFDMRadarParams()
    rng = np.random.default_rng(seed=1)
    Path("results").mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))

    for ax, order in zip(axes, (4, 16)):
        const = qam_constellation(order)

        # Gray check: every pair of nearest neighbours differs in exactly one bit
        dist = np.abs(const[:, None] - const[None, :])
        a, b = np.where(np.isclose(dist, dist[dist > 0].min()))
        gray_ok = all(bin(x ^ y).count("1") == 1 for x, y in zip(a, b))

        bits, X = generate_frame(p, order, rng)
        bits_rx = symbols_to_bits(X.reshape(-1, order="F"), order)
        errors = np.count_nonzero(bits != bits_rx)

        print(f"{order:>2}-QAM: frame {X.shape}, mean power {np.mean(np.abs(X)**2):.3f}, "
              f"Gray-coded: {gray_ok}, round-trip bit errors: {errors}")

        ax.scatter(const.real, const.imag)
        k = int(np.log2(order))
        for label, pt in enumerate(const):
            ax.annotate(format(label, f"0{k}b"), (pt.real, pt.imag),
                        textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8)
        ax.set_title(f"{order}-QAM (Gray-coded, unit power)")
        ax.set_xlabel("In-phase")
        ax.set_ylabel("Quadrature")
        ax.axhline(0, color="grey", lw=0.5)
        ax.axvline(0, color="grey", lw=0.5)
        ax.set_aspect("equal")
        ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig("results/constellations.png", dpi=150)
    print("Saved results/constellations.png")