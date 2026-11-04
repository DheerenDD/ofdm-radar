"""Radar processing: remove the data, window, and compute the range-Doppler map."""
from pathlib import Path

import numpy as np
from scipy.signal.windows import hann

from isac.params import OFDMRadarParams


def remove_data(Y, X):
    """Element-wise division by the known transmit symbols (reciprocal filtering)."""
    return Y / X


def range_doppler_map(params, F, window=True):
    """2D transform of the data-free grid F[n, m].

    IFFT over subcarriers (axis 0) turns the phase ramp into a peak at the target's range.
    FFT over symbols (axis 1) turns the phase rotation into a peak at the target's velocity.
    Returns the range axis [m], velocity axis [m/s] and the complex map (range x velocity).
    """
    if window:
        F = F * hann(params.N, sym=False)[:, None] * hann(params.M, sym=False)[None, :]

    rd = np.fft.ifft(F, axis=0)                          # range profiles
    rd = np.fft.fftshift(np.fft.fft(rd, axis=1), axes=1)  # Doppler, zero velocity centred

    ranges = np.arange(params.N) * params.range_resolution
    velocities = (np.arange(params.M) - params.M // 2) * params.velocity_resolution
    return ranges, velocities, rd


def to_db_above_noise(rd):
    """Power map in dB relative to the estimated noise floor."""
    power = np.abs(rd) ** 2
    noise_floor = np.median(power) / np.log(2)   # median of exponential noise = ln2 * mean
    return 10 * np.log10(power / noise_floor)


def plot_range_doppler(ranges, velocities, rd_db, targets=None, max_range=200, max_velocity=50,
                       title="Range-Doppler map", ax=None):
    """Heat map with optional ground-truth markers."""
    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax
    r_sel = ranges <= max_range
    v_sel = np.abs(velocities) <= max_velocity
    mesh = ax.pcolormesh(velocities[v_sel], ranges[r_sel], rd_db[np.ix_(r_sel, v_sel)],
                         shading="auto", cmap="viridis", vmin=-5, vmax=np.max(rd_db))
    if targets:
        ax.scatter([t.velocity_mps for t in targets], [t.range_m for t in targets],
                   marker="o", s=120, facecolors="none", edgecolors="red", linewidths=1.5,
                   label="ground truth")
        ax.legend(loc="upper right")
    ax.set_xlabel("Velocity [m/s] (positive = approaching)")
    ax.set_ylabel("Range [m]")
    ax.set_title(title)
    plt.colorbar(mesh, ax=ax, label="dB above noise floor")
    return ax


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from isac.frame import generate_frame
    from isac.targets import DEFAULT_SCENE, simulate_echo

    p = OFDMRadarParams()
    rng = np.random.default_rng(seed=3)
    snr_db = -10

    _, X = generate_frame(p, order=4, rng=rng)
    Y = simulate_echo(p, X, DEFAULT_SCENE, snr_db=snr_db, rng=rng)
    F = remove_data(Y, X)
    ranges, velocities, rd = range_doppler_map(p, F, window=True)
    rd_db = to_db_above_noise(rd)

    # Peak height next to each true target (search +/- 2 bins around the true position)
    print(f"Per-cell SNR {snr_db} dB, processing gain {10 * np.log10(p.N * p.M):.1f} dB")
    print(f"{'Target':<28}{'Peak above noise':>18}")
    for t in DEFAULT_SCENE:
        i = int(round(t.range_m / p.range_resolution))
        j = int(round(t.velocity_mps / p.velocity_resolution)) + p.M // 2
        peak = rd_db[max(i - 2, 0):i + 3, max(j - 2, 0):j + 3].max()
        label = f"{t.range_m:.0f} m, {t.velocity_mps:+.1f} m/s, {t.gain_db:+.0f} dB"
        print(f"{label:<28}{peak:>15.1f} dB")

    Path("results").mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_range_doppler(ranges, velocities, rd_db, targets=DEFAULT_SCENE, ax=ax,
                       title=f"OFDM radar range-Doppler map (B = {p.bandwidth / 1e6:.1f} MHz, Hann window)")
    fig.tight_layout()
    fig.savefig("results/range_doppler_map.png", dpi=150)
    print("Saved results/range_doppler_map.png")