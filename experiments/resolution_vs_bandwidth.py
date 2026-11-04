"""Experiment: how bandwidth decides whether two close objects can be separated in range.

Two objects 3 m apart with the same speed are measured with four bandwidths.
The subcarrier spacing stays at 120 kHz; only the number of subcarriers changes.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.signal.windows import hann

from isac.frame import generate_frame
from isac.params import OFDMRadarParams
from isac.processing import remove_data
from isac.targets import Target, simulate_echo

# Same phase for both echoes, so the result depends only on bandwidth, not on random interference
TARGETS = [Target(range_m=80.0, velocity_mps=10.0, phase_rad=0.0),
           Target(range_m=83.0, velocity_mps=10.0, phase_rad=0.0)]
SUBCARRIERS = [256, 512, 1024, 2048]
PAD = 8            # zero-padding factor: smoother curve for plotting, NOT better resolution


def range_profile(params, F, pad=PAD):
    """Range profile at the Doppler bin of the targets, zero-padded for a smooth plot."""
    F = F * hann(params.N, sym=False)[:, None] * hann(params.M, sym=False)[None, :]
    doppler = np.fft.fft(F, axis=1)
    target_bin = np.argmax(np.sum(np.abs(doppler) ** 2, axis=0))   # strongest velocity bin
    profile = np.abs(np.fft.ifft(doppler[:, target_bin], n=pad * params.N)) ** 2
    ranges = np.arange(pad * params.N) * params.range_resolution / pad
    return ranges, 10 * np.log10(profile / profile.max())


if __name__ == "__main__":
    rng = np.random.default_rng(seed=4)
    fig, ax = plt.subplots(figsize=(9, 5))

    print(f"Two objects at {TARGETS[0].range_m:.0f} m and {TARGETS[1].range_m:.0f} m "
          f"({TARGETS[1].range_m - TARGETS[0].range_m:.0f} m apart)\n")
    print(f"{'Subcarriers':>11} {'Bandwidth':>11} {'Resolution':>11}   Result")

    for N in SUBCARRIERS:
        p = OFDMRadarParams(N=N)
        _, X = generate_frame(p, order=4, rng=rng)
        Y = simulate_echo(p, X, TARGETS, snr_db=0, rng=rng)
        ranges, profile_db = range_profile(p, remove_data(Y, X))

        window = (ranges > 70) & (ranges < 95)
        peaks, _ = find_peaks(profile_db[window], height=-10, prominence=3)   # strong peaks, 3 dB dip
        peak_ranges = ranges[window][peaks]
        # Resolved = one peak within 1 m of EACH true object
        resolved = all(np.any(np.abs(peak_ranges - t.range_m) < 1.0) for t in TARGETS)
        result = "separated" if resolved else "merged"
        found = ", ".join(f"{r:.1f} m" for r in peak_ranges)
        print(f"{N:>11} {p.bandwidth / 1e6:>7.2f} MHz {p.range_resolution:>9.2f} m   {result:<10} peaks at {found}")

        ax.plot(ranges[window], profile_db[window],
                label=f"B = {p.bandwidth / 1e6:.0f} MHz (resolution {p.range_resolution:.2f} m): {result}")

    for t in TARGETS:
        ax.axvline(t.range_m, color="red", linestyle="--", linewidth=1)
    ax.set_ylim(-40, 3)
    ax.set_xlabel("Range [m]")
    ax.set_ylabel("Normalised power [dB]")
    ax.set_title("Two objects 3 m apart: more bandwidth separates them (red lines = true positions)")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower center", fontsize=9)
    fig.tight_layout()

    Path("results").mkdir(exist_ok=True)
    fig.savefig("results/resolution_vs_bandwidth.png", dpi=150)
    print("\nSaved results/resolution_vs_bandwidth.png")