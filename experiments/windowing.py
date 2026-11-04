"""Experiment: why the window matters when a weak object sits next to a strong one.

A strong object (car) and a weak object (pedestrian, 30 dB weaker) are 6 m apart.
The same received frame is processed twice: without a window and with a Hann window.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from isac.detection import detect_targets, match_to_truth
from isac.frame import generate_frame
from isac.params import OFDMRadarParams
from isac.processing import range_doppler_map, remove_data, to_db_above_noise
from isac.targets import Target, simulate_echo

TARGETS = [
    Target(range_m=80.0, velocity_mps=10.0, gain_db=0.0),     # strong: car
    Target(range_m=86.0, velocity_mps=10.0, gain_db=-30.0),   # weak: pedestrian, 30 dB weaker
]
SNR_DB = 0          # per-cell SNR of the strong object: high, so sidelobes stand above the noise
THRESHOLD_DB = 15


if __name__ == "__main__":
    p = OFDMRadarParams()
    rng = np.random.default_rng(seed=5)
    _, X = generate_frame(p, order=4, rng=rng)
    F = remove_data(simulate_echo(p, X, TARGETS, snr_db=SNR_DB, rng=rng), X)

    fig, ax = plt.subplots(figsize=(9, 5))
    print(f"Strong object: {TARGETS[0].range_m:.0f} m, weak object: {TARGETS[1].range_m:.0f} m "
          f"({TARGETS[1].gain_db:.0f} dB weaker)\n")
    print(f"{'Window':<14}{'Detections':>11}{'Weak found':>12}{'False alarms':>14}")

    for name, use_window in (("Rectangular", False), ("Hann", True)):
        ranges, velocities, rd = range_doppler_map(p, F, window=use_window)
        rd_db = to_db_above_noise(rd)

        detections = detect_targets(p, rd_db, threshold_db=THRESHOLD_DB)
        matches, false_alarms = match_to_truth(p, detections, TARGETS)
        weak_found = "yes" if matches[1][1] is not None else "no"
        print(f"{name:<14}{len(detections):>11}{weak_found:>12}{len(false_alarms):>14}")

        # Range cut through the objects' velocity
        j = int(np.argmax(np.max(rd_db, axis=0)))
        sel = (ranges > 55) & (ranges < 110)
        ax.plot(ranges[sel], rd_db[sel, j] - rd_db[:, j].max(), label=f"{name} window")

    for t in TARGETS:
        ax.axvline(t.range_m, color="red", linestyle="--", linewidth=1)
    ax.annotate("weak object\n(30 dB weaker)", xy=(TARGETS[1].range_m, -30), xytext=(95, -12),
                arrowprops=dict(arrowstyle="->"), fontsize=9)
    ax.set_ylim(-70, 3)
    ax.set_xlabel("Range [m]")
    ax.set_ylabel("Power relative to strong object [dB]")
    ax.set_title("Sidelobes of a strong object can hide a weak one (red lines = true positions)")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right")
    fig.tight_layout()

    Path("results").mkdir(exist_ok=True)
    fig.savefig("results/windowing.png", dpi=150)
    print("\nSaved results/windowing.png")