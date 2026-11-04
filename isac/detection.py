"""Automatic target detection on the range-Doppler map and validation against ground truth."""
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.ndimage import maximum_filter

from isac.params import OFDMRadarParams


@dataclass(frozen=True)
class Detection:
    range_m: float
    velocity_mps: float
    snr_db: float           # peak height above the noise floor


def _refine(left, centre, right):
    """Sub-bin offset of a peak from three neighbouring dB values (parabolic interpolation)."""
    denom = left - 2 * centre + right
    return 0.0 if denom == 0 else 0.5 * (left - right) / denom


def detect_targets(params, rd_db, threshold_db=15.0, neighbourhood=(5, 5)):
    """Find peaks that are local maxima AND above the threshold (dB above noise).

    Each peak position is refined between FFT bins so the estimate is finer
    than the bin spacing.
    """
    is_peak = (rd_db == maximum_filter(rd_db, size=neighbourhood)) & (rd_db > threshold_db)
    detections = []
    for i, j in zip(*np.nonzero(is_peak)):
        di = _refine(rd_db[i - 1, j], rd_db[i, j], rd_db[i + 1, j]) if 0 < i < params.N - 1 else 0.0
        dj = _refine(rd_db[i, j - 1], rd_db[i, j], rd_db[i, j + 1]) if 0 < j < params.M - 1 else 0.0
        detections.append(Detection(
            range_m=(i + di) * params.range_resolution,
            velocity_mps=(j + dj - params.M // 2) * params.velocity_resolution,
            snr_db=float(rd_db[i, j]),
        ))
    return sorted(detections, key=lambda d: d.range_m)


def match_to_truth(params, detections, targets):
    """Pair each true target with the nearest detection (distance measured in resolution cells)."""
    results = []
    unused = list(detections)
    for t in targets:
        if not unused:
            results.append((t, None))
            continue
        cell_dist = [np.hypot((d.range_m - t.range_m) / params.range_resolution,
                              (d.velocity_mps - t.velocity_mps) / params.velocity_resolution)
                     for d in unused]
        best = int(np.argmin(cell_dist))
        if cell_dist[best] <= 2:                 # within 2 resolution cells counts as a hit
            results.append((t, unused.pop(best)))
        else:
            results.append((t, None))
    return results, unused                       # unused detections are false alarms


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from isac.frame import generate_frame
    from isac.processing import plot_range_doppler, range_doppler_map, remove_data, to_db_above_noise
    from isac.targets import DEFAULT_SCENE, simulate_echo

    p = OFDMRadarParams()
    rng = np.random.default_rng(seed=3)

    _, X = generate_frame(p, order=4, rng=rng)
    Y = simulate_echo(p, X, DEFAULT_SCENE, snr_db=-10, rng=rng)
    ranges, velocities, rd = range_doppler_map(p, remove_data(Y, X))
    rd_db = to_db_above_noise(rd)

    detections = detect_targets(p, rd_db, threshold_db=15.0)
    matches, false_alarms = match_to_truth(p, detections, DEFAULT_SCENE)

    print(f"Detections above 15 dB: {len(detections)}  (true targets: {len(DEFAULT_SCENE)}, "
          f"false alarms: {len(false_alarms)})\n")
    header = (f"{'True range':>11} {'Est. range':>11} {'Error':>8}   "
              f"{'True speed':>11} {'Est. speed':>11} {'Error':>9}   {'Peak':>7}")
    print(header)
    print("-" * len(header))
    range_err, vel_err = [], []
    for t, d in matches:
        if d is None:
            print(f"{t.range_m:>9.2f} m {'missed':>11}")
            continue
        er, ev = d.range_m - t.range_m, d.velocity_mps - t.velocity_mps
        range_err.append(er)
        vel_err.append(ev)
        print(f"{t.range_m:>9.2f} m {d.range_m:>9.2f} m {er:>+6.2f} m   "
              f"{t.velocity_mps:>+7.2f} m/s {d.velocity_mps:>+7.2f} m/s {ev:>+6.2f} m/s   "
              f"{d.snr_db:>5.1f} dB")

    print(f"\nRMS error: range {np.sqrt(np.mean(np.square(range_err))):.2f} m "
          f"(resolution {p.range_resolution:.2f} m), "
          f"velocity {np.sqrt(np.mean(np.square(vel_err))):.2f} m/s "
          f"(resolution {p.velocity_resolution:.2f} m/s)")

    Path("results").mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_range_doppler(ranges, velocities, rd_db, targets=DEFAULT_SCENE, ax=ax,
                       title="Automatic detection (threshold 15 dB above noise)")
    ax.scatter([d.velocity_mps for d in detections], [d.range_m for d in detections],
               marker="x", s=80, color="white", linewidths=2, label="detected")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig("results/detections.png", dpi=150)
    print("Saved results/detections.png")