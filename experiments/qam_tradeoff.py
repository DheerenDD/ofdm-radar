"""Experiment: the trade-off between data rate (communication) and radar quality (sensing).

The SAME frame is used twice:
  1. A communication user decodes it over a noisy link -> data rate and bit error rate.
  2. The base station uses its echo for radar -> how clearly the objects stand out.
Higher-order QAM carries more bits per symbol but makes the radar slightly worse.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc

from isac.detection import detect_targets, match_to_truth
from isac.frame import generate_frame, symbols_to_bits
from isac.params import OFDMRadarParams
from isac.processing import range_doppler_map, remove_data, to_db_above_noise
from isac.targets import DEFAULT_SCENE, simulate_echo

ORDERS = [4, 16, 64]          # QPSK, 16-QAM, 64-QAM
COMMS_SNR_DB = 20             # SNR at the communication user
RADAR_SNR_DB = -20            # per-cell SNR of a 0 dB object at the radar
THRESHOLD_DB = 15


def theoretical_ber(order, snr_db):
    """Approximate BER of Gray-coded square QAM in AWGN (exact for QPSK)."""
    k = np.log2(order)
    es_n0 = 10 ** (snr_db / 10)
    q = 0.5 * erfc(np.sqrt(3 * es_n0 / (order - 1)) / np.sqrt(2))
    return (4 / k) * (1 - 1 / np.sqrt(order)) * q


def peak_near(p, rd_db, target):
    """Peak height (dB above noise) within +/- 2 bins of a true object."""
    i = int(round(target.range_m / p.range_resolution))
    j = int(round(target.velocity_mps / p.velocity_resolution)) + p.M // 2
    return rd_db[max(i - 2, 0):i + 3, max(j - 2, 0):j + 3].max()


if __name__ == "__main__":
    p = OFDMRadarParams()
    rng = np.random.default_rng(seed=6)
    weakest = min(DEFAULT_SCENE, key=lambda t: t.gain_db)
    rows = []

    for order in ORDERS:
        k = int(np.log2(order))
        bits, X = generate_frame(p, order=order, rng=rng)

        # 1) Communication: the user receives the frame over an AWGN link and decodes it
        noise = np.sqrt(10 ** (-COMMS_SNR_DB / 10) / 2) * (rng.standard_normal(X.shape)
                                                          + 1j * rng.standard_normal(X.shape))
        bits_rx = symbols_to_bits((X + noise).reshape(-1, order="F"), order)
        ber = np.mean(bits != bits_rx)
        data_rate = p.N * k / p.T_sym          # bits per second, all subcarriers, no coding

        # 2) Sensing: the base station processes the echo of the same frame
        Y = simulate_echo(p, X, DEFAULT_SCENE, snr_db=RADAR_SNR_DB, rng=rng)
        ranges, velocities, rd = range_doppler_map(p, remove_data(Y, X))
        rd_db = to_db_above_noise(rd)
        matches, false_alarms = match_to_truth(p, detect_targets(p, rd_db, THRESHOLD_DB), DEFAULT_SCENE)
        found = sum(d is not None for _, d in matches)

        # Noise amplification caused by dividing by low-amplitude symbols
        noise_gain_db = 10 * np.log10(np.mean(1 / np.abs(X) ** 2))

        rows.append(dict(name="QPSK" if order == 4 else f"{order}-QAM", k=k, rate=data_rate / 1e6,
                         ber=ber, ber_th=theoretical_ber(order, COMMS_SNR_DB),
                         noise_gain=noise_gain_db, weak_peak=peak_near(p, rd_db, weakest),
                         found=found, false=len(false_alarms)))

    print(f"Communication link SNR {COMMS_SNR_DB} dB | radar per-cell SNR {RADAR_SNR_DB} dB | "
          f"weakest object: {weakest.range_m:.0f} m, {weakest.gain_db:+.0f} dB\n")
    print(f"{'Modulation':<11}{'Bits/sym':>9}{'Data rate':>13}{'BER (sim)':>12}{'BER (theory)':>14}"
          f"{'Radar noise +':>15}{'Weakest peak':>14}{'Found':>8}")
    for r in rows:
        print(f"{r['name']:<11}{r['k']:>9}{r['rate']:>8.0f} Mb/s{r['ber']:>12.1e}{r['ber_th']:>14.1e}"
              f"{r['noise_gain']:>12.1f} dB{r['weak_peak']:>11.1f} dB{r['found']:>5}/{len(DEFAULT_SCENE)}")

    # Plot: data rate goes up, radar quality goes down
    names = [r["name"] for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].bar(names, [r["rate"] for r in rows], color="#1F3864")
    axes[0].set_title("Communication: data rate")
    axes[0].set_ylabel("Data rate [Mbit/s] (uncoded)")
    axes[1].bar(names, [r["weak_peak"] for r in rows], color="#C0504D")
    axes[1].axhline(THRESHOLD_DB, color="black", linestyle="--", linewidth=1, label="detection threshold")
    axes[1].set_title("Sensing: weakest object above noise")
    axes[1].set_ylabel("Peak above noise [dB]")
    axes[1].legend(loc="upper right")
    for ax in axes:
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Same frame, two jobs: more bits per symbol vs. clearer radar")
    fig.tight_layout()

    Path("results").mkdir(exist_ok=True)
    fig.savefig("results/qam_tradeoff.png", dpi=150)
    print("\nSaved results/qam_tradeoff.png")