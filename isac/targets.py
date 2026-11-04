"""Point targets and the monostatic OFDM radar echo model."""
from dataclasses import dataclass
from pathlib import Path
import warnings

import numpy as np
from scipy.constants import c

from isac.params import OFDMRadarParams


@dataclass(frozen=True)
class Target:
    range_m: float                   # distance from the radar [m]
    velocity_mps: float              # radial velocity [m/s], positive = approaching the radar
    gain_db: float = 0.0             # echo power relative to a reference target [dB]
    phase_rad: float | None = None   # echo phase; None = random (realistic)

# Example scene reused by later steps and experiments
DEFAULT_SCENE = [
    Target(range_m=30.0, velocity_mps=-1.5, gain_db=-10.0),   # pedestrian walking away
    Target(range_m=80.0, velocity_mps=20.0, gain_db=0.0),     # car approaching
    Target(range_m=150.0, velocity_mps=-12.0, gain_db=3.0),   # truck driving away
]


def check_targets(params, targets):
    """Warn if a target breaks the assumptions of the simple echo model."""
    for t in targets:
        if t.range_m > params.max_range_cp:
            warnings.warn(f"Target at {t.range_m} m is beyond the CP limit "
                          f"({params.max_range_cp:.0f} m): ISI is ignored by this model.")
        if abs(params.doppler_ratio(t.velocity_mps)) > 0.1:
            warnings.warn(f"Target at {t.velocity_mps} m/s has Doppler above 10% of the "
                          "subcarrier spacing: ICI is ignored by this model.")


def simulate_echo(params, X, targets, snr_db=None, rng=None):
    """Received radar grid Y[n, m] for point targets.

    Each target multiplies the transmitted grid by
      exp(-j 2 pi n delta_f tau)   phase ramp across subcarriers  -> range
      exp(+j 2 pi m T_sym f_D)     phase rotation across symbols  -> velocity
    snr_db is the per-cell SNR of a 0 dB target; None means noise-free.
    """
    rng = np.random.default_rng() if rng is None else rng
    check_targets(params, targets)

    n = np.arange(params.N)[:, None]    # subcarrier index, column vector
    m = np.arange(params.M)[None, :]    # symbol index, row vector

    H = np.zeros((params.N, params.M), dtype=complex)
    for t in targets:
        tau = 2 * t.range_m / c                          # round-trip delay
        f_d = 2 * t.velocity_mps / params.wavelength     # Doppler shift
        phase = rng.uniform(0, 2 * np.pi) if t.phase_rad is None else t.phase_rad
        amplitude = 10 ** (t.gain_db / 20) * np.exp(1j * phase)
        H += (amplitude
              * np.exp(-1j * 2 * np.pi * n * params.delta_f * tau)
              * np.exp(1j * 2 * np.pi * m * params.T_sym * f_d))

    Y = H * X
    if snr_db is not None:
        noise_power = 10 ** (-snr_db / 10)
        Y = Y + np.sqrt(noise_power / 2) * (rng.standard_normal(Y.shape)
                                            + 1j * rng.standard_normal(Y.shape))
    return Y


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from isac.frame import generate_frame

    p = OFDMRadarParams()
    rng = np.random.default_rng(seed=2)
    _, X = generate_frame(p, order=4, rng=rng)

    # --- Check 1: one noise-free target, recover range and speed from the phase steps ---
    truth = Target(range_m=80.0, velocity_mps=20.0)
    F = simulate_echo(p, X, [truth], rng=rng) / X        # divide out the data

    step_n = np.angle(np.mean(F[1:, :] * np.conj(F[:-1, :])))   # phase step per subcarrier
    step_m = np.angle(np.mean(F[:, 1:] * np.conj(F[:, :-1])))   # phase step per symbol
    range_est = -step_n / (2 * np.pi * p.delta_f) * c / 2
    velocity_est = step_m / (2 * np.pi * p.T_sym) * p.wavelength / 2

    print("Check 1: single noise-free target, phase-slope estimate")
    print(f"  Range   : true {truth.range_m:.2f} m,   estimated {range_est:.2f} m")
    print(f"  Velocity: true {truth.velocity_mps:.2f} m/s, estimated {velocity_est:.2f} m/s")

    # --- Check 2: noisy multi-target scene, SNR before processing ---
    snr_db = -10
    Y = simulate_echo(p, X, DEFAULT_SCENE, snr_db=snr_db, rng=rng)
    gain_db = 10 * np.log10(p.N * p.M)
    print("Check 2: default scene with noise")
    print(f"  Grid shape: {Y.shape}, per-cell SNR of a 0 dB target: {snr_db} dB")
    print(f"  Processing gain of the 2D FFT (N*M): {gain_db:.1f} dB "
          f"-> expected peak SNR about {snr_db + gain_db:.1f} dB")

    # --- Plot: the phase patterns that encode range and velocity ---
    Path("results").mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    tau = 2 * truth.range_m / c
    f_d = 2 * truth.velocity_mps / p.wavelength
    n = np.arange(p.N)
    m = np.arange(p.M)

    phase_n = np.unwrap(np.angle(F[:, 0]))
    axes[0].plot(n, phase_n - phase_n[0], label="simulated")
    axes[0].plot(n, -2 * np.pi * n * p.delta_f * tau, "--", label="theory")
    axes[0].set_title(f"Phase across subcarriers (target at {truth.range_m:.0f} m)")
    axes[0].set_xlabel("Subcarrier index n")
    axes[0].set_ylabel("Unwrapped phase [rad]")

    phase_m = np.unwrap(np.angle(F[0, :]))
    axes[1].plot(m, phase_m - phase_m[0], label="simulated")
    axes[1].plot(m, 2 * np.pi * m * p.T_sym * f_d, "--", label="theory")
    axes[1].set_title(f"Phase across symbols (target at {truth.velocity_mps:.0f} m/s)")
    axes[1].set_xlabel("OFDM symbol index m")
    axes[1].set_ylabel("Unwrapped phase [rad]")

    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend()
    fig.tight_layout()
    fig.savefig("results/echo_phase.png", dpi=150)
    print("Saved results/echo_phase.png")