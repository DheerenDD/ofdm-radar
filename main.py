"""Run the whole project: every step and experiment, in order, saving all figures to results/."""
import runpy

import matplotlib

matplotlib.use("Agg")      # save figures to files without opening windows

STEPS = [
    ("1. System parameters and radar limits", "isac.params"),
    ("2. Transmit frame (Gray-coded QAM)", "isac.frame"),
    ("3. Echo model validation", "isac.targets"),
    ("4. Range-Doppler map", "isac.processing"),
    ("5. Automatic detection and validation", "isac.detection"),
    ("6. Experiment: bandwidth vs resolution", "experiments.resolution_vs_bandwidth"),
    ("7. Experiment: windowing", "experiments.windowing"),
    ("8. Experiment: QPSK vs 16-QAM vs 64-QAM", "experiments.qam_tradeoff"),
]

if __name__ == "__main__":
    for title, module in STEPS:
        print("\n" + "=" * 70)
        print(title)
        print("=" * 70)
        runpy.run_module(module, run_name="__main__")
    print("\nAll steps finished. Figures are in the results/ folder.")