"""System parameters and derived radar performance for the OFDM ISAC simulation."""
from dataclasses import dataclass
from scipy.constants import c


@dataclass(frozen=True)
class OFDMRadarParams:
    fc: float = 10e9          # carrier frequency [Hz], FR3 band
    N: int = 1024             # number of subcarriers
    delta_f: float = 120e3    # subcarrier spacing [Hz], 5G NR numerology mu = 3
    M: int = 256              # OFDM symbols per frame
    cp_ratio: float = 0.25    # cyclic prefix length as a fraction of the useful symbol

    # ---------- basic quantities ----------
    @property
    def wavelength(self) -> float:
        return c / self.fc

    @property
    def bandwidth(self) -> float:
        return self.N * self.delta_f

    @property
    def T_useful(self) -> float:
        return 1 / self.delta_f

    @property
    def T_cp(self) -> float:
        return self.cp_ratio * self.T_useful

    @property
    def T_sym(self) -> float:
        return self.T_useful + self.T_cp

    @property
    def frame_duration(self) -> float:
        return self.M * self.T_sym

    # ---------- radar performance ----------
    @property
    def range_resolution(self) -> float:
        return c / (2 * self.bandwidth)

    @property
    def max_range(self) -> float:
        return c / (2 * self.delta_f)

    @property
    def max_range_cp(self) -> float:
        return c * self.T_cp / 2

    @property
    def velocity_resolution(self) -> float:
        return self.wavelength / (2 * self.M * self.T_sym)

    @property
    def max_velocity(self) -> float:
        return self.wavelength / (4 * self.T_sym)

    def doppler_ratio(self, v: float) -> float:
        """Doppler shift of a target at speed v, as a fraction of the subcarrier spacing."""
        return (2 * v / self.wavelength) / self.delta_f

    def summary(self, v_check: float = 50.0) -> None:
        print("OFDM ISAC radar parameters")
        print(f"  Carrier frequency      : {self.fc / 1e9:.2f} GHz (lambda = {self.wavelength * 100:.2f} cm)")
        print(f"  Subcarriers N          : {self.N}")
        print(f"  Subcarrier spacing     : {self.delta_f / 1e3:.0f} kHz")
        print(f"  Bandwidth              : {self.bandwidth / 1e6:.2f} MHz")
        print(f"  Symbol duration (w/ CP): {self.T_sym * 1e6:.2f} us")
        print(f"  Symbols per frame M    : {self.M}")
        print(f"  Frame duration         : {self.frame_duration * 1e3:.3f} ms")
        print("Radar performance")
        print(f"  Range resolution       : {self.range_resolution:.2f} m")
        print(f"  Max unambiguous range  : {self.max_range:.1f} m")
        print(f"  Max range within CP    : {self.max_range_cp:.1f} m")
        print(f"  Velocity resolution    : {self.velocity_resolution:.2f} m/s")
        print(f"  Max unambiguous speed  : +/-{self.max_velocity:.1f} m/s")
        print(f"  Doppler/delta_f at {v_check:.0f} m/s: {self.doppler_ratio(v_check) * 100:.1f} %")


if __name__ == "__main__":
    OFDMRadarParams().summary()