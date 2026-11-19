"""Absolute scattering cross section: the ion-chamber calibration.

    I_abs(Q) = coef * I(Q),   coef = N_frames / (t_exp * F * T * L * dOmega)

following Sheyfer et al., Phys. Rev. Lett. 125, 125504 (2020), Eq. S1.  I(Q) is
in photons per pixel per frame, F is the incident flux (photons/s), T the
sample transmission, t_exp = frame_time * N_frames, L the sample thickness and
dOmega = (pixel / distance)^2.  Lengths are in mm, so coef * I(Q) is in mm^-1
and is plotted without further scaling.

F and T come from each file's own upstream (Up_IC) and downstream (Dn_IC) ion
chambers, averaged over the acquisitions the file contains:
    F = CAL_A * (Up_IC / t_exp) + CAL_B
    T = Dn_IC / Up_IC / AIR_TRANSMISSION
CAL_A, CAL_B and AIR_TRANSMISSION come from the eight attenuation points
below; the repository README describes the calibration.  Used by
series_2021_3.py, temperature_series_2021-3.py and calibration_figure.py.
"""

import numpy as np

# --- fixed beamline constants (from the experiment log / conversation notes) ---
DET_DIST_MM      = 7800.0    # sample-to-detector distance, 7.8 m
PIXEL_MM         = 0.076     # detector pixel size, 76 um (same in x and y)
SAMPLE_THICKNESS = 2.0       # sample thickness, mm
NUM_FRAMES       = 100000    # frames per acquisition
FRAME_TIME       = 20e-6     # s
T_EXP            = FRAME_TIME * NUM_FRAMES          # total exposure time [s]
DELTA_OMEGA      = (PIXEL_MM / DET_DIST_MM) ** 2    # solid angle of one pixel

# --- IC -> photon calibration (no sample; only the air path attenuates) ---
UPIC_DARK   = 83.6           # upstream IC dark reading
DNIC_DARK   = 204.4          # downstream IC dark reading
# Upstream IC, downstream IC and PIN-diode photon counts at 8 attenuation ratios
# (dark already subtracted from the photon counts).
CAL_UPIC    = np.array([257609, 125204, 55586.5, 27036.8, 12296.8, 6028.8, 2714.4, 1361.1])
CAL_DNIC    = np.array([223888, 108873, 48413.4, 23613.7, 10807.5, 5367.4, 2485.5, 1310.9])
CAL_PHOTONS = np.array([1.49e10, 8.12e9, 3.75e9, 1.97e9, 1.06e9, 6.71e8, 4.66e8, 3.83e8]) - 3.04e8
CAL_CROP    = 1              # the unattenuated point lies off the line and is not fitted

# photons = CAL_A * Up_IC + CAL_B, fitted to the seven attenuated points
CAL_A, CAL_B = np.polyfit(CAL_UPIC[CAL_CROP:], CAL_PHOTONS[CAL_CROP:], 1)

# transmission of the air path between the chambers: mean and sample standard
# deviation of (Dn_IC - dark) / (Up_IC - dark) over the same seven points
AIR_TRANS_SERIES = (CAL_DNIC[CAL_CROP:] - DNIC_DARK) / (CAL_UPIC[CAL_CROP:] - UPIC_DARK)
AIR_TRANSMISSION = float(np.mean(AIR_TRANS_SERIES))
AIR_TRANS_STD    = float(np.std(AIR_TRANS_SERIES, ddof=1))

# where each averaged file keeps its ion-chamber readings
INCIDENT_PATH    = '/entry/instrument/incident_beam/incident_beam_intensity'
TRANSMITTED_PATH = '/entry/instrument/incident_beam/transmitted_beam_intensity'


def coef_from_monitors(up_ic, dn_ic):
    """coef for upstream and downstream ion-chamber readings given as numbers."""
    flux = CAL_A * (up_ic / T_EXP) + CAL_B
    transmission = dn_ic / up_ic / AIR_TRANSMISSION
    return NUM_FRAMES / (T_EXP * flux * transmission * SAMPLE_THICKNESS * DELTA_OMEGA)


def abs_xsec_coef(hf):
    """coef for one open averaged file, from the monitors stored in it."""
    up_ic = float(hf[INCIDENT_PATH][()])
    dn_ic = float(hf[TRANSMITTED_PATH][()])
    return coef_from_monitors(up_ic, dn_ic)


def calibration_summary():
    """The calibration in two lines, for the log."""
    return (f'IC->photon calibration : photons = {CAL_A:.3e} * Up_IC {CAL_B:+.3e}\n'
            f'air transmission        : {AIR_TRANSMISSION:.4f} +/- {AIR_TRANS_STD:.4f}')
