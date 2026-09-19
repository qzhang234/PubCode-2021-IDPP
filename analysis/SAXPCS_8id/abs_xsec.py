"""Absolute scattering cross-section calibration (shared helper).

saxpcs.py (Figure 3 and Figure S3) and saxs_evolution.py (Figure S8) both
import from here, so they use IDENTICAL constants, the same IC->photon
calibration and the same per-file coefficient function -- there is only one
place to change if the calibration is revised.

PROVENANCE: this calibration was originally worked out in two beamline Jupyter
notebooks, IC_pind4_convert.ipynb (steps 1 and 2) and Abs_Scatt_Cross.ipynb
(step 3).  Those notebooks are NOT part of this repository; the code below is
the single authoritative implementation and reproduces their numbers.

GOAL
----
Turn the raw azimuthally-averaged SAXS intensity I(Q) (photons per frame, per
pixel) into an absolute differential scattering cross section dSigma/dOmega,
with the buffer subtracted:

    I_abs(Q) = coef_sam * I_sample(Q) - coef_buf * I_buffer(Q)

Theory: Sheyfer et al., Phys. Rev. Lett. 125, 125504 (2020), Eq. S1

    I(Q) = (F_s  T_s  t_exp f dOmega)^-1 I_s(Q)
         - (F_bg T_bg t_exp f dOmega)^-1 I_bg(Q)

For one measurement:
    F      = incident X-ray flux            [photons / s]
    T      = sample (or buffer) transmission [dimensionless]
    t_exp  = total exposure time            [s]  (= frame_time * num_frames)
    f      = sample thickness
    dOmega = solid angle of one detector pixel = (pixel / distance)^2

giving the per-measurement coefficient

    coef = num_frames / (t_exp * F * T * f * dOmega)

F and T are not measured directly; they come from the upstream (Up_IC) and
downstream (Dn_IC) ion chambers stored per run in
    /entry/instrument/incident_beam/incident_beam_intensity     (Up_IC)
    /entry/instrument/incident_beam/transmitted_beam_intensity  (Dn_IC)
These were copied from the raw cluster results into the NeXus files and are
averaged over each frame range by average_ranges.py -- so every averaged file
carries the mean Up_IC / Dn_IC of the frames that went into it.

STEP 1 -- ion chamber -> photon flux
------------------------------------
With NO sample in the beam, Up_IC was recorded at eight X-ray attenuation
ratios together with the photon counts from the downstream PIN diode (pind4,
converted to photons with the standard 8-ID-I pind4 calibration).  A straight
line  photons = CAL_A * Up_IC + CAL_B  is fit to seven of the eight: the
UNATTENUATED point is dropped (CAL_CROP=1) because it is off the line the other
seven define -- keeping it leaves residuals of -2 to -148 % across the set,
dropping it leaves +-1 % over the middle five.  The flux of a real run then follows the
original convention  F = CAL_A * (Up_IC / t_exp) + CAL_B.

STEP 2 -- air transmission
--------------------------
With no sample, (Dn_IC - dark) / (Up_IC - dark) is the transmission of the air
path between the two chambers; it is flat at ~0.868.  The sample/buffer
transmission is therefore  T = Dn_IC / Up_IC / AIR_TRANSMISSION.

STEP 3 -- coefficient
---------------------
abs_xsec_coef() evaluates coef for one open file from its range-averaged
Up_IC / Dn_IC and the constants below.  It is called once per SAXS curve, so
each measurement is scaled by ITS OWN flux and transmission (the reason
average_ranges.py averages the two monitors per range); coef_buf uses the
D0138 buffer's own Up_IC / Dn_IC.

UNITS: lengths follow the original convention (millimetres for pixel, distance
and thickness); dOmega is a pure ratio, so its length units cancel.  The
constants below are the ones in the beamline log for this experiment.

THE SAMPLE THICKNESS IS 2 mm.  That is the cell used in this experiment, it is
what SAMPLE_THICKNESS is set to below, and it is what every number in the paper
is built on: coef_buf = 3.75e4 for the D0138 buffer.

The only reason 1 mm is mentioned anywhere is the cross-check against the
original beamline notebook, which had been written for a 1 mm cell.  Setting
SAMPLE_THICKNESS = 1 here reproduces that notebook to 2 % (7.50e4 here against
7.62e4 there), which is how this implementation was confirmed to be the same
calculation.  Thickness enters only as 1/f, so 1 mm and 2 mm differ by exactly a
factor of two and nothing else.  The 1 mm value is NOT used by this analysis.

    [coef] = 1 / ([T_EXP] [F] [T] [f] [dOmega])
           = 1 / ( s * photons/s * 1 * mm * sr )
           = 1 / (photons mm sr)

so coef * I(Q) comes out in mm^-1 sr^-1.  Every absolute-scale figure in this
work reports that number as it stands: the plotted quantity is exactly
coef * I(Q), with no further scaling.  Absolute cross sections are more often
quoted in cm^-1; SAMPLE_THICKNESS is the only dimensional length that survives
(dOmega's mm cancel), so a reader comparing with the literature multiplies by
ten.  That factor is deliberately NOT applied anywhere in this repository, so
that what is plotted is the coefficient the calibration actually produces.
"""

import numpy as np


# --- fixed beamline constants (from the experiment log / conversation notes) ---
DET_DIST_MM      = 7800.0    # sample-to-detector distance, 7.8 m
PIXEL_MM         = 0.076     # detector pixel size, 76 um (same in x and y)
SAMPLE_THICKNESS = 2.0       # sample / capillary thickness, 2 mm
NUM_FRAMES       = 100000    # frames per measurement
FRAME_TIME       = 20e-6     # acquisition (frame) time, 20 us
T_EXP            = FRAME_TIME * NUM_FRAMES          # total exposure time [s]
DELTA_OMEGA      = (PIXEL_MM / DET_DIST_MM) ** 2    # pixel solid angle

# --- IC -> photon calibration data (no sample; only the air gap attenuates) ---
UPIC_DARK   = 83.6           # upstream IC dark reading
DNIC_DARK   = 204.4          # downstream IC dark reading
# Upstream IC, downstream IC and PIN-diode photon counts at 8 attenuation ratios
# (dark already subtracted from the photon counts).
CAL_UPIC    = np.array([257609, 125204, 55586.5, 27036.8, 12296.8, 6028.8, 2714.4, 1361.1])
CAL_DNIC    = np.array([223888, 108873, 48413.4, 23613.7, 10807.5, 5367.4, 2485.5, 1310.9])
CAL_PHOTONS = np.array([1.49e10, 8.12e9, 3.75e9, 1.97e9, 1.06e9, 6.71e8, 4.66e8, 3.83e8]) - 3.04e8
CAL_CROP    = 1              # drop the first (unattenuated) point before fitting

# STEP 1: linear fit  photons = CAL_A * Up_IC + CAL_B
CAL_A, CAL_B = np.polyfit(CAL_UPIC[CAL_CROP:], CAL_PHOTONS[CAL_CROP:], 1)

# STEP 2: air transmission = mean of (Dn_IC - dark) / (Up_IC - dark).  The
# spread of the per-measurement values (AIR_TRANS_STD) is a fair error-from-mean
# estimate and is reported in the calibration figure legend.
AIR_TRANS_SERIES = (CAL_DNIC[CAL_CROP:] - DNIC_DARK) / (CAL_UPIC[CAL_CROP:] - UPIC_DARK)
AIR_TRANSMISSION = float(np.mean(AIR_TRANS_SERIES))
AIR_TRANS_STD    = float(np.std(AIR_TRANS_SERIES, ddof=1))   # sample standard deviation

# HDF fields holding the (range-averaged) ion-chamber monitors.
INCIDENT_PATH    = '/entry/instrument/incident_beam/incident_beam_intensity'
TRANSMITTED_PATH = '/entry/instrument/incident_beam/transmitted_beam_intensity'


def abs_xsec_coef(hf):
    """Absolute cross-section coefficient for one open XPCS/SAXS file.

    Reads the range-averaged upstream (Up_IC) and downstream (Dn_IC)
    ion-chamber intensities and combines them with the beamline constants and
    the IC->photon calibration (see the module docstring):

        F    = CAL_A * (Up_IC / T_EXP) + CAL_B          incident flux [photons/s]
        T    = Dn_IC / Up_IC / AIR_TRANSMISSION         transmission  [-]
        coef = NUM_FRAMES / (T_EXP * F * T * SAMPLE_THICKNESS * DELTA_OMEGA)

    Returns the scalar coefficient that multiplies this file's I(Q).
    """
    # These are stored as tiny arrays whose shape varies between files -- (1,)
    # in some, (1, 1) in others -- so flatten to 1-D and take the first entry.
    up_ic = float(np.asarray(hf[INCIDENT_PATH][()]).ravel()[0])
    dn_ic = float(np.asarray(hf[TRANSMITTED_PATH][()]).ravel()[0])
    flux = CAL_A * (up_ic / T_EXP) + CAL_B
    transmission = dn_ic / up_ic / AIR_TRANSMISSION
    return NUM_FRAMES / (T_EXP * flux * transmission * SAMPLE_THICKNESS * DELTA_OMEGA)


def calibration_summary():
    """One-line-per-item summary string for logging at import in the scripts."""
    return (f'IC->photon calibration : photons = {CAL_A:.3e} * Up_IC {CAL_B:+.3e}\n'
            f'air transmission        : {AIR_TRANSMISSION:.4f} +/- {AIR_TRANS_STD:.4f}')
