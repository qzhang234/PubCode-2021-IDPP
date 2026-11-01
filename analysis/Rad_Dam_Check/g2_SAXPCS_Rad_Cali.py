"""SAXS + XPCS for diluted H06 at four attenuation conditions, with the legend
labelled by the actual photon flux on the sample (not the nominal attenuation).

The nominal attenuator settings (a7, a9, a11, a14) do NOT give a monotonic flux
because the attenuator setup was wrong.  The true on-sample flux is recovered
from the upstream ion-chamber (IC) reading of each measurement:

    IC readout  --(Pind_calibration.csv)-->  calibrated PIN-diode current
                --(Flux_Cal.py formula)---->  photon flux on sample [ph/s]

The IC readout is stored in each results file at
``entry/instrument/incident_beam/incident_beam_intensity`` and must be scaled by
1e-9 (matches the numbers reported on Slide_Screenshot.png).

Source data: the four result files read here, cluster_results/E0171-E0174,
are committed, byte-identical copies of the copy of record on beamline
storage, /gdata/s8id-dmdtn/2025-1/qzhang202503/flux_control_results (placed
there on 1 Nov 2026). They were correlated in 2025 from raw frames under
/gdata/dm/8ID/8IDI/2025-1/qzhang202503/data with the qmap
rigaku_s360_d36_lin_MC.hdf, and each file records both paths. That DM
experiment was archived on 23 Aug 2025. The same-named files in
.../qzhang202503/Multitau are an older reduction with a different qmap and are
not used.
"""

import os
import sys

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.acs_style import (DOUBLE_COL, MS, MEW, LW_THIN,
                              apply_style, add_minor_grid, label_panels,
                              q_log_ticks, save_fig)

# --- PATHS (all local, relative to this script) ---
# cluster_results/ is the committed copy of the source data (see the docstring)
here = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(here, 'cluster_results')
calib_csv = os.path.join(here, 'Pind_calibration.csv')

# nominal attenuation label -> results file (diluted H06, c6c5-4)
FILES = {
    'a7':  'E0174_H06-c6c5-4_a0007_f100000_r00001_results.hdf',
    'a9':  'E0173_H06-c6c5-4_a0009_f100000_r00001_results.hdf',
    'a11': 'E0172_H06-c6c5-4_a0011_f100000_r00001_results.hdf',
    'a14': 'E0171_H06-c6c5-4_a0014_f100000_r00001_results.hdf',
}
ORDER = ['a7', 'a9', 'a11', 'a14']       # legend / colour order
COLORS = {'a7': 'C0', 'a9': 'C1', 'a11': 'C2', 'a14': 'C3'}

IC_PATH  = 'entry/instrument/incident_beam/incident_beam_intensity'
IC_SCALE = 1e-9                          # apply to the stored IC readout

# --- FLUX CONVERSION CONSTANTS (from Flux_Cal.py) ---
# calibrated PIN-diode current [A]  ->  photon flux on sample [ph/s].
# NOTE: we deliberately drop the attenuation-ratio factor used in Flux_Cal.py
# (there it recovers the *unattenuated* full-beam flux); here the IC/PIN reading
# already corresponds to the attenuated flux hitting the sample.
SPEC_RES      = 0.1882
PHOTON_ENERGY = 12.4e3                   # eV
EV2J          = 1.6e-19                  # electron charge [C]; converts eV -> J
LOSS_MIRROR   = 0.94 * 0.94             # mirror reflectivity
LOSS_DIAMOND  = 0.735                    # BPM + CVD diamond window at 8-ID-E
LOSS_FACTOR   = LOSS_MIRROR * LOSS_DIAMOND


def pind_to_flux(pind_current):
    """Calibrated PIN-diode current [A] -> corrected photon flux [ph/s]."""
    beam_power = pind_current / SPEC_RES
    raw_flux = beam_power / (PHOTON_ENERGY * EV2J)
    return raw_flux / LOSS_FACTOR


def ic_to_pind_slope(csv_path):
    """Fit calibrated-PIN vs up-IC (through the origin) from the calibration CSV.

    Only rows with a real positive IC reading and a numeric (non-saturated)
    calibrated-PIN value are used.  Both detectors respond linearly to flux and
    read ~0 with the shutter closed, so a single proportionality constant maps
    one to the other.
    """
    df = pd.read_csv(csv_path)
    ic = pd.to_numeric(df['Up IC (120 nA, 1e9 scale)'], errors='coerce')
    pind = pd.to_numeric(df['Calibrated pind'], errors='coerce')
    good = ic.notna() & pind.notna() & (ic > 0) & (pind > 0)
    x, y = ic[good].values, pind[good].values
    slope = np.sum(x * y) / np.sum(x * x)          # least squares through origin
    return slope


def fmt_flux(f):
    """Format a flux as '$m.m\\times10^{e}$ ph/s' for the legend."""
    exp = int(np.floor(np.log10(f)))
    mant = f / 10 ** exp
    return rf'${mant:.1f}\times10^{{{exp}}}$ ph/s'


# --- BUILD IC -> flux calibration ---
slope = ic_to_pind_slope(calib_csv)
print(f'IC -> calibrated-PIN slope (through origin): {slope:.4e}')

# --- READ DATA + COMPUTE FLUX ---
data = {}
print('\ncondition   IC [A]      cal-PIN [A]   flux [ph/s]')
for key in ORDER:
    with h5py.File(os.path.join(data_dir, FILES[key]), 'r') as hf:
        ic_raw = hf[IC_PATH][()]
        ic_raw = float(np.asarray(ic_raw).reshape(-1)[0])
        ic = ic_raw * IC_SCALE
        g2 = hf['xpcs/multitau/normalized_g2'][:]
        g2_err = hf['xpcs/multitau/normalized_g2_err'][:]
        saxs_1d = np.squeeze(hf['xpcs/temporal_mean/scattering_1d'][:])
        t_el = hf['xpcs/multitau/delay_list'][:]
        t0 = hf['entry/instrument/detector_1/frame_time'][()]
        ql_sta = hf['xpcs/qmap/static_v_list_dim0'][:]
        ql_dyn = hf['xpcs/qmap/dynamic_v_list_dim0'][:]

    pind = slope * ic
    flux = pind_to_flux(pind)
    data[key] = dict(ic=ic, flux=flux, g2=g2, g2_err=g2_err, saxs_1d=saxs_1d,
                     # static_v_list_dim0 is already in Ainv; the rest of the
                     # paper plots Q in Ainv, so no nm^-1 conversion here
                     q_saxs=ql_sta, q_g2=ql_dyn[0], t_g2=t_el * t0)
    print(f'{key:>7}   {ic:.3e}   {pind:.3e}    {flux:.3e}')

# --- NUMBERS QUOTED IN THE SI CAPTION AND THE TEXT ---
MAIN_FLUX = 1.2e10   # ph/s on the main-experiment sample (SI, SA-XPCS at 8-ID-I)
by_flux = sorted(ORDER, key=lambda k: data[k]['flux'])
low = data[by_flux[0]]
print('\nordered by flux:  ' + ''.join(f'{k:>7}' for k in by_flux))
print('monitor ratio:    ' + ''.join(f"{data[k]['ic'] / low['ic']:7.2f}" for k in by_flux))
print('mean I(Q) ratio:  ' + ''.join(
    f"{np.nanmean(data[k]['saxs_1d']) / np.nanmean(low['saxs_1d']):7.2f}" for k in by_flux))
print('flux / main expt: ' + ''.join(f"{data[k]['flux'] / MAIN_FLUX:7.1f}" for k in by_flux))
q_shown = np.concatenate([d['q_saxs'][d['saxs_1d'] > 0] for d in data.values()])
print(f'Q range of panel (a): {q_shown.min():.4f} to {q_shown.max():.3f} A^-1')
print(f"Q of panel (b): {low['q_g2']:.5f} A^-1, "
      f"2pi/Q = {2 * np.pi / low['q_g2'] / 10:.0f} nm")

# --- FIGURE ---
apply_style()
fig, axs = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.9))
label_panels(axs)

# marker shape as well as colour, so the four flux conditions stay separable
# without colour (ACS: "avoid relying on color alone")
MARKERS = {'a7': 'o', 'a9': 's', 'a11': '^', 'a14': 'D'}

for key in ORDER:
    d = data[key]
    color = COLORS[key]
    label = fmt_flux(d['flux'])
    # left: SAXS I(q)
    axs[0].plot(d['q_saxs'], d['saxs_1d'], marker=MARKERS[key], mfc='none', mec=color,
                mew=MEW, ms=MS, color=color, ls='none', label=label)
    # right: g2 at the first q bin
    axs[1].errorbar(d['t_g2'], d['g2'][:, 0], d['g2_err'][:, 0], marker=MARKERS[key],
                    mfc='none', mec=color, mew=MEW, ms=MS, capsize=1.5,
                    elinewidth=LW_THIN, capthick=LW_THIN,
                    color=color, ls='none', label=label)

axs[0].set_xlabel(r'$Q$ ($\AA^{-1}$)')
axs[0].set_ylabel('Intensity (photon/pixel/frame)')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
q_log_ticks(axs[0])   # after set_xscale: setting a scale resets the locators
add_minor_grid(axs[0])
# I(Q) falls from the top left, so open a decade of headroom for the key to
# sit in rather than letting it cover the high-Q end of the curves
_i_hi = max(d['saxs_1d'].max() for d in data.values())
_i_lo = min(d['saxs_1d'][d['saxs_1d'] > 0].min() for d in data.values())
# The top of this axis used to sit 12x above the data purely to clear a legend
# in the upper right.  The curves run diagonally, so the lower-left corner is
# empty whatever the limits are; putting the key there lets the axis close up
# on the data and stretches the vertical scale by about a quarter, which is
# what makes the separation between the four fluxes readable.  No data moved.
axs[0].set_ylim(_i_lo * 0.7, _i_hi * 1.6)
axs[0].legend(title='Flux on sample', loc='lower left')

axs[1].set_xlabel(r'Delay Time, $\Delta t$ (s)')
axs[1].set_ylabel(r'$g_2$')
axs[1].set_xscale('log')
add_minor_grid(axs[1])
axs[1].legend(title='Flux on sample')

fig.tight_layout(pad=0.4, w_pad=1.4)
save_fig(fig, os.path.join(here, 'FigureS6_Flux_Control.pdf'))
plt.show()
