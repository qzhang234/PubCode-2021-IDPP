"""SI Fig. S3 (file FigureS6_Flux_Control.pdf): SAXS and g2 of a 20 mg/mL
(VPAVG)30 sample at four attenuator settings, keyed by the photon flux on the
sample.  Prints the numbers the SI caption and the text quote.

The nominal attenuator settings (a7, a9, a11, a14) do not give a monotonic
flux, so the flux of each measurement is recovered from its own upstream
ion-chamber (IC) reading:

    IC reading x 1e-9  --(Pind_calibration.csv, fit through the origin)-->  PIN current [A]
                       --(the Flux_Cal.py formula)----------------------->  flux [ph/s]

Reads the four committed result files in cluster_results/; their provenance is
in README.md of this folder.
"""

import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MEW, LW_THIN,               # noqa: E402
                              apply_style, add_minor_grid, label_panels,
                              q_log_ticks, save_fig)

DATA_DIR = os.path.join(HERE, 'cluster_results')
CALIB_CSV = os.path.join(HERE, 'Pind_calibration.csv')

# nominal attenuation label -> result file
FILES = {
    'a7':  'E0174_H06-c6c5-4_a0007_f100000_r00001_results.hdf',
    'a9':  'E0173_H06-c6c5-4_a0009_f100000_r00001_results.hdf',
    'a11': 'E0172_H06-c6c5-4_a0011_f100000_r00001_results.hdf',
    'a14': 'E0171_H06-c6c5-4_a0014_f100000_r00001_results.hdf',
}
ORDER = ['a7', 'a9', 'a11', 'a14']       # legend order
COLORS = {'a7': 'C0', 'a9': 'C1', 'a11': 'C2', 'a14': 'C3'}
MARKERS = {'a7': 'o', 'a9': 's', 'a11': '^', 'a14': 'D'}   # shape as well as colour

IC_PATH = 'entry/instrument/incident_beam/incident_beam_intensity'
IC_SCALE = 1e-9                          # stored IC reading -> A

# PIN current [A] -> flux on the sample [ph/s]: the formula and constants of
# Flux_Cal.py.  The IC reading already corresponds to the attenuated beam.
SPEC_RES      = 0.1882                   # Si PIN-diode spectral responsivity [A/W]
PHOTON_ENERGY = 12.4e3                   # eV
EV2J          = 1.6e-19                  # electron charge [C]; converts eV -> J
LOSS_MIRROR   = 0.94 * 0.94              # mirror reflectivity
LOSS_DIAMOND  = 0.735                    # BPM + CVD diamond window at 8-ID-E
LOSS_FACTOR   = LOSS_MIRROR * LOSS_DIAMOND

MAIN_FLUX = 1.2e10                       # ph/s on the sample in the main experiment


def pind_to_flux(pind_current):
    """PIN-diode current [A] -> photon flux on the sample [ph/s]."""
    beam_power = pind_current / SPEC_RES
    raw_flux = beam_power / (PHOTON_ENERGY * EV2J)
    return raw_flux / LOSS_FACTOR


def ic_to_pind_slope(csv_path):
    """Slope of calibrated PIN current against upstream IC, fitted through the
    origin over the rows where both are positive numbers (not saturated)."""
    df = pd.read_csv(csv_path)
    ic = pd.to_numeric(df['Up IC (120 nA, 1e9 scale)'], errors='coerce')
    pind = pd.to_numeric(df['Calibrated pind'], errors='coerce')
    good = ic.notna() & pind.notna() & (ic > 0) & (pind > 0)
    x = ic[good].values
    y = pind[good].values
    return np.sum(x * y) / np.sum(x * x)          # least squares through the origin


def flux_label(f):
    """A flux as '$m.m\\times10^{e}$ ph/s' for the legend."""
    exponent = int(np.floor(np.log10(f)))
    mantissa = f / 10 ** exponent
    return rf'${mantissa:.1f}\times10^{{{exponent}}}$ ph/s'


slope = ic_to_pind_slope(CALIB_CSV)
print(f'IC -> calibrated-PIN slope (through origin): {slope:.4e}')

data = {}
print('\ncondition   IC [A]      cal-PIN [A]   flux [ph/s]')
for key in ORDER:
    with h5py.File(os.path.join(DATA_DIR, FILES[key]), 'r') as hf:
        ic = float(np.asarray(hf[IC_PATH][()]).reshape(-1)[0]) * IC_SCALE
        g2 = hf['xpcs/multitau/normalized_g2'][:]
        g2_err = hf['xpcs/multitau/normalized_g2_err'][:]
        saxs_1d = np.squeeze(hf['xpcs/temporal_mean/scattering_1d'][:])
        delays = hf['xpcs/multitau/delay_list'][:]
        frame_time = hf['entry/instrument/detector_1/frame_time'][()]
        q_static = hf['xpcs/qmap/static_v_list_dim0'][:]        # A^-1
        q_dynamic = hf['xpcs/qmap/dynamic_v_list_dim0'][:]
    pind = slope * ic
    flux = pind_to_flux(pind)
    data[key] = {'ic': ic, 'flux': flux, 'g2': g2, 'g2_err': g2_err, 'saxs_1d': saxs_1d,
                 'q_saxs': q_static, 'q_g2': q_dynamic[0], 't_g2': delays * frame_time}
    print(f'{key:>7}   {ic:.3e}   {pind:.3e}    {flux:.3e}')

# --- the numbers quoted in the SI caption and the text ---
fluxes = [data[key]['flux'] for key in ORDER]
by_flux = []
for k in np.argsort(fluxes, kind='stable'):
    by_flux.append(ORDER[k])
low = data[by_flux[0]]
line_keys = 'ordered by flux:  '
line_monitor = 'monitor ratio:    '
line_mean = 'mean I(Q) ratio:  '
line_flux = 'flux / main expt: '
for key in by_flux:
    d = data[key]
    line_keys += f'{key:>7}'
    line_monitor += f"{d['ic'] / low['ic']:7.2f}"
    line_mean += f"{np.nanmean(d['saxs_1d']) / np.nanmean(low['saxs_1d']):7.2f}"
    line_flux += f"{d['flux'] / MAIN_FLUX:7.1f}"
print('\n' + line_keys)
print(line_monitor)
print(line_mean)
print(line_flux)
q_shown = []
for key in ORDER:
    d = data[key]
    q_shown.append(d['q_saxs'][d['saxs_1d'] > 0])
q_shown = np.concatenate(q_shown)
print(f'Q range of panel (a): {q_shown.min():.4f} to {q_shown.max():.3f} A^-1')
print(f"Q of panel (b): {low['q_g2']:.5f} A^-1, "
      f"2pi/Q = {2 * np.pi / low['q_g2'] / 10:.0f} nm")

# --- figure ---
apply_style()
fig, axs = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.9))
label_panels(axs)

for key in ORDER:
    d = data[key]
    color = COLORS[key]
    label = flux_label(d['flux'])
    # (a) SAXS I(Q)
    axs[0].plot(d['q_saxs'], d['saxs_1d'], marker=MARKERS[key], mfc='none', mec=color,
                mew=MEW, ms=MS, color=color, ls='none', label=label)
    # (b) g2 at the lowest Q bin
    axs[1].errorbar(d['t_g2'], d['g2'][:, 0], d['g2_err'][:, 0], marker=MARKERS[key],
                    mfc='none', mec=color, mew=MEW, ms=MS, capsize=1.5,
                    elinewidth=LW_THIN, capthick=LW_THIN,
                    color=color, ls='none', label=label)

axs[0].set_xlabel(r'$Q$ ($\AA^{-1}$)')
axs[0].set_ylabel('Intensity (photon/pixel/frame)')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
q_log_ticks(axs[0])   # after set_xscale, which resets the locators
add_minor_grid(axs[0])
# Close the axis up on the data; the curves run diagonally, so the key fits in
# the empty lower-left corner.
I_hi = -np.inf
I_lo = np.inf
for key in ORDER:
    s = data[key]['saxs_1d']
    I_hi = max(I_hi, s.max())
    I_lo = min(I_lo, s[s > 0].min())
axs[0].set_ylim(I_lo * 0.7, I_hi * 1.6)
axs[0].legend(title='Flux on sample', loc='lower left')

axs[1].set_xlabel(r'Delay Time, $\Delta t$ (s)')
axs[1].set_ylabel(r'$g_2$')
axs[1].set_xscale('log')
add_minor_grid(axs[1])
axs[1].legend(title='Flux on sample')

fig.tight_layout(pad=0.4, w_pad=1.4)
save_fig(fig, os.path.join(HERE, 'FigureS6_Flux_Control.pdf'))
