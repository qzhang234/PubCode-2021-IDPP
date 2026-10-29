"""Figure S2: the ion-chamber calibration behind every absolute intensity here.

Two panels:
  (a) the air-path transmission across the attenuation series, and its mean;
  (b) the PIN-diode photon count against the upstream ion-chamber readout,
      with the linear fit whose slope and intercept are CAL_A and CAL_B.

READS NO DATA FILES.  Everything plotted is a constant in abs_xsec.py:
CAL_UPIC, CAL_DNIC and CAL_PHOTONS are the eight attenuation points recorded
during the calibration, and AIR_TRANS_SERIES, AIR_TRANSMISSION, CAL_A and
CAL_B are derived from them at import.  That matters here: the calibration was
measured in the 2022-1 beamtime, whose averaged files moved to the companion
repository with the rest of the structural work, but the calibration itself
never lived in those files.  The same CAL_A, CAL_B and AIR_TRANSMISSION feed
abs_xsec_coef(), so this figure documents the conversion the isothermal series
actually uses.

Split out of saxpcs.py, which drew it as one of three figures and did read the
2022-1 averages for the other two.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MEW, LW_THIN, LW_DATA,      # noqa: E402
                              apply_style, add_minor_grid, label_panels, save_fig)
from abs_xsec import (AIR_TRANS_SERIES, AIR_TRANSMISSION, AIR_TRANS_STD,  # noqa: E402
                      CAL_UPIC, CAL_PHOTONS, CAL_CROP, CAL_A, CAL_B,
                      calibration_summary)

FIG_SIZE = (DOUBLE_COL, 2.90)

print(calibration_summary())

apply_style()
fig, (axt, axc) = plt.subplots(1, 2, figsize=FIG_SIZE)

# (a) air-path transmission.  The unattenuated point is dropped from the fit
# (CAL_CROP), so it is dropped here too and the series shown is the one the
# mean is taken over.
axt.plot(np.arange(len(AIR_TRANS_SERIES)), AIR_TRANS_SERIES, 'ko:',
         ms=MS, mew=MEW, lw=LW_THIN)
axt.axhline(AIR_TRANSMISSION, color='b', ls='-', lw=LW_DATA,
            label=f'mean = {AIR_TRANSMISSION:.4f} $\\pm$ {AIR_TRANS_STD:.4f}')
axt.set_xlabel('Measurement #')
axt.set_ylabel('Transmission of air')
add_minor_grid(axt)
axt.legend(loc='lower left', borderaxespad=0.4)

# (b) ion chamber to photon flux.
_up = CAL_UPIC[CAL_CROP:]
axc.plot(_up, CAL_PHOTONS[CAL_CROP:], 'ko', ms=MS, mew=MEW, label='measured')
_uline = np.linspace(0, _up.max() * 1.05, 100)
axc.plot(_uline, CAL_A * _uline + CAL_B, 'b-', lw=LW_DATA,
         label=f'fit: {CAL_A:.2e}$\\cdot$Up\\_IC {CAL_B:+.2e}')
axc.set_xlabel('Upstream ion chamber')
axc.set_ylabel('Number of photons')
add_minor_grid(axc)
axc.legend(loc='upper left')

label_panels([axt, axc])
fig.tight_layout()
save_fig(fig, 'FigureS2_Calibration.pdf')
