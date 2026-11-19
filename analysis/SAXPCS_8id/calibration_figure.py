"""SI Fig. S2 (file FigureS2_Calibration.pdf): the ion-chamber calibration
behind every absolute intensity.

(a) the transmission of the air path across the attenuation series, and its
mean; (b) the PIN-diode photon count against the upstream ion-chamber reading,
with the linear fit whose slope and intercept are CAL_A and CAL_B.

Reads no data file: everything plotted is a constant of abs_xsec.py, the same
numbers abs_xsec_coef() uses.  The calibration was recorded in a separate
beamtime on the same ion chambers.
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

# (a) air-path transmission, over the seven points the mean is taken over
axt.plot(np.arange(len(AIR_TRANS_SERIES)), AIR_TRANS_SERIES, 'ko:',
         ms=MS, mew=MEW, lw=LW_THIN)
axt.axhline(AIR_TRANSMISSION, color='b', ls='-', lw=LW_DATA,
            label=f'mean = {AIR_TRANSMISSION:.4f} $\\pm$ {AIR_TRANS_STD:.4f}')
axt.set_xlabel('Measurement #')
axt.set_ylabel('Transmission of air')
add_minor_grid(axt)
axt.legend(loc='lower left', borderaxespad=0.4)

# (b) ion chamber to photon count, over the seven attenuated points
up_ic = CAL_UPIC[CAL_CROP:]
axc.plot(up_ic, CAL_PHOTONS[CAL_CROP:], 'ko', ms=MS, mew=MEW, label='measured')
up_line = np.linspace(0, up_ic.max() * 1.05, 100)
axc.plot(up_line, CAL_A * up_line + CAL_B, 'b-', lw=LW_DATA,
         label=f'fit: {CAL_A:.2e}$\\cdot$Up\\_IC {CAL_B:+.2e}')
axc.set_xlabel('Upstream ion chamber')
axc.set_ylabel('Number of photons')
add_minor_grid(axc)
axc.legend(loc='upper left')

label_panels([axt, axc])
fig.tight_layout()
save_fig(fig, os.path.join(HERE, 'FigureS2_Calibration.pdf'))
