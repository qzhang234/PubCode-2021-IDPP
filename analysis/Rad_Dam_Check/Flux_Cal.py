"""Photon flux at 8-ID-I from the calibrated PIN-diode current.

The flux formula behind SI Fig. S3.  Running the file prints the raw flux for
the single PIN-diode reading recorded on 2025-02-01 (see
PINdiode_Calibration_DP00429_2024-02-01.pdf).  g2_SAXPCS_Rad_Cali.py uses the
same formula and constants, driven by each measurement's own ion-chamber
reading.

    beam power = PIN current / spectral responsivity
    raw flux   = beam power / (photon energy in joules)
    corrected  = raw flux / (mirror and window losses)

The constants describe the post-APS-Upgrade 8-ID configuration of this control
measurement (2025), not the 2021 main experiment.
"""

Pind_Current = 32.7e-6*5.10              # PIN current [A]: reading x gain
Photon_Energy = 12.4e3                   # photon energy [eV]
Spec_Res = 0.1882                        # Si PIN-diode spectral responsivity [A/W]
loss_mirror = 0.94*0.94                  # mirror reflectivity
loss_diamond = 0.735                     # BPM and CVD diamond window at 8-ID-E
norm_ring_current = 1.0                  # ring-current normalisation (1.0 = none)
loss_factor = loss_mirror*loss_diamond*norm_ring_current

# The electron charge converts eV to J.  1.6e-19 is e rounded to two digits,
# 0.14% below the exact value, far inside the uncertainty of the responsivity
# and the loss factors.
ev2J = 1.6e-19                           # [C]

Beam_Power = Pind_Current/Spec_Res
Raw_Flux = Beam_Power/(Photon_Energy*ev2J)       # photons per second
Corrected_Flux = Raw_Flux/loss_factor


print(f"\n The measured raw flux of the 8-ID-I beamline at {Photon_Energy/1e3:.1f} keV is: {Raw_Flux:.2e} \
photons/sec")
