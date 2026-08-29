"""Incident photon flux at 8-ID-I from the calibrated PIN-diode current.

This is the flux-conversion formula behind Figure S5.  Running the file prints
the unattenuated full-beam flux for the single PIN-diode reading recorded on
2025-02-01 (see PINdiode_Calibration_DP00429_2024-02-01.pdf); the per-condition
fluxes in the figure legend are computed by g2_SAXPCS_Rad_Cali.py, which
reimplements the same formula and drives it from each measurement's own
ion-chamber reading instead.

    beam power   = PIN current / spectral responsivity
    raw flux     = beam power / (photon energy in joules)
    corrected    = raw flux / (mirror and window losses)

NOTE: the constants below describe the post-APS-Upgrade 8-ID configuration used
for this control measurement (2025), not the 2022 main experiment.
"""

Pind_Current = 32.7e-6*5.10              # PIN current [A]: reading x gain
Photon_Energy = 12.4e3                      # Photon energy in eV
Spec_Res = 0.1882                        # Si PIN-diode spectral responsivity [A/W]
loss_mirror = 0.94*0.94                  # Loss from mirror reflectivity
loss_diamond = 0.735                     # Loss from BPM and CVD diamond window at 8-ID-E
norm_ring_current = 1.0                  # ring-current normalisation (1.0 = none)
loss_factor = loss_mirror*loss_diamond*norm_ring_current

# Constants for flux calculation
ev2J = 1.6e-19                           # Conversion from eV to Joule  
E_0p1bw = 9.89e-05/1e-3                  # From Xianbo's calculation on 07/22/2025  

# Calculation of raw and corrected flux
Beam_Power = Pind_Current/Spec_Res
Raw_Flux = Beam_Power/(Photon_Energy*ev2J)       # Flux in photons per second
Corrected_Flux = Raw_Flux/loss_factor


print(f"\n The measured raw flux of the 8-ID-I beamline at {Photon_Energy/1e3:.1f} keV is: {Raw_Flux:.2e} \
photons/sec")
