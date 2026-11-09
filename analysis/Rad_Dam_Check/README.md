# Radiation-damage / flux-dependence control

Generates SI Fig. S3, the flux-control figure (file `FigureS6_Flux_Control.pdf`): SAXS and SA-XPCS of a diluted (VPAVG)30 sample
(H06, 20 mg/mL, held at 33.0 °C) measured at four attenuator settings after the
APS Upgrade. If the reported dynamics were a beam-damage artifact, the
correlation functions would change with flux; they do not.

| File | Role |
|---|---|
| `g2_SAXPCS_Rad_Cali.py` | Figure source. Reads the four result files in `cluster_results/`, converts each one's ion-chamber readout to the on-sample photon flux, and plots SAXS `I(Q)` and `g2(τ)` keyed by that flux. Run by `make figures`. |
| `Flux_Cal.py` | The PIN-diode-current → photon-flux formula and its constants. `g2_SAXPCS_Rad_Cali.py` reimplements the same formula; this file documents it and prints the full-beam flux for the recorded calibration reading. |
| `Pind_calibration.csv` | Calibration sheet: upstream ion chamber against calibrated PIN-diode current. Fit through the origin gives the IC → PIN slope. |
| `PINdiode_Calibration_DP00429_2024-02-01.pdf` | Vendor/beamline calibration certificate for the PIN diode. |
| `cluster_results/` | The four result files, E0171–E0174, committed so the figure rebuilds from this repository alone. Byte-identical to the copy of record in `/gdata/s8id-dmdtn/2025-1/qzhang202503/flux_control_results`. |

The nominal attenuator labels (a7, a9, a11, a14) are not monotonic in flux
because the attenuator setup was wrong, so the true on-sample flux is recovered
per measurement:

```
IC readout (entry/instrument/incident_beam/incident_beam_intensity, x 1e-9)
    --(Pind_calibration.csv, fit through origin)--> calibrated PIN current [A]
    --(Flux_Cal.py formula)-----------------------> photon flux on sample [ph/s]
```

`g2_SAXPCS_Rad_Cali.py` also prints every number the SI caption and the text
quote: the monitor and mean-intensity ratios, the fluxes relative to the main
experiment, and the Q values of both panels.

## Provenance

The docstring of `g2_SAXPCS_Rad_Cali.py` records the same. The four result
files were correlated in 2025 from raw frames under
`/gdata/dm/8ID/8IDI/2025-1/qzhang202503/data` with the qmap
`rigaku_s360_d36_lin_MC.hdf`; each file records both paths. That DM experiment
was archived on 23 Aug 2025, so the raw frames are no longer online, and
restoring them needs the 8-ID DM administrator. The copy of record of the
results has been `/gdata/s8id-dmdtn/2025-1/qzhang202503/flux_control_results`
since 1 Nov 2026. The same-named files in `.../qzhang202503/Multitau` are an
older reduction with a different qmap (`rigaku3m_qmap_S360_D36_Lin_refined.h5`)
and are not used.
