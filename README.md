# HASTools.jl

This is a set of tools to use while running helium atom scattering (HAS) experiments and for subsequent data analysis.

The intended workflow looks like:

1) Make an `Experiment` to handle He beam information and all collected data
2) Measure diffraction peaks and store the data
3) Fit the diffraction peaks with a Gaussian, Lorentzian, or pseudo-Voigt model
4) Use fit to determine surface orientation / structural information
5) Look at scan curves to decide on angles to measure TOF data
6) Measure TOF data and store
7) Fit the inelastic events in TOF data with Gaussian, Lorentzian, or pseudo-Voigt models
8) Use fits to determine surface phonon (**q**, ω<sub>**q**</sub>)
