# Causality-in-QGP
Programs to realize the causality analysis in the quark-gluon plasma (QGP) in the hydrodynamics stage running with MUSIC

# Save the coefficients with MUSIC:
- After the usual installation of the MUSIC, we make some changes in the MUSIC to save the hydrodynamics coefficients necessary for a causality analysis.
- Make these substitutions: 

  grid_info.cpp: in the source folder of MUSIC. It serves for MUSIC to produce two new files called evolution_Wmunu.dat (containing the components of the shear stress tensor, pressure, speed of sound squared, entropy, energy, temperature, and the two viscosities) and evolution_bulk_pressure.dat (containing the bulk);

  evolve.cpp: in the same folder of MUSIC. It defines the flag 4, which is equivalent to the part wroten in the file grid_info.cpp.

- Compile MUSIC.
- Make the substitution:

  music_input_mode_2: In the folder that you will run MUSIC, you have a parameters folder. Here, you can control the configurations using flag 4 and you can define the often the data will be saved along the hydrodynamic evolution or other specific configurations.

- The result will be two files: evolution_Wmunu.dat and evolution_bulk_pressure.dat

# Make the causality analysis:
- causality.py: Read the result files from MUSIC and apply them to the causality equations. This program can be changed as the final result wanted 

# Example of result:
<div align="center">
<img src="https://user-images.githubusercontent.com/117451854/200028241-6284e2a4-8902-4790-b839-511ae7042d88.png" width="700px" />
</div>
Plots for a collision Pb-Pb at 5.02 TeV, different instants of the hydrodynamics evolution are showed. Points in red are acausal, in purple are indeterminate and in blue are causal.

# Theoric references:
- Bayesian analysis used: 

  MORELAND, J. Scott; BERNHARD, Jonah E.; BASS, Steffen A. Bayesian calibration of a hybrid nuclear collision model using p− Pb and Pb-Pb data at energies available at the CERN Large Hadron Collider. Physical Review C, v. 101, n. 2, p. 024911, 2020.
- Calculation of coefficients:
  
  https://webhome.phy.duke.edu/~jp401/music_manual/hydro.html
- Calculation of causality violation equations:
  
  BEMFICA, Fábio S. et al. Nonlinear constraints on relativistic fluids far from equilibrium. Physical review letters, v. 126, n. 22, p. 222301, 2021.
