# Causality-in-QGP
Programs to realize the causality analysis in the quark-gluon plasma (QGP) in the hydrodynamics stage running with MUSIC

# Save the coefficients with MUSIC:
- After the usual installation of the MUSIC, we make some changes in the MUSIC to save the hydrodynamics coefficients necessary for a causality analysis.
- Make this substitutions: 

  grid_info.cpp: in the source folder of MUSIC. It serves for MUSIC to produce two new files called evolution_Wmunu.dat (contain the components of the shear stress tensor, pressure, speed of sound squared, entropy, energy, temperature, and the two viscosities) and evolution_bulk_pressure.dat (contain the bulk);

  evolve.cpp: in the same folder of MUSIC. It defines the flag 4, which is equivalent to the part wroten in the file grid_info.cpp.

- Compile MUSIC.
- Make the substitution:

  music_input_mode_2: In the folder that you will run MUSIC you have a parameters folder. Here, you can control the configurations for using flag 4, you can define the often for the data will be saved along the hydrodynamic evolution and other specific configurations.

- The result will be two files: evolution_Wmunu.dat and evolution_bulk_pressure.dat

# Make the causality analysis:
- causality.py: Read the result files from MUSIC and apply in the causality equations. This program can be changed as the final result wanted 
