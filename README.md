# SUTRA-MS-FT
Coupled 2D hydro-thermal-solute modeling in frozen soils

SUTRA-MS-FT is extended SUTRA-MS1.1 considering freeze-thaw processes (FT).

This code repository provides the complete source code for research on soil water-heat-salt transport under freeze-thaw conditions, along with 4 examples dedicated to verifying the core capabilities of four different models respectively.

# Quick Run Guide
1. Locate the 'bin' directory in the repository and copy the executable files (.exe) within it.
2. Paste the copied '.exe' files into the 'root directory of code execution' (specifically, the folder containing the corresponding .inp and .ics files), then double-click the files to launch the program (Note that a Fortran compilation and runtime environment is required).
3. After the program runs successfully, call the corresponding MATLAB script files (`.m` files) for each example to reproduce the visualization results of model verification (Note the file storage location).

# Source Code Modification and Compilation (Fortran Specific)
All source code files of the models are stored in the `source` directory. This code is written in Fortran.
1. Open the source code projects/files with Visual Studio (or other compatible Fortran compilers).
2. Modify the source code as needed, then recompile it to generate new `.exe` executable files after debugging.
