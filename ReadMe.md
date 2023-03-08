# Reference Dislocation-density based Crystal Plasticity Finite Element(CPFE) Model
**Copyrights: Damin Xia and Caglar Oskay**

This folder contains the codes and input files for reference DD-CPFE simulation. The code is developed based on DDEHM code, the idea is to apply the consititutive relation from DDEHM for every qudrature point. See formulation.pdf for theory details.

The code is developed under small deformation assumption. Please be careful when using it under large deformation (i.e., strain > 5%).

**What are the necessary input files?**
1. ABAQUS inp file (e.g. Direct_9grain.inp). Note that for a quasi-2D microstructure, the DOF of the out-of-plane direction is constrained.
2. elasmod.dat: elastic coefficients for the material, which could contain the elastic coefficients for multiple phases.
3. phase1.dat and phase2.dat: record which grain belongs to which phase.
4. NBgrain.dat: record neighboring grain information. In the context of direct simulation, it will be used in the postprocessing stage.
5. test.xtali: contains material paramters and model solver data.
6. Texture.txti: contains grain orientation information.
7. (optional) ABAQUS cae file (e.g. 9RVE_HEX.cae).

**Which files are the codes?**
1. UMAT: IntegrateCrystalEqns.f90, ROMUtility.f90, SetUpCrystalProps.f90, StressSolve.f90, uexternaldb.f90, Utility.f90, Modules.f90 and umat.f.
2. Postprocessing code: eng_stress_strain.m, ReadData.py.

**How to excute the code and obtain the results**
run the bash file (./run_tension.sh), then excute ReadData.py to extract the stress-strain data. Run eng_stress_strain.m to plot the data.

**Future development**
1. Code strcture can be improved, it is not a very efficient code right now.
2. Bug may still exist, the error comparing with DDEHM model is over 5%.

*Make sure to commit all changes*





