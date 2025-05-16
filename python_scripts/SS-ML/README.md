# PRM_4_113603

Below is an illustration of the MLOps workflow in terms of a Generate+ETL (GETL) framework used in Phys. Rev. Materials 4, 113603 (DOI: https://doi.org/10.1103/PhysRevMaterials.4.113603):

![MLOPs workflow used in PRM_4_113603](../img/PRM_4_113603_MLOps.drawio.png)

Next, a description of the scripts that make up the implemented MLOps workflow used in that work. I'm working to make this README file self-contained. The scripts were used at different stages of the research (tests, development, and application) and are operational at a high technical level manual execution of the Python/Shell-based GETL pipelines depicted in the figure above:

The directory **SS-ML** contains the simulations of the stress-strain curves depicted in **Fig. 4**, in which I applied the GPR ML model to 10k atoms cells. It contains a subdirectory **scripts** with:

- **fcsci.py**: this is a script that reads a step from a LAMMPS **.dump** file to compute their FCSI indexes (see on **19/09/2019**).

- **ML-ICOHP_MPI.py**: this is a script that reads a set of steps from a LAMMPS **.dump** file to compute the ML predicted ICOHP values for all detected bonds. This is script was never mentioned in any diary.

- **ML-ICOHP.py**: this is a script that reads a set of steps from a LAMMPS **.dump** file to compute the ML predicted ICOHP values for all detected bonds. Mentioned so many times in the diary of **Aug-2019** and also on **27/09/2019-(4)** and on **06/02/2020-(3)**.

- **ML-ICOHP_RMSD.py**: this is a script that reads a LAMMPS **.dump** file and writes a corresponding file with ML predicted -ICOHP values for all detected bonds. The user can define which steps will be read as well as the cutoff radius which in principle has to be the same used in the DFT ICOHP simulations. This code also reports the RMSD values for the last step based on its ICOHPLIST.lobster file. See on **08/08/2019-(3)**

- **nucleation.py**: see on **09/09/2019-(1)**.

- **part-ICOHP.py**: see on **12/08/2019-(9)**. Mentioned so many times in the diary of **Aug-2019**.

- **plot_2D-bin.py**: this is script was never mentioned in any diary.

- **prepare_R_data_frame.py**: this is a script that reads a **\<STEP>_<TRUE_STEP>_ICOHP.lobster** file generated with the script **ML-ICOHP.py** and creates a per-interaction "data frame" as described on **07/08/2019-(5)**. See also **27/09/2019-(6)**.

- **3DI_ICOHP.py**: this is a script that simply interpolates the ICOHP values localized between two atoms using true values computed with LOBSTER for the 100-atoms cells. Mentioned so many times in the diary of **Feb-2020**.
