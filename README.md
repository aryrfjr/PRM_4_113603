# PRM_4_113603
Scripts used in Phys. Rev. Materials 4, 113603 (DOI: https://doi.org/10.1103/PhysRevMaterials.4.113603)

The directory **ML** contains the database of structures generated to train the ML model. It contains the following subdirectories with scripts:

- **misc**:

  - **1st_test_quippy/1st_test_quippy.py**: first tests with SOAP descriptors discussed in diaries of **Aug-2018** and **Sep-2018**.

  - **2nd_test_quippy/2nd_test_quippy.py**: more tests with SOAP descriptors discussed in diary of **Sep-2018**.

  - **3rd_test_quippy/3rd_test_quippy.py**: more tests with SOAP descriptors discussed in diaries of **Sep-2018** and **Oct-2018**.

  - **4th_test_quippy/4th_test_quippy.py**: more tests with SOAP descriptors discussed in diary of **Oct-2018**.

  - **4th_test_quippy/setup_jobs.py**: more tests with SOAP descriptors discussed in diary of **Oct-2018**.

- **big-data-full**:

  - **scripts/mix_SSDBs.py**: This is a script that creates mixed databases (see on **03/08/2019-(1)**) reading from the per-bond single-SOAP databases (PBSSDBs; see on **07/07/2019-(3)**). An example of outputs generated can be found in the subdirectory **DB1**.

  - **scripts/PBSSDB-kernel_fit.py**: This script fits a per-bond (PB) single-SOAP (SS) kernel for -ICOHP values. It reads the files **.bnd** and **.vec** from the per-bond single-SOAP databases (PBSSDBs; see on **07/07/2019-(3)**) and writes the target vs. predicted -ICOHP values (with RMSE information), like those depicted in **Fig. 2** (see on **02/08/2019-(4)**). The set of RMSE values for different runs of this script with same ML model parameters were used to plot the convergence tests in **Fig. 1**.

  - **scripts/check_QE_ICOHP.py**:

  - **scripts/avgsoaps.py**:

  - **scripts/setup_ICOHP_jobs.py**:

  - **scripts/compare_clusters_soap.py**:

  - **scripts/PCSSDB-kernel_fit.py**:

  - **scripts/create_SSDB.py**:

  - **scripts/create_PBMSDB.py**:

  - **scripts/check_ICOHPLIST.py**:

  - **scripts/compare_global_soaps.py**:

  - **scripts/compare.py**:

The directory **SS-ML** contains the simulations of the stress-strain curves depicted in **Fig. 4**, in which I applied the GPR ML model to 10k atoms cells. It contains a subdirectory **scripts** with:

- **3DI_ICOHP.py**: This is a script that simply interpolates the ICOHP values localized between two atoms using true values computed with LOBSTER for the 100-atoms cells. Mentioned so many times in the diary of **Feb-2020**.

- **fcsci.py**: This is a script that reads a step from a LAMMPS .dump file to compute their FCSI indexes (see on **19/09/2019**).

- **ML-ICOHP_MPI.py**: This is a script that reads a set of steps from a LAMMPS .dump file to compute the ML predicted ICOHP values for all detected bonds. This is script was never mentioned in any diary.

- **ML-ICOHP.py**: This is a script that reads a set of steps from a LAMMPS .dump file to compute the ML predicted ICOHP values for all detected bonds. Mentioned so many times in the diary of **Aug-2019** and also on **27/09/2019-(4)** and on **06/02/2020-(3)**.

- **ML-ICOHP_RMSD.py**: This is a script that reads a LAMMPS .dump file and writes a corresponding file with ML predicted ICOHP values for all detected bonds. The user can define which steps will be read as well as the cutoff radius which in principle has to be the same used in the DFT ICOHP simulations. This code also reports the RMSD values for the last step based on its ICOHPLIST.lobster file. See on **08/08/2019-(3)**

- **nucleation.py**: See on **09/09/2019-(1)**.

- **part-ICOHP.py**: See on **12/08/2019-(9)**. Mentioned so many times in the diary of **Aug-2019**.

- **plot_2D-bin.py**: This is script was never mentioned in any diary.

- **prepare_R_data_frame.py**: This is a script that reads a <STEP>_<TRUE_STEP>_ICOHP.lobster file generated with the script ML-ICOHP.py and creates a per-interaction "data frame" as described on **07/08/2019-(5)**. See also **27/09/2019-(6)**.
