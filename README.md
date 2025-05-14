# PRM_4_113603

Below is an illustration of the MLOps workflow in terms of a Generate+ETL pipeline used in Phys. Rev. Materials 4, 113603 (DOI: https://doi.org/10.1103/PhysRevMaterials.4.113603):

![MLOPs workflow used in PRM_4_113603](img/PRM_4_113603_MLOps.drawio.png)

Next, a description of the scripts used in that work for a high technical level operation with a manual Python/Shell-based GETL+MLOps pipeline:

The directory **ML** contains the database of structures generated to train the ML model. It contains the following subdirectories with scripts:

- **misc**:

  - **1st_test_quippy/1st_test_quippy.py**: first tests with SOAP descriptors discussed in diaries of **Aug-2018** and **Sep-2018**.

  - **2nd_test_quippy/2nd_test_quippy.py**: more tests with SOAP descriptors discussed in diary of **Sep-2018**.

  - **3rd_test_quippy/3rd_test_quippy.py**: more tests with SOAP descriptors discussed in diaries of **Sep-2018** and **Oct-2018**.

  - **4th_test_quippy/4th_test_quippy.py**: more tests with SOAP descriptors discussed in diary of **Oct-2018**.

  - **4th_test_quippy/setup_jobs.py**: more tests with SOAP descriptors discussed in diary of **Oct-2018**.

- **big-data-full**:

  - **scripts/setup_ICOHP_jobs.py**: This script will setup the Slurm job files to run in parallel 15 ICOHP calculations, 15 for each **<ID_RUN>**. Recalling that only the last **<STEP>** will be considered, from which 14 **<SUB_STEP>** will  be derived to generate the corresponding distorted structures (see on **05/10/2018-(10)**). For instance, the directory **./big-data-full/Zr49Cu49Al2/c/md/lammps/100/42** contains a MD simulation for a **random** distribution of 100 atoms (**<ID_RUN>** is 42). In that MD simulation that initial structure was thermalized (thermal equilibrium) at 2000 K for a while (see more on **03/11/2016-(10)**) and then quenched to 300 K and thermalized. That final structure (the last MD **\<STEP>** labeled as 2000) was taken to generate the 15 structures (**<SUB_STEP>** from 0 to 14) for which the simulations with Lobster will run. **NOTE**: the work done by this script can be done by some **GenAI**.

  - **scripts/check_QE_ICOHP.py**: The electronic structure calculation and post-processing for computing the -ICOHP values may not converge or get finished successfully. This script checks which steps and sub-steps of a given run have been finished according to the following directory structure (see on **05/10/2018-(10)**):
    - **<CHEM_COMPOSITION>/c/md/lammps/100/<ID_RUN>/\<STEP>/<SUB_STEP>**.

  - **scripts/compare_clusters_soap.py**: This is script was never mentioned in any diary. It basically compares two clusters (central atoms plus first coordination shells) in their individual **.xyz** files.

  - **scripts/compare_global_soaps.py**: This is script is mentioned in diary of **05/11/2018**. It basically compares two structures using **QUIP** directly instead of **theo4m**.

  - **scripts/compare.py**: This is script was never mentioned in any diary. Basically in its current state, it compares individual SOAPs like a file **SOAPS-Al.vec**.

  - **scripts/create_SSDB.py**: This is script was never mentioned in any diary. It creates a single-SOAP (SS) database for all **<ID_RUN>s** in a given nominal composition (see on **05/10/2018-(10)**).

  - **scripts/create_PBMSDB.py**: This script creates a per-bond (PB) multiple-SOAP (MS) database for all **<ID_RUN>s** in a given nominal composition. The code follows the directory structure below (see on **05/10/2018-(10)**), with **<SUB_STEP>s** ranging from 0 to 14:
    - **<CHEM_COMPOSITION>/c/md/lammps/100/<ID_RUN>/2000/<SUB_STEP>**.
    - NOTE: As per **07/07/2019-(3)**, the structure of the files in the per-bond single-SOAP database (PBSSDB) and in the new per-bond multiple-SOAP database (PBMSDB) will be the same. The only difference is that the former will use the file **SOAPS.vec** of each **<SUB_STEP>** whereas the latter will use the individual ones **SOAPS-Zr.vec**, **SOAPS-Cu.vec**, and **SOAPS-Al.vec**.

  - **scripts/avgsoaps.py**: Just some tests with the class **SOAPList** defined in the file **atom/desc/soap/soaplist.py** (see on **19/12/2018-(2)**).

  - **scripts/check_ICOHPLIST.py**: This script checks for a given run, step, and sub-step if all the interactions in the reference **lobsterin-quippy** file exist in the respective **ICOHPLIST.lobster** (see on **17/06/2019-(5)**). It uses the corresponding **.xyz** file to write individual **.xyz** files for each cluster (central atom plus first coordination shell).

  - **scripts/mix_SSDBs.py**: This is a script that creates mixed databases (see on **03/08/2019-(1)**) reading from the per-bond single-SOAP databases (PBSSDBs; see on **07/07/2019-(3)**). An example of outputs generated can be found in the subdirectory **DB1**.

  - **scripts/PBSSDB-kernel_fit.py**: This script fits a per-bond (PB) single-SOAP (SS) kernel for -ICOHP values. It reads the files **.bnd** and **.vec** from the per-bond single-SOAP databases (PBSSDBs; see on **07/07/2019-(3)**) and writes the target vs. predicted -ICOHP values (with RMSE information), like those depicted in **Fig. 2** (see on **02/08/2019-(4)**). The set of RMSE values for different runs of this script with same ML model parameters were used to plot the convergence tests in **Fig. 1**.

  - **scripts/PCSSDB-kernel_fit.py**: This is script was never mentioned in any diary. Basically it fits a per-cluster (central atoms plus first coordination shells) GPR kernel for sum(-ICOHP). Based on diary of **03/07/2019**, it was before Prof. Gábor suggested the **bond feature vector**. **NOTE**: the per-cluster databases are in the **-PCDB** directories whilst the per-bond ones will be saved in the **-PBDB** directories.

The directory **SS-ML** contains the simulations of the stress-strain curves depicted in **Fig. 4**, in which I applied the GPR ML model to 10k atoms cells. It contains a subdirectory **scripts** with:

- **fcsci.py**: This is a script that reads a step from a LAMMPS .dump file to compute their FCSI indexes (see on **19/09/2019**).

- **ML-ICOHP_MPI.py**: This is a script that reads a set of steps from a LAMMPS .dump file to compute the ML predicted ICOHP values for all detected bonds. This is script was never mentioned in any diary.

- **ML-ICOHP.py**: This is a script that reads a set of steps from a LAMMPS .dump file to compute the ML predicted ICOHP values for all detected bonds. Mentioned so many times in the diary of **Aug-2019** and also on **27/09/2019-(4)** and on **06/02/2020-(3)**.

- **ML-ICOHP_RMSD.py**: This is a script that reads a LAMMPS .dump file and writes a corresponding file with ML predicted ICOHP values for all detected bonds. The user can define which steps will be read as well as the cutoff radius which in principle has to be the same used in the DFT ICOHP simulations. This code also reports the RMSD values for the last step based on its ICOHPLIST.lobster file. See on **08/08/2019-(3)**

- **nucleation.py**: See on **09/09/2019-(1)**.

- **part-ICOHP.py**: See on **12/08/2019-(9)**. Mentioned so many times in the diary of **Aug-2019**.

- **plot_2D-bin.py**: This is script was never mentioned in any diary.

- **prepare_R_data_frame.py**: This is a script that reads a \<STEP>_<TRUE_STEP>_ICOHP.lobster file generated with the script ML-ICOHP.py and creates a per-interaction "data frame" as described on **07/08/2019-(5)**. See also **27/09/2019-(6)**.

- **3DI_ICOHP.py**: This is a script that simply interpolates the ICOHP values localized between two atoms using true values computed with LOBSTER for the 100-atoms cells. Mentioned so many times in the diary of **Feb-2020**.
