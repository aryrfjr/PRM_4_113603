# PRM_4_113603

Below is an illustration of the MLOps workflow in terms of a Generate+ETL (GETL) framework used in Phys. Rev. Materials 4, 113603 (DOI: https://doi.org/10.1103/PhysRevMaterials.4.113603):

![MLOPs workflow used in PRM_4_113603](../img/PRM_4_113603_MLOps.drawio.png)

Next, a description of the shell scripts that make up the implemented MLOps workflow used in that work. They were used at different stages of the research (tests, development, and application) and are operational at a high technical level manual execution of the Python/Shell-based GETL pipelines depicted in the figure above:

- **zca-bd-full-SD-cpu.sh**: this script creates a LAMMPS input file for a Classical Molecular Dynamics like the [example for the 100-atom cell of Zr₄₉Cu₄₉Al₂](../data_examples/G/big-data-full/Zr49Cu49Al2/c/md/lammps/100/21/Zr49Cu49Al2.lmp.inp). It is a job script for the queue **cpu_long** of the supercomputer SDumont. It is used in the **Generate** step of the MLOps workflow depicted in the figure above.

- **zca-bd-full-SD.sh**: this script performs the same task as **zca-bd-full-SD-cpu.sh**, but in the queue **cpu_long** of the supercomputer SDumont.

- **zca-bd-full.sh**: this script performs the same task as **zca-bd-full-SD-cpu.sh**, but in the supercomputer CENAPAD-SP.
