# PRM_4_113603

Below is an illustration of the MLOps workflow in terms of the Generate+ETL (GETL) framework used in **Phys. Rev. Materials 4, 113603** (DOI: https://doi.org/10.1103/PhysRevMaterials.4.113603; or the [preprint](https://www.researchgate.net/publication/345634787_Chemical_bonding_in_metallic_glasses_from_machine_learning_and_crystal_orbital_Hamilton_population)):

![MLOPs workflow used in PRM_4_113603](img/PRM_4_113603_MLOps.drawio.png)

Below a description of the directories in this root folder:

- **scripts**: this folder contains the essential Python and Bash scripts which make up the implemented MLOps workflow used in that work. The scripts were used at different stages of the research (tests, development, and application) and are operational at a high technical level manual execution of the Python/Shell-based GETL pipelines depicted in the figure above.

- **data_examples**: this folder contains samples of the data generated for the nominal composition Zr₄₉Cu₄₉Al₂.

- **atomistic_models**: this folder contains simulation-specific model inputs used for the generation of **automated/reproducible raw data** from atomistic simulations using LAMMPS (a classical molecular dynamics simulator) and Quantum ESPRESSO (a first-principles electronic structure code based on density functional theory, DFT). These models inputs are EAM potentials and pseudopotentials.

- **img**: this folder contains some figures used in the README files of this repository.

**NOTE**: Under each one of the directories above, the original folder structure used to produce the results reported in **Phys. Rev. Materials 4, 113603** are replicated. The two root folders **ML** and **SS-ML** (and its subfolders) are the same as those used in the original Python/Shell-based GETL pipelines depicted in the figure above.
