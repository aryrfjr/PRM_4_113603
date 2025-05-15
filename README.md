# PRM_4_113603

Below is an illustration of the MLOps workflow in terms of the Generate+ETL (GETL) framework used in Phys. Rev. Materials 4, 113603 (DOI: https://doi.org/10.1103/PhysRevMaterials.4.113603):

![MLOPs workflow used in PRM_4_113603](img/PRM_4_113603_MLOps.drawio.png)

The structures of the two directories **ML** and **SS-ML** are the same as those used in the original Python/Shell-based GETL pipelines depicted in the figure above. Below a description of all directories in this repository:

- **ML**: this directory contains all Python scripts that make up the **pre-Deployment** pipeline, which consists of the following steps:
  
  1. **Generate**: this is the key data acquisition step. It is optional in regular MLOps workflows where raw data already exists and is simply collected. In the case of the research work described here, the raw ifnormation of the materials under study (structures and bond data) were generated from scratch using physics-based models. The corresponding simulations are automated and reproducible, where Classical Molecular Dynamics (CMD) simulations are used to generate ensembles of 100-atom cells that are statistically equivalent (in terms of local environments) to the actual structures of the materials studied. The atomistic structure in these 100-atom cells resulting from these CMD simulations are used for the generation of the following synthetic training data:
     
     - the **bond strengths between pairs of atoms** using the Crystal Orbital Hamilton Population (COHP) method via electronic structure calculations within the Density Functional Theory (DFT) framework. An integrated COHP curve yields a value (labeled as -ICOHP) that represents the target variable for the Machine Learning (ML) model (i.e., bond strength).
       
     - the **atomic-centered feature vectors** for the ML model using the Smooth Overlap of Atomic Positions (SOAP) descriptors.
       
  2. **Extract**: Bond pairs, –ICOHP values ... Derive –ICOHP values, compute SOAP descriptors ... Raw data is parsed from files/databases/simulations. ... Feature extraction ... Structured feature engineering from raw inputs, part of a reproducible data process ...
     
  3. **Transform**: SOAP features, normalize, etc. ... Data is cleaned, featurized (e.g., SOAP), scaled, etc. ...
     
  4. **Load**: Build a Data Base of Interactions (DBI) for each interaction type ... Transformed data is saved or passed to the next stage (e.g., model) ...
     
  5. **Train**: GPR model using DBI ... Modeling The actual model is trained on the prepared features ... Train Gaussian Process Regression (GPR) with hyperparameter tuning ... Model development and training ...
      
  6. **Tune**: RMSEs, tune hyperparams, data augmentation ... Cross-validation & RMSE analysis: Evaluated prediction accuracy and determined optimal training sizes (~600 interactions per type). ... Data augmentation: Applied geometric transformations (shear, tension, compression) to generate more structural diversity and improve ML generalization. ... Includes experimentation, validation, and tuning ... 

- **SS-ML**: this directory contains all Python scripts that make up the **Deployment/Production** pipeline for the ML model ...  inference and analysis ... Apply trained ML model to large CMD structures and analyze mechanical behavior and bond exchange ... Model deployment & inference ... Batch inference on large-scale data, output analysis, and visualization ... 

- **data_examples**:

- **shell_scripts**: this directory contains all Shell scripts used to execute the Python scripts in the directories above in both local and remote HPC environments.

- **img**: some images for the README files.
