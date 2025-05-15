# Examples of structured data in terms of G+ETL

Below is an illustration of the MLOps workflow in terms of a Generate+ETL framework used in Phys. Rev. Materials 4, 113603 (DOI: https://doi.org/10.1103/PhysRevMaterials.4.113603):

![MLOPs workflow used in PRM_4_113603](img/PRM_4_113603_MLOps.drawio.png)

Next, a description of the content of each folder under this **data_examples** directory:

- **G**: this directory contains an example of the set of files generated in the first step (**Generate**) of the Pre-Deployment pipeline for the composition Zr₄₉Cu₄₉Al₂ with:

  - A cubic cell with 100 atoms resulted from a Molecular Dynamics (MD) simulation (**<ID_RUN>** = **21**) which started with a liquid at 2000 K and finished with the solid at 300 K.
    
  - The electronic structure simulation for the last step of that MD simulation (**\<SUB_STEP>** = **0**) and corresponding ICOHP values computed for each interaction.

- **ETL_model**: 

- **TT**: this directory contains the set of files generated for the same example in the step **Generate** for the step **Train/Tune**.

- **ETL_inference**: 
