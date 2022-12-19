# SPADAN
SPADAN stands for Systematic Protein Association Dynamic ANalyzer which is an old name of a beautiful city in Iran!
SPADAN is a systematic pipeline developed in MATLAB that generates a large-scale ordinary differential equation model from protein interactions in a biological phenomenon by getting time-series experimental concentrations from it. 


1. [Prerequisites](#1-prerequisites)
2. [Description](#2-description)
3. [Usage](#3-usage)
4. [Run the pipeline](#4-run)

------
## 1 Prerequisites

In order to run the pipeline, it is needed to install MATLAB 2015b.
| NOTE                                                         |
| :----------------------------------------------------------- |
| Althugh the SPADAN is optimized to use minimum value of memory, in order to use the pipeline for large-scale networks, the computational hardware should have enough RAM capacity. |


------
## 2 Description
The details about databases used in the pipeline, imortant functions, and important variables are described below.

### 2.1 Databases

| Database name| Pipeline step | Description |
| ------ | ----- | ----- |
| Biomart | STEP I | The database is used to link Ensemble protein IDs, Ensemble gene IDs, Ensemble transcription IDs, and gene names |
| SIGNOR2.0 | STEP I   | The database is used to extract experimentally validated interactions between input genes |

### 2.2 Important functions

| Function name name| Pipeline step | Description |
| ------ | ----- | ----- |
| phos_afct, dephos_afct, ubiq_afct, cmplx_afct, exp_afct, deg_afct, prod_afct | STEP II | Expands the interaction edge to biochemical reactions and adds the recations to the biochemical reactions network.|
| mass2, deg, product | STEP III | automatically writes the list of ODE equations describing the kinetics of each reaction based on mass-action law.|
|optimizer_Ifmin5, optimizer_simplex5_2, optimizer_ILMA5| STEP IV | The modified optimizer code to perform large-scale parameter approximation.|
|dynamics_auto.m*| STEPIV|Systematically generated ODE model that includes the list of model ODEs and is able to be simulated|
|dynamics_auto_v20.m*| STEP IV |Optimized systematically generated ODE model that is able to be simulated|

### 2.3 Important variables

| Variable name| description |
| ------ | ----- |
| signor_edges | The list of all edges existing in SIGNORE databas|
| prtns_nodes_adress | The ID number assigned by SPADAN to each Protein, gene, RNA, PTmodified protein, and protein complexes.  |
| reactions | The 3D matrix that includes all biochemical reactions by the ID of reactants and products.|
| reaction_type | A complementary variable for "reactions" including the information about the type of each interaction.|
|nodes | The list of node names of biochemical reactions network|
| prtns_activity | A matrix that shows whether different PTmod proteins are active or not|
| X_rel_param  | A boolean matrix that shows which parameters are related to each state variable.|
|X| State variables of the model|
|X_dot| Derivative of state variables of the model |
|K* | Values of estimated parameters |
|A| The boolean matrix which is used to calculate ODEs|
|V| The string matrix inculding automatically generated ODEs for each state variable.|
| edges | The list of PPI edges considered in modelling process |
| C | The boolean matrix which is used for calculating model outputs|






------
## 3 Usage
SPADAN inputs must be excel files including experimentally measured concentrations/intensities of biomolecules and their Ensemble IDs. The input files should be placed in the "SPADAN/data files/input data" directory. After running SPADAN, it will automatically put outpus in the "SPADAN/data files/output data" directory.
| NOTE                                                         |
| :----------------------------------------------------------- |
| Some SPADAN outputs such as estimated parameter values, list of ODEs, and simulatable model are available in MATLAB environment after finishing the pipeline. |
------
## 4 Run the pipeline
After navigating to the repository's directory in MATLAB, write the below command in the command window.
```{r}
SPADAN_v21
```
Thus, the pipeline will be runned and a short description about the performing step will be shown in the command window.
| NOTE                                                         |
| :----------------------------------------------------------- |
| Performing pipeline procedure needs no human interactions. |
------
