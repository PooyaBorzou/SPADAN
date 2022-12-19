# SPADAN
SPADAN stands for Systematic Protein Association Dynamic ANalyzer which is an old name of a beautiful city in Iran!
SPADAN is a systematic pipeline developed in MATLAB that generates a large-scale ordinary differential equation model from protein interactions in a biological phenomenon by getting time-series experimental concentrations from it. 


1. [Prerequisites](#1-prerequisites)
2. [Description](#2-description)
3. [Usage](#3-usage)
4. [Run the pipeline](#4-run)
5. [Publication](#5-publication)

------
## 1 prerequisites

In order to run the pipeline, it is needed to install MATLAB 2015b.
| NOTE                                                         |
| :----------------------------------------------------------- |
| Althugh the SPADAN is optimized to use minimum value of memory, in order to use the pipeline for large-scale networks, the computational hardware should have enough ram memory. |


------
## 2 description
The details about databases used in the pipeline, imortant functions, and important variablres are described below.

### 2.1 Databases

| database name| the pipeline step used it | description |
| ------ | ----- | ----- |
| Biomart | STEP I | The database is used to link Ensemble protein IDs, Ensemble gene IDs, Ensemble transcription IDs, and gene names |
| SIGNOR2.0 | STEP I   | The database is used to extract experimentally validated interactions between input genes |

### 2.2 Important functions

| function name name| the pipeline step used it | description |
| ------ | ----- | ----- |
| phos_afct, dephos_afct, ubiq_afct, cmplx_afct, exp_afct, deg_afct, prod_afct | STEP II | Expands the interaction edge to biochemical reactions and adds the recations to the biochemical reactions network.|
| mass2, deg, product | STEP III | automatically writes the list of ODE equations describing the kinetics of each reaction based on mass-action law.|
|optimizer_Ifmin5, optimizer_simplex5_2, optimizer_ILMA5| STEP IV | The modified optimizer code to perform large-scale parameter approximation.|

### 2.3 Important variables

| function name name| the pipeline step used it | description |
| ------ | ----- | ----- |
| phos_afct, dephos_afct, ubiq_afct, cmplx_afct, exp_afct, deg_afct, prod_afct | STEP II | Expands the interaction edge to biochemical reactions and adds the recations to the biochemical reactions network.|
| mass2, deg, product | STEP III | automatically writes the list of ODE equations describing the kinetics of each reaction based on mass-action law.|
|optimizer_Ifmin5, optimizer_simplex5_2, optimizer_ILMA5| STEP IV | The modified optimizer code to perform large-scale parameter approximation.|


------
## 3 usage

## 4 run
```{r}
optimtools
```
## 5 publication
