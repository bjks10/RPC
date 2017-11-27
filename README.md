# RPC
The Robust Profile Clustering (RPC) model is a population-based clustering technique that adjusts for differences that may exist within different subpopulations. Here, participants are clustered at two levels: (1) globally, where subjects are assignd to an overall population-level cluster using an overfitted mixture model, and (2) locally, where variations to global patterns are accommodated via a Beta-Bernoulli process dependent on subpopulation differences.

## Getting Started
The code and supporting materials are run using MATLAB software. To run the example data, you will need:
* RPCexample.m - run RPC and create ouput data
* drchrnd.m - Dirichlet random generator function
* heatmap.m - function file to generate desired heatmap figures
* ExampleData.mat - input data source 

## Example Data
The example dataset is a MAT-file that contains the following variables:
* subdata: 1800x50 matrix. This matrix is the input dataset containing subject level data for 50 variables. Each variable assigned a single categorical value (1,2,3,4). 
*	sub_nu: binary 50x3 matrix. This matrix is used as a reference to illustrate the true probability of allocation for each variable within each subpopulation to global (ν= 1) or local (ν= 0). 
* subpop_i: 1800x1 vector. This vector contains subpopulation ID for the 1800 subjects included in the dataset.
* Subpop1_true: 50x3 matrix. This matrix contains the 3 modal cluster patterns modally expected from subpopulation 1. Each column represents a different global pattern with a deviation reflected in 13 of the 50 variables. 
* Subpop2_true: 50x3 matrix. This matrix contains the 3 modal cluster patterns modally expected from subpopulation 2. Each column represents a different global pattern with a deviation reflected in 24 of the 50 variables.
* Subpop3_true: 50x3 matrix. This matrix contains the 3 modal cluster patterns modally expected from subpopulation 3. Each column represents a different global pattern with no deviation in any of the 50 variables. 

## Example Output
Execution of the RPCexample.m file with Example dataset should generate the following results:
* Posterior median estimates of all RPC model parameters (π,λ,θ_0,θ_1,ν,β) saved as RPCparm_ex.mat
* Text file outputs of π,ν,β (pis.txt, nus.txt, betas.txt).
* Text file output of predicted modal pattern of global clusters (GlobalPattern_predicted.txt)
* Trace plot figure of π (pis.png)
* Trace plot figure of β (betas.png)
* Dendrogram plot from Posterior pairwise label switching step (dendrogram.png)
* Heat map comparison side-by-side plots of predicted deviations by subpopulation and actual deviations by subpopulation. (G_deviations.png)

## Authors
**Briana Stephenson**, Amy Herring, Andrew Olshan


