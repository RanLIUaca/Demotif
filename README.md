# Demotif
Demotif is an approach to decompose hierarchical motifs.


## Installation
You can install the released version of demotif from GitHub with:
```
devtools::install_github("RanLIUaca/demotif")
library(demotif)
```

## Prior settings
We set the non-informative prior distributions as the default settings for users: $\boldsymbol{\alpha_0}$, $\boldsymbol{\alpha_j}$ and $\widetilde{\boldsymbol{\alpha_j}}$ are three vector whose entries are all one. $\boldsymbol{\pi_{0,a_i}}$ and $\boldsymbol{\pi_{0,b_i}}$ are set to be equal. $p_0$ is set to be 0.5.

## Input
For  `Data_process` function, the input is
* Data: A data frame containing the following columns:
         - Column 1: Sequence data (R), in string format.
         - Column 2: Observed W data (W_obs), numeric, may contain NA values.         
         - Column 3: Observed G data (G_obs), numeric, may contain NA values.

For `de_motif` function, the inputs are
* Data: A data frame containing the sequence, the label vector W, and the label vector G.
* motif_len_w: An integer specifying the length of the first binding motif.
* motif_len_g: An integer specifying the length of the second binding motif.
* mh_step_w: An integer specifying the number of steps to skip between each execution of the MH algorithm for the first binding motif.
* mh_step_g: An integer specifying the number of steps to skip between each execution of the MH algorithm for the second binding motif.
* stop_jump_step: An integer specifying the time step at which to stop the MH algorithm for the first and second binding motifs.
* N: The number of sampling iterations to perform.

For `res_ana` function, the inputs are
*  motif_len_w: The length of the first binding motif.
*  motif_len_g: The length of the second binding motif.
*  N: The number of sampling iterations.



## Output
For  `Data_process` function, the output is
- a list containing:
  - `Data`: The updated data frame with W_obs and G_obs columns added.
  - `total_n`: The total number of sequences.
  - `dict`: A dictionary of all unique characters in the sequences.
  - `len_dict`: The number of unique characters in the dictionary.
  - `Len_seq`: The length of each sequence.
  - `UW_loc`: The indices where W_obs is NA (unobserved).
  - `UG_loc`: The indices where G_obs is NA (unobserved).

For `de_motif` function, the output is
* NULL. Results are saved to files in specified directories.

For `res_ana` function, the output is
- a list containing:
  - Samples for `theta_0` (background probabilities),
  - Samples for `theta` (the first binding motif probabilities),
  - Samples for `tilde{theta}` (the second binding motif probabilities),
  - Samples for `W` and `G` (indicators for motif presence),
  - Samples for `A` and `B` (positions of motifs),
  - Log-likelihood values,
  - Estimated parameters (`theta_0`, `theta`, `tilde{theta}`),
  - The number of jumps for both motifs,
  - a plot for the curve of the loglikelihood, the estimated theta logo and the estimated tilde{theta} logo,
  - a plot for the logos before and after a jump.



## Example
This is a toy example for decomposing hierarchical motifs. We shall use the data contained in the package.
```R
library(demotif)
# Set the proportion of missing data to 50%
rate_missg = 0.5
# Set the length of the first motif to 9
motif_len_w = 9
# Set the length of the second motif to 5
motif_len_g = 5
# Load example data from the 'demotif' package
# The dataset 'rate_missg_0.5_motif_len_g_5' corresponds to sequences with 50% missing data
data("rate_missg_0.5_motif_len_g_5")

#### Process Data ####
# Process the data using the `Data_process` function and unpack the result into individual variables
# The function  performs preprocessing and returns:
# - Data: Processed sequence data
# - total_n: Total number of sequences
# - dict: Dictionary of valid characters (e.g., A, C, G, T for DNA)
# - len_dict: Length of the dictionary
# - Len_seq: Lengths of the sequences
# - UW_loc: Indices of sequences with unknown W values
# - UG_loc: Indices of sequences with unknown G values
Data_row = data
Data_row_list = Data_process(Data_row)
names(Data_row_list) = c("Data", "total_n", "dict", "len_dict", "Len_seq", "UW_loc", "UG_loc")
list2env(Data_row_list, envir = .GlobalEnv)

##### Obtain the all samples #####
# Generate samples using the `de_motif` function (store them in the './result/temp_data/')
# Parameters:
# - Data: Processed data
# - motif_len_w: Length of the first motif
# - motif_len_g: Length of the second motif
# - 10: the number of steps to skip between each execution of the MH algorithm for the first binding motif.
# - 10: the number of steps to skip between each execution of the MH algorithm for the second binding motif. 
# - 20: the time step at which to stop the MH algorithm for the first and second binding motifs.
# - 50: The number of sampling iterations to perform.
de_motif(Data,motif_len_w,motif_len_g,10,10,20,50)

#### Present all results #####
# Analyze the results in the './result/temp_data/' using the `res_ana` function
# Parameters:
# - motif_len_w: Length of the first motif
# - motif_len_g: Length of the second motif
# - 20: the number of burn-in
# - 50: the number of iterations
# Return: 
#  Res: A list containing:
# - Samples for `theta_0` (background probabilities),
# - Samples for `theta` (the first binding motif probabilities),
# - Samples for `tilde{theta}` (the second binding motif probabilities),
# - Samples for `W` and `G` (indicators for motif presence),
# - Samples for `A` and `B` (positions of motifs),
# - Log-likelihood values,
# - Estimated parameters (`theta_0`, `theta`, `tilde{theta}`),
# - The number of jumps for both motifs,
# - a plot for the curve of the loglikelihood, the estimated theta logo and the estimated tilde{theta} logo,
# - a plot for the logos before and after a jump.
Res = res_ana(motif_len_w, motif_len_g, 20, 50)

#### Predict motif locations #####
# Perform predictive residue analysis using the `pred_res_ana` function.
# This function estimates motif locations and computes posterior probabilities.
#
# Parameters:
# - Data_row_list$Data: A list containing observed sequence data and related information.
# - Res[[9]]: Estimation for `theta_0` (background probabilities).
# - Res[[10]]: Estimation for `theta` (the first binding motif probabilities).
# - Res[[11]]: Estimation for `tilde{theta}` (the second binding motif probabilities).
# - motif_len_w: Length of the first motif (W).
# - motif_len_g: Length of the second motif (G).
# - 20: The number of burn-in iterations for MCMC sampling.
# - 50: The total number of MCMC iterations.
#
# Return:
# - pred_frame: A data frame containing:
#   - Estimated values for W and G,
#   - Posterior probabilities for motif presence,
#   - The most likely motif positions based on log-likelihood.
pred_frame = pred_res_ana(Data_row_list$Data, Res[[9]], Res[[10]], Res[[11]], motif_len_w, motif_len_g, 20, 50)
```

## Reference
-   Tang, X., & Liu, R. (2025). De-motif sampling: an approach to decompose hierarchical motifs with applications in T cell recognition. Briefings in Bioinformatics, 26(3), bbaf221. https://doi.org/10.1093/bib/bbaf221

## Contact
Xinyi Tang: xytang@link.cuhk.edu.hk; Ran Liu: ranliu@bnu.edu.cn
