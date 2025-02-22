#' Predict Motif 
#'
#' This function predicts motif locations in a given sequence using MCMC sampling.
#' It processes the input data, generates samples using `generate_samp()`, 
#' reads the sampled data, and determines the most probable motif locations.
#'
#' @param Data_name A string representing the dataset name.
#' @param Data_one_row A dataframe containing one row of sequence data with columns: 
#'   1. `seq` (sequence string)
#'   2. `W_obs` (observed W values, may contain NAs)
#'   3. `G_obs` (observed G values, may contain NAs)
#' @param dict A dictionary containing motif-related information.
#' @param theta_0 Estimated parameter for theta_0.
#' @param theta Estimated parameter for theta.
#' @param theta_til Estimated parameter for theta_til.
#' @param motif_len_w Length of motif W.
#' @param motif_len_g Length of motif G.
#' @param burnin Number of initial iterations to discard as burn-in.
#' @param num_max Total number of MCMC iterations.
#' 
#' @return A list containing:
#'   - `res_frame`: A dataframe with predicted motif states and probabilities.
#'   - `est_A`: Estimated start and end positions for motif W.
#'   - `est_B`: Estimated start and end positions for motif G.
#'
#' @export
predict_demotif = function(Data_name, Data_one_row, dict,  theta_0, theta, theta_til, motif_len_w,motif_len_g, burnin, num_max){
  ### read all samples ###
  dir_name = './result'
  if(dir.exists(dir_name)==0){dir.create(dir_name)}
  dir_name = paste0(dir_name,'/temp_data')
  if(dir.exists(dir_name)==0){dir.create(dir_name)}
  model_name = paste0(dir_name,'/',Data_name,'_motif_len_w=',motif_len_w,'_motif_len_g=',motif_len_g)
  if(dir.exists(dir_name)==0){dir.create(dir_name)}
  
  ## Process input data ####
  len_dict = length(dict)
  Data_pre = data.frame(seq=NA, W_obs = NA, G_obs = NA)
  Data_pre$seq = Data_one_row[,1]; Data_pre$W_obs = Data_one_row[,2]; Data_pre$G_obs = Data_one_row[,3]
  Len_seq = str_length(Data_pre$seq[1])
  UG_loc = which(is.na(Data_pre$G_obs))
  UW_loc = which(is.na(Data_pre$W_obs))
  total_n = 1
  
  
  ##### Generate MCMC samples and save them in './result/temp_data' ####
  generate_samp(Data_name, Data_pre, dict, theta_0, theta, theta_til, motif_len_w, motif_len_g, num_max)
  
  ###### Read MCMC samples from saved files ####
  Logllk = rep(NA,num_max)
  W_samp = rep(NA,num_max)
  G_samp = rep(NA,num_max)
  A_samp = rep(NA,num_max)
  B_samp = rep(NA,num_max)
  
  for (rep_N in 1:num_max){
    #rep_N=1
    # print(rep_N)
    Logllk[rep_N] = as.numeric(read.csv(paste0(model_name,'/logllk_rep_N=',rep_N,'.csv'),header = T)[1,2])
    
    if(length(UW_loc)>0){
      W_samp[rep_N] = as.numeric(read.csv(paste0(model_name,'/W_rep_N=',rep_N,'.csv'),header = T)[1,2])
    }
    if(length(UG_loc)>0){
      G_samp[rep_N] = as.numeric(read.csv(paste0(model_name,'/G_rep_N=',rep_N,'.csv'),header = T)[1,2])
    }
    A_samp[rep_N] = as.matrix(read.csv(paste0(model_name,'/A_rep_N=',rep_N,'.csv'),header = T)[1,2])
    B_samp[rep_N] = as.matrix(read.csv(paste0(model_name,'/B_rep_N=',rep_N,'.csv'),header = T)[1,2])
  }
  
  #### Predict motif locations ####
  # Identify the iteration with the highest log-likelihood after the burn-in period
  max_loc = burnin  + which.max(Logllk[(burnin+1):num_max])
  # Extract the estimated W and G values from the best iteration
  est_W = W_samp[max_loc]; est_G = G_samp[max_loc]
  # Estimate motif start and end positions
  est_A = A_samp[max_loc]:(A_samp[max_loc]+motif_len_w-1); 
  est_B = B_samp[max_loc]:(B_samp[max_loc]+motif_len_g-1)
  # Compute posterior probabilities of W and G being 1 across MCMC samples
  prob_w = sum(W_samp[(burnin+1):num_max]==1)/length((burnin+1):num_max)
  prob_g = sum(G_samp[(burnin+1):num_max]==1)/length((burnin+1):num_max)
  res_frame = data.frame(Seq = Data_pre$seq, est_W = est_W, est_G = est_G, prob_w = prob_w, prob_g = prob_g)
  
  return(list(res_frame, est_A, est_B))
}