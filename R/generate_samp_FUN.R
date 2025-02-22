#' Generate Sample Data for Motif Analysis
#'
#' This function generates sample motif data by iteratively updating 
#' the binding locations (A, B) and binding labels (W, G) using probabilistic sampling.
#'
#' @param Data_name A string indicating the dataset name.
#' @param Data_pre A data frame containing the input sequence (`seq`), observed `W` (`W_obs`), and observed `G` (`G_obs`).
#' @param dict A dictionary containing motif-related information.
#' @param theta_0 Estimated parameter for theta_0.
#' @param theta Estimated parameter for theta.
#' @param theta_til Estimated parameter for theta_til.
#' @param motif_len_w The length of motif W.
#' @param motif_len_g The length of motif G.
#' @param num_max The number of sampling iterations.
#' 
#' @return This function does not return a value but saves sampled data as CSV files in './result/temp_data'.
#' 
#' @export
#' 
generate_samp = function(Data_name, Data_pre, dict, theta_0, theta, theta_til, motif_len_w, motif_len_g, num_max){
  ####  Create directories for saving the results ####
  dir_name = './result'
  if(dir.exists(dir_name)==0){dir.create(dir_name)}
  dir_name = paste0(dir_name,'/temp_data')
  if(dir.exists(dir_name)==0){dir.create(dir_name)}
  dir_name = paste0(dir_name,'/',Data_name,'_motif_len_w=',motif_len_w,'_motif_len_g=',motif_len_g)
  if(dir.exists(dir_name)==0){dir.create(dir_name)}
  
  
  #### Initialize observed data and missing value indices ####
  R = Data_pre$seq
  W_obs = Data_pre$W_obs
  G_obs = Data_pre$G_obs
  
  # Identify missing values in W and G
  UG_loc = which(is.na(Data_pre$G_obs))
  UW_loc = which(is.na(Data_pre$W_obs))
  total_n = 1
  
  #### Set initial values ####
  ##prior parameters
  pri_wp = .5
  pri_gp = .5
  
  # Initialize W (filling missing values with sampled probabilities)
  W_ini = W_obs; 
  if(length(UW_loc)>0){
    W_ini[UW_loc] = sample(c(1,0),size = length(UW_loc),replace=T,prob = c(pri_wp,(1-pri_wp)))
  }
  
  # Initialize G (filling missing values with sampled probabilities)
  G_ini = G_obs; 
  if(length(UG_loc)>0){
    UG_loc = which(is.na(G_ini))
    G_ini[UG_loc] = sample(c(1,0),size = length(UG_loc),replace=T,prob = c(pri_gp,(1-pri_gp)))
  }
  
  # Initialize B (G motif start position)
  B_ini = rep(NA,1)
  for (i in 1:total_n) {
    B_ini[i] = sample(1:(Len_seq[i]-motif_len_g+1),size = 1,replace=T)
  }
  
  # Initialize A (W motif start position)
  A_ini = rep(NA,1)
  for (i in 1:total_n) {
    A_ini[i] = sample(1:(Len_seq[i]-motif_len_w+1),size = 1,replace=T)
  }
  
  #### Start iterative sampling process ####
  for (rep_N in 1:num_max) {
    if(rep_N==1){
      W_samp = W_ini
      G_samp = G_ini
      A_samp = A_ini
      B_samp = B_ini
    }else{
      if(length(UW_loc)>0){
        W_samp[UW_loc] = w_samp_fun(R,G_samp,which(is.na(Data_pre$W_obs)),A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0,theta,theta_til,pri_wp)[[1]]
      }
      if(length(UG_loc)>0){
        G_samp[UG_loc] = g_samp_fun(R,W_samp,which(is.na(Data_pre$G_obs)),A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0,theta,theta_til,pri_gp)[[1]]
      }
      A_samp = pred_A_samp_fun(R,W_samp,G_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0,theta,theta_til)[[1]]
      B_samp = pred_B_samp_fun(R,W_samp,G_samp,A_samp,motif_len_w,motif_len_g,dict,theta_0,theta,theta_til)[[1]]
      
    }
    
    logllk = logllk_fun(R,W_samp,G_samp,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0,theta,theta_til)
    
    write.csv(W_samp,file=paste0(dir_name,'/W_rep_N=',rep_N,'.csv'))
    write.csv(G_samp,file=paste0(dir_name,'/G_rep_N=',rep_N,'.csv'))
    write.csv(A_samp,file=paste0(dir_name,'/A_rep_N=',rep_N,'.csv'))
    write.csv(B_samp,file=paste0(dir_name,'/B_rep_N=',rep_N,'.csv'))
    write.csv(logllk,file=paste0(dir_name,'/logllk_rep_N=',rep_N,'.csv'))
  }
  return()
}


