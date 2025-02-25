#' Predictive Residue Analysis
#'
#' This function performs predictive residue analysis using an MCMC-based approach.
#'
#' @param Data A list containing observed sequence data and related information.
#' @param theta_0_samp A parameter matrix for theta_0.
#' @param theta_samp A parameter matrix for theta.
#' @param theta_til_samp A parameter matrix for theta_til.
#' @param motif_len_w The length of the W motif.
#' @param motif_len_g The length of the G motif.
#' @param burnin The number of burn-in iterations for MCMC sampling.
#' @param N The total number of MCMC iterations.
#'
#' @return A data frame containing estimated values for W and G along with posterior probabilities.
#' @export
#'
pred_res_ana = function(Data,theta_0_samp,theta_samp,theta_til_samp, motif_len_w,motif_len_g,burnin, N){
  ################################################
  ### Generate parameters and latent varibel ###
  ###############################################
  R = Data$seq
  W_obs = Data$W_obs
  G_obs = Data$G_obs
  UG_loc = which(is.na(G_obs))
  UW_loc = which(is.na(W_obs))
  #### initial value ###
  ##prior parameters
  pri_al0=1
  pri_alw=1
  pri_alg=1
  pri_wp = sum(Data$W[!is.na(Data$W)]==1)/sum(!is.na(Data$W))
  pri_gp = sum(Data$G[!is.na(Data$G)]==1)/sum(!is.na(Data$G))
  pri_wp = .5
  pri_gp = .5

  alpha_0 = rep(pri_al0,length(dict))
  alpha_j = rep(pri_alw,length(dict))
  alpha_j_til = rep(pri_alg,length(dict))
  
  
  
  W_ini = W_obs;
  if(length(UW_loc)>0){
  UW_loc = which(is.na(W_ini))
  W_ini[UW_loc] = sample(c(1,0),size = length(UW_loc),replace=T,prob = c(pri_wp,(1-pri_wp)))
  }
  
  G_ini = G_obs; 
  if(length(UG_loc)>0){
    UG_loc = which(is.na(G_ini))
    G_ini[UG_loc] = sample(c(1,0),size = length(UG_loc),replace=T,prob = c(pri_gp,(1-pri_gp)))
  }
  
  B_ini = rep(NA,total_n)
  for (i in 1:total_n) {
    B_ini[i] = sample(1:(Len_seq[i]-motif_len_g+1),size = 1,replace=T)
  }
  
  A_ini = rep(NA,total_n)
  for (i in 1:total_n) {
    A_ini[i] = sample(1:(Len_seq[i]-motif_len_w+1),size = 1,replace=T)
  }
  
  
  
  Logllk = rep(NA,N)
  W_SAMP = matrix(NA, nrow=total_n, ncol = N)
  G_SAMP = matrix(NA, nrow=total_n, ncol = N)
  A_SAMP = matrix(NA, nrow=total_n, ncol = N)
  B_SAMP = matrix(NA, nrow=total_n, ncol = N)

  for (rep_N in 1:N){
    # print(rep_N)
    if(rep_N==1){
      W_samp = W_ini
      G_samp = G_ini
      A_samp = A_ini
      B_samp = B_ini
      
      W_SAMP[,rep_N] = W_samp
      G_SAMP[,rep_N] = G_samp
      A_SAMP[,rep_N] = A_samp
      B_SAMP[,rep_N] = B_samp
    }else{
      if(length(UW_loc)>0){
      W_samp[UW_loc] = w_samp_fun(R,G_samp,UW_loc,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp,pri_wp)[[1]]
      W_SAMP[,rep_N] = W_samp
      }else{
        W_SAMP[,rep_N] = W_samp
      }
      if(length(UG_loc)>0){
        G_samp[UG_loc] = g_samp_fun(R,W_samp,UG_loc,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp,pri_gp)[[1]]
        G_SAMP[,rep_N] = G_samp
      }else{
        G_SAMP[,rep_N] = G_samp
      }
      A_samp = A_samp_fun(R,W_samp,G_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp)[[1]]
      A_SAMP[,rep_N] = A_samp
      B_samp = B_samp_fun(R,W_samp,G_samp,A_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp)[[1]]
      B_SAMP[,rep_N] = B_samp
    }
    
    Logllk[rep_N] = logllk_fun(R,W_samp,G_samp,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp)
  }
  
  #### Predict motif locations ####
  # Identify the iteration with the highest log-likelihood after the burn-in period
  max_loc = burnin  + which.max(Logllk[(burnin+1):N])
  # Extract the estimated W and G values from the best iteration
  est_W = W_SAMP[,max_loc];
  est_G = G_SAMP[,max_loc]
  # Compute posterior probabilities of W and G being 1 across MCMC samples
  prob_00 = rowSums((W_SAMP[, (burnin+1):N] == 0)&(G_SAMP[, (burnin+1):N] == 0))/length((burnin+1):N)
  prob_01 = rowSums((W_SAMP[, (burnin+1):N] == 0)&(G_SAMP[, (burnin+1):N] == 1))/length((burnin+1):N)
  prob_10 = rowSums((W_SAMP[, (burnin+1):N] == 1)&(G_SAMP[, (burnin+1):N] == 0))/length((burnin+1):N)
  prob_11 = rowSums((W_SAMP[, (burnin+1):N] == 1)&(G_SAMP[, (burnin+1):N] == 1))/length((burnin+1):N)
  
  res_frame = data.frame(Seq = Data$seq, est_W = est_W, est_G = est_G, prob_00 = prob_00, prob_01 = prob_01,prob_10 = prob_10,prob_11 = prob_11)
  
return(res_frame)
}
