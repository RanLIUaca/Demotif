#' de_motif
#'
#' A function to perform motif detection with Metropolis-Hastings (MH) sampling.
#'
#' @param Data A data frame containing the sequence, the label vector W, and the label vector G.
#' @param motif_len_w An integer specifying the length of the first binding motif.
#' @param motif_len_g An integer specifying the length of the second binding motif.
#' @param mh_step_w An integer specifying the number of steps to skip between each execution of the MH algorithm for the first binding motif.
#' @param mh_step_g An integer specifying the number of steps to skip between each execution of the MH algorithm for the second binding motif.
#' @param stop_jump_step An integer specifying the time step at which to stop the MH algorithm for the first and second binding motifs.
#' @param N The number of sampling iterations to perform.
#'
#' @import MCMCpack
#' @import stringr
#' @import MASS
#' @return NULL. Results are saved to files in specified directories.
#' @export
#'

de_motif = function(Data,motif_len_w,motif_len_g,mh_step_w,mh_step_g,stop_jump_step,N){
  #### Direction for saving ###
  dir_name = './result'
  if(dir.exists(dir_name)==0){dir.create(dir_name)}
  dir_name = paste0(dir_name,'/temp_data')
  if(dir.exists(dir_name)==0){dir.create(dir_name)}
  dir_name = paste0(dir_name,'/motif_len_w=',motif_len_w,'_motif_len_g=',motif_len_g)
  if(dir.exists(dir_name)==0){dir.create(dir_name)}
  
  change_name = paste0(dir_name,'/change')
  if(dir.exists(change_name)==0){dir.create(change_name)}
  
  
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
  # pri_wp = sum(Data$W[!is.na(Data$W)]==1)/sum(!is.na(Data$W))
  # pri_gp = sum(Data$G[!is.na(Data$G)]==1)/sum(!is.na(Data$G))
  pri_wp = .5
  pri_gp = .5
  # sum(Data$G[!is.na(Data$G)]==1)/sum(!is.na(Data$G))
  
  alpha_0 = rep(pri_al0,length(dict))
  alpha_j = rep(pri_alw,length(dict))
  alpha_j_til = rep(pri_alg,length(dict))
  
  # initial values
  theta_0_ini = c(rdirichlet(1,alpha_0))
  
  theta_ini = matrix(nrow = length(dict), ncol = motif_len_w)
  for (i in 1:motif_len_w) {
    temp = rdirichlet(1,alpha_j)
    theta_ini[,i] = temp
  }
  
  theta_til_ini = matrix(nrow = length(dict), ncol = motif_len_g)
  for (i in 1:motif_len_g) {
    temp = rdirichlet(1,alpha_j_til)
    theta_til_ini[,i] = temp
  }
  
  W_ini = Data$W_obs; 
  if(length(UW_loc)>0){
  UW_loc = which(is.na(W_ini))
  W_ini[UW_loc] = sample(c(1,0),size = length(UW_loc),replace=T,prob = c(pri_wp,(1-pri_wp)))
  }
  
  G_ini = Data$G_obs; 
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
  
  
  no_jump_w=0;no_jump_g=0;
  loc_jump_w=c(NA);loc_jump_g=c(NA)
  for (rep_N in 1:N) {
    # rep_N=2
    # print(rep_N)
    if(rep_N==1){
      theta_0_samp = theta_0_ini
      theta_samp = theta_ini
      theta_til_samp = theta_til_ini
      
      W_samp = W_ini
      G_samp = G_ini
      A_samp = A_ini
      B_samp = B_ini
    }else{
      theta_0_samp = theta_0_samp_fun(R,W_samp,G_samp,A_samp,B_samp,motif_len_w,motif_len_g,dict,pri_al0)
      theta_samp = theta_samp_fun(R,W_samp,G_samp,A_samp,B_samp,motif_len_w,motif_len_g,dict,pri_alw)[[1]]
      theta_til_samp = theta_til_samp_fun(R,W_samp,G_samp,B_samp,motif_len_g,dict,pri_alg)[[1]]
      
      temp_theta_til_samp = as.matrix(theta_til_samp)
      temp_theta_samp = as.matrix(theta_samp)
      row.names(temp_theta_til_samp) = dict
      row.names(temp_theta_samp) = dict
      if(length(UW_loc)>0){
      W_samp[UW_loc] = w_samp_fun(R,G_samp,UW_loc,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp,pri_wp)[[1]]
      }
      if(length(UG_loc)>0){
      G_samp[UG_loc] = g_samp_fun(R,W_samp,UG_loc,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp,pri_gp)[[1]]
      }
      A_samp = A_samp_fun(R,W_samp,G_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp)[[1]]
      B_samp = B_samp_fun(R,W_samp,G_samp,A_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp)[[1]]
      # temp_B_samp = B_samp
      
      if(((rep_N%%mh_step_w)==0)&(rep_N<stop_jump_step)){
        A_theta_samp = A_theta_samp_fun(R,W_samp,G_samp,UW_loc,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp,pri_alw,pri_wp)
        A_samp = A_theta_samp[[1]]
        theta_samp = A_theta_samp[[2]]
        if(length(UW_loc)>0){
        W_samp = A_theta_samp[[3]]
        }
        no_jump_w= no_jump_w +A_theta_samp[[4]]
        if(A_theta_samp[[4]]==1){
          loc_jump_w = rbind(loc_jump_w,rep_N)
          write.csv(temp_theta_samp,file=paste0(dir_name,'/change_theta_pre_rep_N=',rep_N,'.csv'))
          theta_post = as.data.frame(theta_samp); row.names(theta_post) = dict
          write.csv(theta_post,file=paste0(dir_name,'/change_theta_post_rep_N=',rep_N,'.csv'))
          
        }
      }
      
      
      if(((rep_N%%mh_step_g)==0)&(rep_N<stop_jump_step)){
        B_theta_til_samp = B_theta_til_samp_fun(R,W_samp,G_samp,UG_loc,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp,pri_alg,pri_gp)
        B_samp = B_theta_til_samp[[1]]
        theta_til_samp = B_theta_til_samp[[2]]
        if(length(UW_loc)>0){
        G_samp = B_theta_til_samp[[3]]
        }
        no_jump_g= no_jump_g +B_theta_til_samp[[4]]
        if(B_theta_til_samp[[4]]==1){
          loc_jump_g = rbind(loc_jump_g,rep_N)
          write.csv(temp_theta_til_samp,file=paste0(dir_name,'/change_theta_til_pre_rep_N=',rep_N,'.csv'))
          theta_til_post = as.data.frame(theta_til_samp); row.names(theta_til_post) = dict
          write.csv(theta_til_post,file=paste0(dir_name,'/change_theta_til_post_rep_N=',rep_N,'.csv'))
        } 
      }
    }
    
    logllk = logllk_fun(R,W_samp,G_samp,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp)
    
    write.csv(theta_0_samp,file=paste0(dir_name,'/theta_0_rep_N=',rep_N,'.csv'))
    write.csv(theta_samp,file=paste0(dir_name,'/theta_rep_N=',rep_N,'.csv'))
    write.csv(theta_til_samp,file=paste0(dir_name,'/theta_til_rep_N=',rep_N,'.csv'))
    
    write.csv(W_samp,file=paste0(dir_name,'/W_rep_N=',rep_N,'.csv'))
    write.csv(G_samp,file=paste0(dir_name,'/G_rep_N=',rep_N,'.csv'))
    write.csv(A_samp,file=paste0(dir_name,'/A_rep_N=',rep_N,'.csv'))
    write.csv(B_samp,file=paste0(dir_name,'/B_rep_N=',rep_N,'.csv'))
    write.csv(logllk,file=paste0(dir_name,'/logllk_rep_N=',rep_N,'.csv'))
  }
  # print(c(no_jump_w,no_jump_g))
  num_change = data.frame(no_jump_w=no_jump_w, no_jump_g=no_jump_g)
  Loc_jump_w = data.frame(Loc_jump_w=loc_jump_w);Loc_jump_g = data.frame(Loc_jump_g=loc_jump_g)
  write.csv(num_change,file=paste0(dir_name,'/num_change','.csv'))
  write.csv(Loc_jump_w,file=paste0(dir_name,'/Loc_jump_w','.csv'))
  write.csv(Loc_jump_g,file=paste0(dir_name,'/Loc_jump_g','.csv'))
  return()
}


