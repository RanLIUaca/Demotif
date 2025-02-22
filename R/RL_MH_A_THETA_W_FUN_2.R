#' pi_A_theta_fun
#'
#' This function calculates the posterior probability for a given configuration of `A`, `theta`, and other parameters.
#'
#' @param R A vector representing the observed sequences.
#' @param W A vector indicating the presence of the first motif.
#' @param G_samp A vector indicating the presence of the second motif.
#' @param A A vector of positions for the first motif.
#' @param B_samp A vector of positions for the second motif.
#' @param dict A dictionary of possible characters in the sequences.
#' @param theta_0_samp Background probabilities for each character.
#' @param theta Probabilities for the first motif.
#' @param theta_til_samp Probabilities for the second motif.
#' @param pri_alw Prior Dirichlet parameter for `theta`.
#'
#' @return A numeric value representing the posterior probability.
#' @export
#' 
pi_A_theta_fun = function(R,W,G_samp,A,B_samp,dict,theta_0_samp,theta,theta_til_samp,pri_alw){
  res = logllk_fun(R,W,G_samp,A,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta,theta_til_samp)
  for(rep_j in 1:length(theta[1,])){
    res = res +log(ddirichlet(theta[,rep_j], rep(pri_alw,len_dict)))
  }
  for (rep_i in UW_loc) {
    res = res + log(dbinom(W[rep_i], size = 1, prob = pri_Wp))
  }
  return(res)
}

#' prob_theta_fun
#'
#' This function calculates the posterior probability of `theta` based on the observed sequences and other parameters.
#'
#' @param R A vector representing the observed sequences.
#' @param W A vector indicating the presence of the first motif.
#' @param G_samp A vector indicating the presence of the second motif.
#' @param A A vector of positions for the first motif.
#' @param B_samp A vector of positions for the second motif.
#' @param dict A dictionary of possible characters in the sequences.
#' @param theta Probabilities for the first motif.
#' @param pri_alw Prior Dirichlet parameter for `theta`.
#'
#' @return A numeric value representing the posterior probability of `theta`.
#' @export
#' 
prob_theta_fun = function(R,W,G_samp,A,B_samp,dict,theta,pri_alw){
  # Post_alpha_w = theta_samp_fun(R,W,G_samp,A,B_samp,dict,pri_alw)[[2]]
  # sample parameter
  Post_alpha_w = matrix(nrow = length(dict), ncol = motif_len_w)
  for (rep_j in 1:motif_len_w) {
    #rep_j=1
    # post_parameter
    H_a=rep(0,length(dict))
    names(H_a)=dict
    for (rep_i in 1:total_n) {
      #rep_i=1
      Ri = R[rep_i]
      tmp_index = 1:str_length(Ri)
      Ri = str_sub(Ri,tmp_index,tmp_index)
      
      aij = A[rep_i]+rep_j-1;
      bi = (B_samp[rep_i]:(B_samp[rep_i]+motif_len_g-1))
      if(((W[rep_i]==1)&(G_samp[rep_i]==1))){
        aij = setdiff(aij,bi)
      }
      
      if(((W[rep_i]==1)&(G_samp[rep_i]==0))||((W[rep_i]==1)&(G_samp[rep_i]==1))){
        if(length(aij)==1){
          Ri_aij = Ri[aij]
          tmp_index = 1:str_length(Ri_aij)
          Ri_aij = str_sub(Ri_aij,tmp_index,tmp_index)
          Ri_aij = factor(Ri_aij, levels = dict)
          H_a=H_a+table(Ri_aij)
        }
      }
    }
    
    post_alpha_w = pri_alw+H_a
    Post_alpha_w[,rep_j] = post_alpha_w
  }
  
  res = 0
  for (rep_j in 1:length(theta[1,])) {
    res = res + log(ddirichlet(theta[,rep_j], Post_alpha_w[,rep_j]))
  }
  return(res)
}

#' prob_w_fun
#'
#' This function calculates the posterior probability of `W` (presence of the first motif).
#'
#' @param R A vector representing the observed sequences.
#' @param W A vector indicating the presence of the first motif.
#' @param G_samp A vector indicating the presence of the second motif.
#' @param UW_loc Locations where `W` is unknown.
#' @param A A vector of positions for the first motif.
#' @param B_samp A vector of positions for the second motif.
#' @param dict A dictionary of possible characters in the sequences.
#' @param theta_0_samp Background probabilities for each character.
#' @param theta Probabilities for the first motif.
#' @param theta_til_samp Probabilities for the second motif.
#' @param pri_wp Prior probability for `W`.
#'
#' @return A numeric value representing the posterior probability of `W`.
#' @export
#'  
prob_w_fun = function(R,W,G_samp,UW_loc,A,B_samp,dict,theta_0_samp,theta,theta_til_samp,pri_wp){
  # Post_prob = g_samp_fun(R,W_samp,UW_loc,A_samp,B,dict,theta_0_samp,theta_samp,theta_til,pri_gp)[[2]]
  Post_prob = matrix(NA, nrow = total_n, ncol = 2)
  for (rep_i in UW_loc) {
    # rep_i=4
    Rui = R[rep_i]
    tmp_index = 1:str_length(Rui)
    Rui = str_sub(Rui,tmp_index,tmp_index)
    
    gui = G_samp[rep_i]
    aui = (A[rep_i]:(A[rep_i]+motif_len_w-1))
    bui = (B_samp[rep_i]:(B_samp[rep_i]+motif_len_g-1))
    llk_one_seq_w1 = exp(logllk_fun_one_seq(Rui,1,gui,aui,bui,dict,theta_0_samp,theta,theta_til_samp))
    llk_one_seq_w0 = exp(logllk_fun_one_seq(Rui,0,gui,aui,bui,dict,theta_0_samp,theta,theta_til_samp))
    
    post_wp = (pri_wp*llk_one_seq_w1)/(pri_wp*llk_one_seq_w1+(1-pri_wp)*llk_one_seq_w0)
    Post_prob[rep_i,] = c(post_wp,(1-post_wp))
    # print(post_gp)
  }
  
  res = 0
  for (rep_i in UW_loc) {
    res = res + log(dbinom(W[rep_i], size = 1, prob = Post_prob[rep_i,1]))
  }
  return(res)
}

#' A_theta_samp_fun
#'
#' This function performs Metropolis-Hastings sampling for the positions of the first motif (`A`) and its probabilities (`theta`).
#'
#' @param R A vector representing the observed sequences.
#' @param W_samp A vector indicating the presence of the first motif (sampled).
#' @param G_samp A vector indicating the presence of the second motif (sampled).
#' @param UW_loc Locations where `W` is unknown.
#' @param A_samp A vector of positions for the first motif (sampled).
#' @param B_samp A vector of positions for the second motif (sampled).
#' @param motif_len_w Length of the first motif.
#' @param motif_len_g Length of the second motif.
#' @param dict A dictionary of possible characters in the sequences.
#' @param theta_0_samp Background probabilities for each character.
#' @param theta_samp Probabilities for the first motif (sampled).
#' @param theta_til_samp Probabilities for the second motif (sampled).
#' @param pri_alw Prior Dirichlet parameter for `theta`.
#' @param pri_wp Prior probability for `W`.
#'
#' @return A list containing updated values for `A`, `theta`, and `W`, along with the jump indicator.
#' @export
#' 
A_theta_samp_fun = function(R,W_samp,G_samp,UW_loc,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp,pri_alw,pri_wp){
  del = sample(c(1,-1),size = 1,replace=T,prob = c(1/2,1/2))
  A_samp_star = A_samp + del
  if(del==-1) A_samp_star[which(A_samp==1)] = A_samp[which(A_samp==1)]
  if(del==1) A_samp_star[which((A_samp+motif_len_w-1) == Len_seq)] = A_samp[which((A_samp+motif_len_w-1) == Len_seq)]
  
  res_theta = theta_samp_fun(R,W_samp,G_samp,A_samp_star,B_samp,motif_len_w,motif_len_g,dict,pri_alw)
  theta_samp_star = res_theta[[1]]

  W_samp_star = W_samp
  if(length(UW_loc)>0){
  res_w = w_samp_fun(R,G_samp,UW_loc,A_samp_star,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp_star,theta_til_samp,pri_wp)
  W_samp_star[UW_loc] = res_w[[1]]
  }
  
  if(length(UW_loc)>0){
  nume = (pi_A_theta_fun(R,W_samp_star,G_samp,A_samp_star,B_samp,dict,theta_0_samp,theta_samp_star,theta_til_samp,pri_alw)+
            prob_theta_fun(R,W_samp_star,G_samp,A_samp,B_samp,dict,theta_samp,pri_alw)+
            prob_w_fun(R,W_samp,G_samp,UW_loc,A_samp,B_samp,dict,theta_0_samp,theta_samp,theta_til_samp,pri_wp))
  
  deno = (pi_A_theta_fun(R,W_samp,G_samp,A_samp,B_samp,dict,theta_0_samp,theta_samp,theta_til_samp,pri_alw)+
            prob_theta_fun(R,W_samp,G_samp,A_samp_star,B_samp,dict,theta_samp_star,pri_alw)+
            prob_w_fun(R,W_samp_star,G_samp,UW_loc,A_samp_star,B_samp,dict,theta_0_samp,theta_samp_star,theta_til_samp,pri_wp))
  }else{
    nume = (pi_A_theta_fun(R,W_samp_star,G_samp,A_samp_star,B_samp,dict,theta_0_samp,theta_samp_star,theta_til_samp,pri_alw)+
              prob_theta_fun(R,W_samp_star,G_samp,A_samp,B_samp,dict,theta_samp,pri_alw))
    deno = (pi_A_theta_fun(R,W_samp,G_samp,A_samp,B_samp,dict,theta_0_samp,theta_samp,theta_til_samp,pri_alw)+
              prob_theta_fun(R,W_samp,G_samp,A_samp_star,B_samp,dict,theta_samp_star,pri_alw))
            }
  acc_rate = min(1, exp(nume-deno))
  acc_rn = runif(1)
  if((acc_rn<acc_rate)&(nume!=deno)){
    # print(c('MHC_jump',acc_rn,acc_rate));
    no_jump_w=1;
    return(list(A_samp_star,theta_samp_star,W_samp_star,no_jump_w))
  }else{no_jump_w=0;return(list(A_samp,theta_samp,W_samp,no_jump_w))}
}
