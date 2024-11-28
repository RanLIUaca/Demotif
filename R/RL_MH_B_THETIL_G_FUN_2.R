#' pi_B_theta_til_G_fun
#'
#' This function calculates the posterior probability for a given configuration of `B`, `theta_til`, and `G`.
#'
#' @param R A vector representing the observed sequences.
#' @param W_samp A vector indicating the presence of the first motif (sampled).
#' @param G A vector indicating the presence of the second motif.
#' @param UG_loc Locations where `G` is unknown.
#' @param A_samp A vector of positions for the first motif (sampled).
#' @param B A vector of positions for the second motif.
#' @param dict A dictionary of possible characters in the sequences.
#' @param theta_0_samp Background probabilities for each character.
#' @param theta_samp Probabilities for the first motif (sampled).
#' @param theta_til Probabilities for the second motif.
#' @param pri_alg Prior Dirichlet parameter for `theta_til`.
#' @param pri_gp Prior probability for `G`.
#'
#' @return A numeric value representing the posterior probability.
#' @export
#' 
pi_B_theta_til_G_fun = function(R,W_samp,G,UG_loc,A_samp,B,dict,theta_0_samp,theta_samp,theta_til,pri_alg,pri_gp){
  res = logllk_fun(R,W_samp,G,A_samp,B,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til)
  for(rep_j in 1:length(theta_til[1,])){
    res = res +log(ddirichlet(theta_til[,rep_j], rep(pri_alg,len_dict)))
  }
  for (rep_i in UG_loc) {
    res = res + log(dbinom(G[rep_i], size = 1, prob = pri_gp))
  }
  return(res)
}


#' prob_theta_til_fun
#'
#' This function calculates the posterior probability of `theta_til` based on the observed sequences and other parameters.
#'
#' @param R A vector representing the observed sequences.
#' @param W_samp A vector indicating the presence of the first motif (sampled).
#' @param G A vector indicating the presence of the second motif.
#' @param A_samp A vector of positions for the first motif (sampled).
#' @param B A vector of positions for the second motif.
#' @param dict A dictionary of possible characters in the sequences.
#' @param theta_til Probabilities for the second motif.
#' @param pri_alg Prior Dirichlet parameter for `theta_til`.
#'
#' @return A numeric value representing the posterior probability of `theta_til`.
#' @export
#' 
prob_theta_til_fun = function(R,W_samp,G,A_samp,B,dict,theta_til,pri_alg){
  # sample parameter
  Post_alpha_g = matrix(nrow = length(dict), ncol = motif_len_g)
  for (rep_j in 1:motif_len_g) {
    #rep_j=1
    # post_parameter
    H_b=rep(0,length(dict))
    names(H_b)=dict
    for (rep_i in 1:total_n) {
      #rep_i=1
      Ri = R[rep_i]
      tmp_index = 1:str_length(Ri)
      Ri = str_sub(Ri,tmp_index,tmp_index)
      
      bij = B[rep_i]+rep_j-1
      
      
      if(((W_samp[rep_i]==1)&(G[rep_i]==1))||((W_samp[rep_i]==0)&(G[rep_i]==1))){
        Ri_bij = Ri[bij]
        tmp_index = 1:str_length(Ri_bij)
        Ri_bij = str_sub(Ri_bij,tmp_index,tmp_index)
        Ri_bij = factor(Ri_bij, levels = dict)
        H_b=H_b+table(Ri_bij)
      }
    }
    
    post_alpha_g = pri_alg+H_b
    Post_alpha_g[,rep_j] = post_alpha_g
  }
  
  res = 0
  for (rep_j in 1:length(theta_til[1,])) {
    res = res + log(ddirichlet(theta_til[,rep_j], Post_alpha_g[,rep_j]))
    # print(res)
    # print(Post_alpha_g[,rep_j])
    # print(sum(Post_alpha_g[,rep_j]))
    
  }
  return(res)
}

#' prob_g_fun
#'
#' This function calculates the posterior probability of `G` (presence of the second motif).
#'
#' @param R A vector representing the observed sequences.
#' @param W_samp A vector indicating the presence of the first motif (sampled).
#' @param G A vector indicating the presence of the second motif.
#' @param UG_loc Locations where `G` is unknown.
#' @param A_samp A vector of positions for the first motif (sampled).
#' @param B A vector of positions for the second motif.
#' @param dict A dictionary of possible characters in the sequences.
#' @param theta_0_samp Background probabilities for each character.
#' @param theta_samp Probabilities for the first motif (sampled).
#' @param theta_til Probabilities for the second motif.
#' @param pri_gp Prior probability for `G`.
#'
#' @return A numeric value representing the posterior probability of `G`.
#' @export
#' 
prob_g_fun = function(R,W_samp,G,UG_loc,A_samp,B,dict,theta_0_samp,theta_samp,theta_til,pri_gp){
  # Post_prob = g_samp_fun(R,W_samp,UG_loc,A_samp,B,dict,theta_0_samp,theta_samp,theta_til,pri_gp)[[2]]
  Post_prob = matrix(NA, nrow = total_n, ncol = 2)
  for (rep_i in UG_loc) {
    # rep_i=4
    Rui = R[rep_i]
    tmp_index = 1:str_length(Rui)
    Rui = str_sub(Rui,tmp_index,tmp_index)
    
    wui = W_samp[rep_i]
    aui = (A_samp[rep_i]:(A_samp[rep_i]+motif_len_w-1))
    bui = (B[rep_i]:(B[rep_i]+motif_len_g-1))
    llk_one_seq_g1 = exp(logllk_fun_one_seq(Rui,wui,1,aui,bui,dict,theta_0_samp,theta_samp,theta_til))
    llk_one_seq_g0 = exp(logllk_fun_one_seq(Rui,wui,0,aui,bui,dict,theta_0_samp,theta_samp,theta_til))
    
    post_gp = (pri_gp*llk_one_seq_g1)/(pri_gp*llk_one_seq_g1+(1-pri_gp)*llk_one_seq_g0)
    Post_prob[rep_i,] = c(post_gp,(1-post_gp))
    # print(post_gp)
  }
  
  res = 0
  for (rep_i in UG_loc) {
    res = res + log(dbinom(G[rep_i], size = 1, prob = Post_prob[rep_i,1]))
  }
  return(res)
}

#' B_theta_til_samp_fun
#'
#' This function performs Metropolis-Hastings sampling for the positions of the second motif (`B`), its probabilities (`theta_til`), and `G`.
#'
#' @param R A vector representing the observed sequences.
#' @param W_samp A vector indicating the presence of the first motif (sampled).
#' @param G_samp A vector indicating the presence of the second motif (sampled).
#' @param UG_loc Locations where `G` is unknown.
#' @param A_samp A vector of positions for the first motif (sampled).
#' @param B_samp A vector of positions for the second motif (sampled).
#' @param motif_len_w Length of the first motif.
#' @param motif_len_g Length of the second motif.
#' @param dict A dictionary of possible characters in the sequences.
#' @param theta_0_samp Background probabilities for each character.
#' @param theta_samp Probabilities for the first motif (sampled).
#' @param theta_til_samp Probabilities for the second motif (sampled).
#' @param pri_alg Prior Dirichlet parameter for `theta_til`.
#' @param pri_gp Prior probability for `G`.
#'
#' @return A list containing updated values for `B`, `theta_til`, and `G`, along with the jump indicator.
#' @export
#' 
B_theta_til_samp_fun = function(R,W_samp,G_samp,UG_loc,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp,pri_alg,pri_gp){
  del = sample(c(1,-1),size = 1,replace=T,prob = c(1/2,1/2))
  B_samp_star = B_samp + del
  if(del==-1) B_samp_star[which(B_samp==1)] = B_samp[which(B_samp==1)]
  if(del==1) B_samp_star[which((B_samp+motif_len_g-1) == Len_seq)] = B_samp[which((B_samp+motif_len_g-1) == Len_seq)]
 
  res_theta_til = theta_til_samp_fun(R,W_samp,G_samp,B_samp_star,motif_len_g,dict,pri_alg)
  theta_til_samp_star = res_theta_til[[1]]
  G_samp_star = G_samp
  
  if(length(UG_loc)>0){
  res_g = g_samp_fun(R,W_samp,UG_loc,A_samp,B_samp_star,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp_star,pri_gp)
  G_samp_star[UG_loc] = res_g[[1]]
  
  
  nume = (pi_B_theta_til_G_fun(R,W_samp,G_samp_star,UG_loc,A_samp,B_samp_star,dict,theta_0_samp,theta_samp,theta_til_samp_star,pri_alg,pri_gp)+
            prob_theta_til_fun(R,W_samp,G_samp_star,A_samp,B_samp,dict,theta_til_samp,pri_alg)+
            prob_g_fun(R,W_samp,G_samp,UG_loc,A_samp,B_samp,dict,theta_0_samp,theta_samp,theta_til_samp,pri_gp))
 
  
  deno = (pi_B_theta_til_G_fun(R,W_samp,G_samp,UG_loc,A_samp,B_samp,dict,theta_0_samp,theta_samp,theta_til_samp,pri_alg,pri_gp)+
            prob_theta_til_fun(R,W_samp,G_samp,A_samp,B_samp_star,dict,theta_til_samp_star,pri_alg)+
            prob_g_fun(R,W_samp,G_samp_star,UG_loc,A_samp,B_samp_star,dict,theta_0_samp,theta_samp,theta_til_samp_star,pri_gp))
  }else{
    nume = (pi_B_theta_til_G_fun(R,W_samp,G_samp_star,UG_loc,A_samp,B_samp_star,dict,theta_0_samp,theta_samp,theta_til_samp_star,pri_alg,pri_gp)+
              prob_theta_til_fun(R,W_samp,G_samp_star,A_samp,B_samp,dict,theta_til_samp,pri_alg))
    
    deno = (pi_B_theta_til_G_fun(R,W_samp,G_samp,UG_loc,A_samp,B_samp,dict,theta_0_samp,theta_samp,theta_til_samp,pri_alg,pri_gp)+
              prob_theta_til_fun(R,W_samp,G_samp,A_samp,B_samp_star,dict,theta_til_samp_star,pri_alg))
  }
  acc_rate = min(1, exp(nume-deno))
  # print(acc_rate)
  acc_rn = runif(1)
  if((acc_rn<acc_rate)&(nume!=deno)){
    # print(c('TCR_jump',acc_rn,acc_rate));
    no_jump=1;
    return(list(B_samp_star,theta_til_samp_star,G_samp_star,no_jump))
  }else{no_jump=0;return(list(B_samp,theta_til_samp,G_samp,no_jump))}
}