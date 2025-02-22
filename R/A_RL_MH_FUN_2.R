#' A_samp_fun
#'
#' A function to sample the position of the first binding motif using posterior probabilities.
#'
#' @param R A vector of sequences.
#' @param W_samp A vector of sampled W labels.
#' @param G_samp A vector of sampled G labels.
#' @param B_samp A vector of sampled positions for the second motif.
#' @param motif_len_w An integer specifying the length of the first binding motif.
#' @param motif_len_g An integer specifying the length of the second binding motif.
#' @param dict A dictionary of possible motif elements.
#' @param theta_0_samp A vector of sampled background probabilities.
#' @param theta_samp A matrix of sampled probabilities for the first motif.
#' @param theta_til_samp A matrix of sampled probabilities for the second motif.
#'
#' @return A list containing the sampled positions of the first binding motif.
#' @export
#'

A_samp_fun = function(R,W_samp,G_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp){
  A_samp_unobs = rep(NA, total_n)
  for (rep_i in 1:total_n) {
    # rep_i=1
    Ri = R[rep_i]
    len_seq = str_length(Ri)
    tmp_index = 1:str_length(Ri)
    Ri = str_sub(Ri,tmp_index,tmp_index)
    bi = (B_samp[rep_i]:(B_samp[rep_i]+motif_len_g-1))
    
    A_potential = matrix(NA,nrow = (len_seq-motif_len_w+1), ncol = motif_len_w)
    for (rep_Ai in 1:(len_seq-motif_len_w+1)) {
      A_potential[rep_Ai,] = seq(rep_Ai,length.out = motif_len_w, by=1)
    }
    
    res=rep(1,(len_seq-motif_len_w+1))
    for (rep_Ai in 1:(len_seq-motif_len_w+1)) {
      if((W_samp[rep_i]==1)&(G_samp[rep_i]==0)){
        res[rep_Ai]=exp(logllk_fun_10(Ri,A_potential[rep_Ai,],theta_0_samp,theta_samp,dict))
      }else if((W_samp[rep_i]==1)&(G_samp[rep_i]==1)){
        res[rep_Ai]=exp(logllk_fun_11(Ri,A_potential[rep_Ai,],bi,theta_0_samp,theta_samp,theta_til_samp,dict))
      }
    }
    
    post_prob=rep(NA,(len_seq-motif_len_w+1))
    for (rep_Ai in 1:(len_seq-motif_len_w+1)) {
      post_prob[rep_Ai] = res[rep_Ai]/sum(res)
    }
    
    A_samp_unobs[rep_i] = sample(1:(len_seq-motif_len_w+1),size = 1, prob = post_prob)
    # print(post_prob)
  }
  return(list(A_samp_unobs))
}
