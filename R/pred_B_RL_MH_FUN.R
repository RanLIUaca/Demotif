#' pred_B_samp_fun
#'
#' A function to sample the position of the second binding motif using posterior probabilities.
#'
#' @param R A sequence.
#' @param W_samp A sampled W label.
#' @param G_samp A sampled G label.
#' @param A_samp A sampled position for the first motif.
#' @param motif_len_w An integer specifying the length of the first binding motif.
#' @param motif_len_g An integer specifying the length of the second binding motif.
#' @param dict A dictionary of possible motif elements.
#' @param theta_0_samp A vector of sampled background probabilities.
#' @param theta_samp A matrix of sampled probabilities for the first motif.
#' @param theta_til_samp A matrix of sampled probabilities for the second motif.
#'
#' @return A list containing the sampled positions of the second binding motif.
#' @export
#'
pred_B_samp_fun = function(R,W_samp,G_samp,A_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp){
  B_samp_unobs = rep(NA, 1)
  # B_loc_last = rep(NA, total_n)
  for (rep_i in 1:1) {
    # rep_i=1
    Ri = R[rep_i]
    len_seq = str_length(Ri)
    tmp_index = 1:str_length(Ri)
    Ri = str_sub(Ri,tmp_index,tmp_index)
    ai = (A_samp[rep_i]:(A_samp[rep_i]+motif_len_w-1))
    
    B_potential = matrix(NA,nrow = (len_seq-motif_len_g+1), ncol = motif_len_g)
    for (rep_Bi in 1:(len_seq-motif_len_g+1)) {
      B_potential[rep_Bi,] = seq(rep_Bi,length.out = motif_len_g, by=1)
    }
    
    res=rep(1,(len_seq-motif_len_g+1))
    for (rep_Bi in 1:(len_seq-motif_len_g+1)) {
      if((W_samp[rep_i]==0)&(G_samp[rep_i]==1)){
        res[rep_Bi]=exp(logllk_fun_01(Ri,B_potential[rep_Bi,],theta_0_samp,theta_til_samp,dict))
      }else if((W_samp[rep_i]==1)&(G_samp[rep_i]==1)){
        res[rep_Bi]=exp(logllk_fun_11(Ri,ai,B_potential[rep_Bi,],theta_0_samp,theta_samp,theta_til_samp,dict))
      }
    }

    post_prob=rep(NA,(len_seq-motif_len_g+1))
    for (rep_Bi in 1:(len_seq-motif_len_g+1)) {
      post_prob[rep_Bi] = res[rep_Bi]/sum(res)
    }
   # print(post_prob)
    
    B_samp_unobs[rep_i] = sample(1:(len_seq-motif_len_g+1),size = 1, prob = post_prob)
    # B_loc_last[rep_i] = B_samp_unobs[rep_i,length(B_samp_unobs[rep_i,])]
    # names(B_loc_last) = 'B_loc_last'
  }
  # print(table(B_loc_last))
  return(list(B_samp_unobs))
}