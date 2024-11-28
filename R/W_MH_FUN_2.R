#' w_samp_fun
#'
#' @param R A vector of string sequences.
#' @param G_samp A vector indicating the presence of the second motif (sampled).
#' @param UW_loc Indices of sequences with unknown `W` values (unobserved presence of the first motif).
#' @param A_samp A vector of positions for the first motif (sampled).
#' @param B_samp A vector of positions for the second motif (sampled).
#' @param motif_len_w Length of the first motif.
#' @param motif_len_g Length of the second motif.
#' @param dict A dictionary of possible characters in the sequences.
#' @param theta_0_samp Background probabilities.
#' @param theta_samp Probabilities for the first motif.
#' @param theta_til_samp Probabilities for the second motif.
#' @param pri_wp Prior probability of the presence of the first motif (`W = 1`).
#'
#' @return A list containing:
#'   - `W_samp_unobs`: A vector of sampled values for the unobserved `W` (1 or 0).
#'   - `Post_prob`: A matrix of posterior probabilities for `W = 1` and `W = 0` for each sequence in `UW_loc`.
#' @export
#' 
w_samp_fun = function(R,G_samp,UW_loc,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp,pri_wp){
  W_samp_unobs = rep(NA, length(UW_loc))
  Post_prob = matrix(NA, nrow = total_n, ncol = 2)
  for (rep_i in UW_loc) {
    # rep_i=4
    Rui = R[rep_i]
    tmp_index = 1:str_length(Rui)
    Rui = str_sub(Rui,tmp_index,tmp_index)
    
    gui = G_samp[rep_i]
    aui = (A_samp[rep_i]:(A_samp[rep_i]+motif_len_w-1))
    bui = (B_samp[rep_i]:(B_samp[rep_i]+motif_len_g-1))
    llk_one_seq_w1 = exp(logllk_fun_one_seq(Rui,1,gui,aui,bui,dict,theta_0_samp,theta_samp,theta_til_samp))
    llk_one_seq_w0 = exp(logllk_fun_one_seq(Rui,0,gui,aui,bui,dict,theta_0_samp,theta_samp,theta_til_samp))
    
    post_wp = (pri_wp*llk_one_seq_w1)/(pri_wp*llk_one_seq_w1+(1-pri_wp)*llk_one_seq_w0)
    Post_prob[rep_i,] = c(post_wp,(1-post_wp))
    W_samp_unobs[which(UW_loc==rep_i)] =  sample(c(1,0),size = 1,replace=T,prob = c(post_wp,(1-post_wp)))
    # print(post_wp)
  }
  return(list(W_samp_unobs,Post_prob))
}