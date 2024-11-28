#' g_samp_fun
#'
#' @param  R A vector containing sequence data.
#' @param  W_samp A vector containing sampled W values.
#' @param  UG_loc A vector of indices where G is unobserved (NA values).
#' @param  A_samp A vector containing sampled A values.
#' @param  B_samp A vector containing sampled B values.
#' @param  motif_len_w The length of the first motif.
#' @param  motif_len_g The length of the second motif.
#' @param  dict A dictionary containing unique characters from the sequences.
#' @param  theta_0_samp Sampled values for theta_0.
#' @param  theta_samp Sampled values for theta.
#' @param  theta_til_samp Sampled values for theta_til.
#' @param  pri_gp The prior probability of G being 1.
#'
#' @return A list containing:
#'         - `G_samp_unobs`: The sampled values for unobserved G.
#'         - `Post_prob`: A matrix of posterior probabilities for each unobserved G.
#' @export
#'
g_samp_fun = function(R,W_samp,UG_loc,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp,pri_gp){
  G_samp_unobs = rep(NA, length(UG_loc))
  Post_prob = matrix(NA, nrow = total_n, ncol = 2)
  for (rep_i in UG_loc) {
    # rep_i=4
    Rui = R[rep_i]
    tmp_index = 1:str_length(Rui)
    Rui = str_sub(Rui,tmp_index,tmp_index)
    
    wui = W_samp[rep_i]
    aui = (A_samp[rep_i]:(A_samp[rep_i]+motif_len_w-1))
    bui = (B_samp[rep_i]:(B_samp[rep_i]+motif_len_g-1))
    llk_one_seq_g1 = exp(logllk_fun_one_seq(Rui,wui,1,aui,bui,dict,theta_0_samp,theta_samp,theta_til_samp))
    llk_one_seq_g0 = exp(logllk_fun_one_seq(Rui,wui,0,aui,bui,dict,theta_0_samp,theta_samp,theta_til_samp))
    
    post_gp = (pri_gp*llk_one_seq_g1)/(pri_gp*llk_one_seq_g1+(1-pri_gp)*llk_one_seq_g0)
    Post_prob[rep_i,] = c(post_gp,(1-post_gp))
    G_samp_unobs[which(UG_loc==rep_i)] =  sample(c(1,0),size = 1,replace=T,prob = c(post_gp,(1-post_gp)))
    # print(post_gp)
  }
  return(list(G_samp_unobs,Post_prob))
}