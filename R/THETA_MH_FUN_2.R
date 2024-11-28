#' theta_samp_fun
#'
#' @param R A vector of string sequences.
#' @param W_samp A vector indicating the presence of the first motif (sampled).
#' @param G_samp A vector indicating the presence of the second motif (sampled).
#' @param A_samp A vector of positions for the first motif (sampled).
#' @param B_samp A vector of positions for the second motif (sampled).
#' @param motif_len_w Length of the first motif.
#' @param motif_len_g Length of the second motif.
#' @param dict A dictionary of possible characters in the sequences.
#' @param pri_alw A prior Dirichlet parameter for the probabilities of the first motif.
#'
#' @return A list containing:
#'   - `theta_samp`: A matrix of sampled probabilities for each position in the first motif.
#'   - `Post_alpha_w`: A matrix of posterior Dirichlet parameters for each position in the first motif.
#' @export
#' 
theta_samp_fun = function(R,W_samp,G_samp,A_samp,B_samp,motif_len_w,motif_len_g,dict,pri_alw){
  total_n=length(R)
  
  # sample parameter
  theta_samp = matrix(nrow = length(dict), ncol = motif_len_w)
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
      
      aij = A_samp[rep_i]+rep_j-1
      bi = (B_samp[rep_i]:(B_samp[rep_i]+motif_len_g-1))
      # A_samp[rep_i,]

      if(((W_samp[rep_i]==1)&(G_samp[rep_i]==1))){
        aij = setdiff(aij,bi)
      }
      
      if(((W_samp[rep_i]==1)&(G_samp[rep_i]==0))||((W_samp[rep_i]==1)&(G_samp[rep_i]==1))){
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
    
    temp = rdirichlet(1,post_alpha_w)
    theta_samp[,rep_j] = temp
  }
  
  return(list(theta_samp,Post_alpha_w))
}