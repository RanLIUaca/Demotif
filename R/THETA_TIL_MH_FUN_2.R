#' theta_til_samp_fun
#'
#' @param R A vector of string sequences.
#' @param W_samp A vector indicating the presence of the first motif (sampled).
#' @param G_samp A vector indicating the presence of the second motif (sampled).
#' @param B_samp A vector of positions for the second motif (sampled).
#' @param motif_len_g Length of the second motif.
#' @param dict A dictionary of possible characters in the sequences.
#' @param pri_alg A prior Dirichlet parameter for the probabilities of the second motif.
#'
#' @return A list containing:
#'   - `theta_til_samp`: A matrix of sampled probabilities for each position in the second motif.
#'   - `Post_alpha_g`: A matrix of posterior Dirichlet parameters for each position in the second motif.
#' @export
#' 
theta_til_samp_fun = function(R,W_samp,G_samp,B_samp,motif_len_g,dict,pri_alg){
  total_n=length(R)
  
  # sample parameter
  theta_til_samp = matrix(nrow = length(dict), ncol = motif_len_g)
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
      
      bij = B_samp[rep_i]+rep_j-1
      
      
      if(((W_samp[rep_i]==1)&(G_samp[rep_i]==1))||((W_samp[rep_i]==0)&(G_samp[rep_i]==1))){
        Ri_bij = Ri[bij]
        tmp_index = 1:str_length(Ri_bij)
        Ri_bij = str_sub(Ri_bij,tmp_index,tmp_index)
        Ri_bij = factor(Ri_bij, levels = dict)
        H_b=H_b+table(Ri_bij)
      }
    }
    
    post_alpha_g = pri_alg+H_b
    Post_alpha_g[,rep_j] = post_alpha_g
    
    # print(post_alpha_g)
    temp = rdirichlet(1,post_alpha_g)
    theta_til_samp[,rep_j] = temp

  }
  # print(theta_til_samp)
  return(list(theta_til_samp,Post_alpha_g))
}