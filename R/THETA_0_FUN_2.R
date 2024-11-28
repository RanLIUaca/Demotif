#' h_00
#'
#' @param Ri A string representing a sequence of characters.
#' @param dict A dictionary of possible characters in the sequences.
#'
#' @return A vector representing the count of each character in `Ri`.
#' @export
#' 
h_00 = function(Ri,dict){
  Ri_factor = factor(Ri, levels = dict)
  return(table(Ri_factor))
}

#' h_01
#'
#' @param Ri A string representing a sequence of characters.
#' @param bi A vector of positions to exclude from `Ri`.
#' @param dict A dictionary of possible characters in the sequences.
#'
#' @return A vector representing the count of each character in `Ri` excluding positions in `bi`.
#' @export
#' 
h_01 = function(Ri,bi,dict){
  Ri_els = Ri[-bi]
  Ri_els_factor = factor(Ri_els, levels = dict)
  return(table(Ri_els_factor))
}

#' h_10
#'
#' @param Ri A string representing a sequence of characters.
#' @param ai A vector of positions to exclude from `Ri`.
#' @param dict A dictionary of possible characters in the sequences.
#'
#' @return A vector representing the count of each character in `Ri` excluding positions in `ai`.
#' @export
#' 
h_10 = function(Ri,ai,dict){
  Ri_els = Ri[-ai]
  Ri_els_factor = factor(Ri_els, levels = dict)
  return(table(Ri_els_factor))
}

#' h_11
#'
#' @param Ri A string representing a sequence of characters.
#' @param ai A vector of positions to exclude from `Ri`.
#' @param bi A vector of positions to exclude from `Ri`.
#' @param dict A dictionary of possible characters in the sequences.
#'
#' @return A vector representing the count of each character in `Ri` excluding positions in `ai` and `bi`.
#' @export
#' 
h_11 = function(Ri,ai,bi,dict){
  Ri_els = Ri[-union(ai,bi)]
  Ri_els_factor = factor(Ri_els, levels = dict)
  return(table(Ri_els_factor))
}


#' theta_0_samp_fun
#'
#' @param R A vector of string sequences.
#' @param W_samp A vector indicating the presence of the first motif (sampled).
#' @param G_samp A vector indicating the presence of the second motif (sampled).
#' @param A_samp A vector of positions for the first motif (sampled).
#' @param B_samp A vector of positions for the second motif (sampled).
#' @param motif_len_w Length of the first motif.
#' @param motif_len_g Length of the second motif.
#' @param dict A dictionary of possible characters in the sequences.
#' @param pri_al0 A prior Dirichlet parameter for the background probabilities.
#'
#' @return A vector representing the sampled background probabilities (`theta_0`).
#' @export
#' 
theta_0_samp_fun = function(R,W_samp,G_samp,A_samp,B_samp,motif_len_w,motif_len_g,dict,pri_al0){
  total_n=length(R)
  # post_parameter
  H_ab=rep(0,length(dict))
  names(H_ab)=dict
  for (rep_i in 1:total_n) {
    Ri = R[rep_i]
    tmp_index = 1:str_length(Ri)
    Ri = str_sub(Ri,tmp_index,tmp_index)
    
    ai = (A_samp[rep_i]:(A_samp[rep_i]+motif_len_w-1));
    bi = (B_samp[rep_i]:(B_samp[rep_i]+motif_len_g-1))
    if((W_samp[rep_i]==0)&(G_samp[rep_i]==0)){
      H_ab=H_ab+h_00(Ri,dict)
    }else if((W_samp[rep_i]==0)&(G_samp[rep_i]==1)){
      H_ab=H_ab+h_01(Ri,bi,dict)
    }else if((W_samp[rep_i]==1)&(G_samp[rep_i]==0)){
      H_ab=H_ab+h_10(Ri,ai,dict)
    }else{
      H_ab=H_ab+h_11(Ri,ai,bi,dict)
    }
  }
  
  post_alpha_0 = pri_al0+H_ab

  # sample parameter
  theta_0_samp = c(rdirichlet(1,post_alpha_0))
  
  return(theta_0_samp)
}