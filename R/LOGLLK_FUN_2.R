#' logllk_fun_00
#'
#' This function calculates the log-likelihood of a sequence when both motifs are absent (wi = 0, gi = 0).
#'
#' @param  Ri A string representing the sequence.
#' @param  theta_0_samp A vector of sampled probabilities for the background model.
#' @param  dict A dictionary of unique characters in the sequence.
#'
#' @return A numeric value representing the log-likelihood of the sequence.
#' @export
#' 
logllk_fun_00 = function(Ri,theta_0_samp,dict){
  Ri_factor = factor(Ri, levels = dict)
  res=sum((table(Ri_factor))*log(theta_0_samp))
  return(res)
}

#' logllk_fun_01
#'
#' This function calculates the log-likelihood of a sequence when the second motif is present (gi = 1) and the first motif is absent (wi = 0).
#'
#' @param  Ri A string representing the sequence.
#' @param  bi The positions of the second motif within the sequence.
#' @param  theta_0_samp A vector of sampled probabilities for the background model.
#' @param  theta_til_samp A matrix of sampled probabilities for the second motif.
#' @param  dict A dictionary of unique characters in the sequence.
#'
#' @return A numeric value representing the log-likelihood of the sequence.
#' @export
#' 
logllk_fun_01 = function(Ri,bi,theta_0_samp,theta_til_samp,dict){
  Ri_g = Ri[bi]; Ri_els = Ri[-bi]
  Ri_els_factor = factor(Ri_els, levels = dict)
  res=sum((table(Ri_els_factor))*log(theta_0_samp))
  for (rep_j in 1:length(bi)) {
    Ri_bij =  Ri_g[rep_j]
    Ri_bij_factor = factor(Ri_bij, levels = dict)
    res=res+sum((table(Ri_bij_factor))*log(theta_til_samp[,rep_j]))
  }
  return(res)
}

#' logllk_fun_10
#'
#' This function calculates the log-likelihood of a sequence when the first motif is present (wi = 1) and the second motif is absent (gi = 0).
#'
#' @param  Ri A string representing the sequence.
#' @param  ai The positions of the first motif within the sequence.
#' @param  theta_0_samp A vector of sampled probabilities for the background model.
#' @param  theta_samp A matrix of sampled probabilities for the first motif.
#' @param  dict A dictionary of unique characters in the sequence.
#'
#' @return A numeric value representing the log-likelihood of the sequence.
#' @export
#' 
logllk_fun_10 = function(Ri,ai,theta_0_samp,theta_samp,dict){
  Ri_w = Ri[ai]; Ri_els = Ri[-ai]
  Ri_els_factor = factor(Ri_els, levels = dict)
  res=sum((table(Ri_els_factor))*log(theta_0_samp))
  for (rep_j in 1:length(ai)) {
    Ri_aij =  Ri_w[rep_j]
    Ri_aij_factor = factor(Ri_aij, levels = dict)
    res=res+sum((table(Ri_aij_factor))*log(theta_samp[,rep_j]))
  }
  return(res)
}

#' logllk_fun_11
#'
#' This function calculates the log-likelihood of a sequence when both motifs are present (wi = 1, gi = 1).
#'
#' @param  Ri A string representing the sequence.
#' @param  ai The positions of the first motif within the sequence.
#' @param  bi The positions of the second motif within the sequence.
#' @param  theta_0_samp A vector of sampled probabilities for the background model.
#' @param  theta_samp A matrix of sampled probabilities for the first motif.
#' @param  theta_til_samp A matrix of sampled probabilities for the second motif.
#' @param  dict A dictionary of unique characters in the sequence.
#'
#' @return A numeric value representing the log-likelihood of the sequence.
#' @export
#' 
logllk_fun_11 = function(Ri,ai,bi,theta_0_samp,theta_samp,theta_til_samp,dict){
  Ri_g = Ri[bi]; Ri_w = Ri[setdiff(ai, bi)]
  Ri_els = Ri[-union(ai,bi)]
  Ri_els_factor = factor(Ri_els, levels = dict)
  res=sum((table(Ri_els_factor))*log(theta_0_samp))
  for (rep_j in 1:length(bi)) {
    #rep_j=1
    Ri_bij =  Ri_g[rep_j]
    Ri_bij_factor = factor(Ri_bij, levels = dict)
    res=res+sum((table(Ri_bij_factor))*log(theta_til_samp[,rep_j]))
  }
  if(length(setdiff(ai, bi))>0){
  for (rep_j in 1:length(setdiff(ai, bi))) {
    Ri_aij =  Ri_w[rep_j]
    Ri_aij_factor = factor(Ri_aij, levels = dict)
    res=res+sum((table(Ri_aij_factor))*log(theta_samp[,rep_j]))
  }
  }
  return(res)
}


#' logllk_fun_one_seq
#'
#' This function calculates the log-likelihood of a single sequence based on the presence or absence of motifs.
#'
#' @param  Ri A string representing the sequence.
#' @param  wi Indicator for the first motif (1 if present, 0 otherwise).
#' @param  gi Indicator for the second motif (1 if present, 0 otherwise).
#' @param  ai The positions of the first motif within the sequence.
#' @param  bi The positions of the second motif within the sequence.
#' @param  dict A dictionary of unique characters in the sequence.
#' @param  theta_0_samp A vector of sampled probabilities for the background model.
#' @param  theta_samp A matrix of sampled probabilities for the first motif.
#' @param  theta_til_samp A matrix of sampled probabilities for the second motif.
#'
#' @return A numeric value representing the log-likelihood of the sequence.
#' @export
#' 
logllk_fun_one_seq = function(Ri,wi,gi,ai,bi,dict,theta_0_samp,theta_samp,theta_til_samp){
  if((wi==0)&(gi==0)){
    res=logllk_fun_00(Ri,theta_0_samp,dict)
  }else if((wi==0)&(gi==1)){
    res=logllk_fun_01(Ri,bi,theta_0_samp,theta_til_samp,dict)
  }else if((wi==1)&(gi==0)){
    res=logllk_fun_10(Ri,ai,theta_0_samp,theta_samp,dict)
  }else{
    res=logllk_fun_11(Ri,ai,bi,theta_0_samp,theta_samp,theta_til_samp,dict)
  }
  # print(c(rep_i,res))
  
  return(res)
}

#' logllk_fun
#'
#' This function calculates the total log-likelihood for a set of sequences based on the sampled motifs.
#'
#' @param  R A vector of sequences.
#' @param  W_samp A vector indicating the presence/absence of the first motif (W).
#' @param  G_samp A vector indicating the presence/absence of the second motif (G).
#' @param  A_samp A vector of positions for the first motif in each sequence.
#' @param  B_samp A vector of positions for the second motif in each sequence.
#' @param  motif_len_w The length of the first motif.
#' @param  motif_len_g The length of the second motif.
#' @param  dict A dictionary of unique characters in the sequences.
#' @param  theta_0_samp A vector of sampled probabilities for the background model.
#' @param  theta_samp A matrix of sampled probabilities for the first motif.
#' @param  theta_til_samp A matrix of sampled probabilities for the second motif.
#'
#' @return A numeric value representing the total log-likelihood across all sequences.
#' @export
#' 
logllk_fun = function(R,W_samp,G_samp,A_samp,B_samp,motif_len_w,motif_len_g,dict,theta_0_samp,theta_samp,theta_til_samp){
  res=0
  total_n=length(R)
  for (rep_i in 1:total_n) {
    #rep_i=17
    Ri = R[rep_i]
    tmp_index = 1:str_length(Ri)
    Ri = str_sub(Ri,tmp_index,tmp_index)
    
    ai = (A_samp[rep_i]:(A_samp[rep_i]+motif_len_w-1));
    bi = (B_samp[rep_i]:(B_samp[rep_i]+motif_len_g-1))
    if((W_samp[rep_i]==0)&(G_samp[rep_i]==0)){
      res=res+logllk_fun_00(Ri,theta_0_samp,dict)
    }else if((W_samp[rep_i]==0)&(G_samp[rep_i]==1)){
      res=res+logllk_fun_01(Ri,bi,theta_0_samp,theta_til_samp,dict)
    }else if((W_samp[rep_i]==1)&(G_samp[rep_i]==0)){
      res=res+logllk_fun_10(Ri,ai,theta_0_samp,theta_samp,dict)
    }else{
      res=res+logllk_fun_11(Ri,ai,bi,theta_0_samp,theta_samp,theta_til_samp,dict)
    }
    # print(c(rep_i,res))
  }
  return(res)
}