#' Data_process
#'
#' @param  Data A data frame containing the following columns:
#'         - Column 1: Sequence data (R), in string format.
#'         - Column 2: Observed W data (W_obs), numeric, may contain NA values.
#'         - Column 3: Observed G data (G_obs), numeric, may contain NA values.
#'
#' @return A list containing:
#'         - `Data`: The updated data frame with W_obs and G_obs columns added.
#'         - `total_n`: The total number of sequences.
#'         - `dict`: A dictionary of all unique characters in the sequences.
#'         - `len_dict`: The number of unique characters in the dictionary.
#'         - `Len_seq`: The length of each sequence.
#'         - `UW_loc`: The indices where W_obs is NA (unobserved).
#'         - `UG_loc`: The indices where G_obs is NA (unobserved).
#' @export
#'
Data_process = function(Data){
  total_n = length(Data[,1])
  R = Data[,1]; Data$W_obs = Data[,2]; Data$G_obs = Data[,3]
  Len_seq = rep(NA, total_n)
  for (rep_n in 1:total_n) {
    Ri = R[rep_n]
    Len_seq[rep_n] = str_length(Ri)
  }
  
  UG_loc = which(is.na(Data$G_obs))
  UW_loc = which(is.na(Data$W_obs))
  
  
  chars = unlist(strsplit(Data[,1], ""))
  dict = unique(chars)
  digit_dict = as.character(1:length(dict))
  names(digit_dict) = dict
  names(dict) = digit_dict
  len_dict = length(dict)
  return(list(Data,total_n,dict,len_dict,Len_seq,UW_loc,UG_loc))
}