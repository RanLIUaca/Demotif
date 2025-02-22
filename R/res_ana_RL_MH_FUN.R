#' res_ana
#'
#' This function performs the analysis of motif estimation and sampling results. It includes reading the sampled data, estimating parameters, plotting results, and analyzing motif jumps.
#'
#' @param motif_len_w The length of the first binding motif.
#' @param motif_len_g The length of the second binding motif.
#' @param N The number of sampling iterations.
#' @import ggplot2
#' @import ggseqlogo
#' @import ggpubr
#' @return A list containing:
#'   - Samples for `theta_0` (background probabilities),
#'   - Samples for `theta` (the first binding motif probabilities),
#'   - Samples for `tilde{theta}` (the second binding motif probabilities),
#'   - Samples for `W` and `G` (indicators for motif presence),
#'   - Samples for `A` and `B` (positions of motifs),
#'   - Log-likelihood values,
#'   - Estimated parameters (`theta_0`, `theta`, `tilde{theta}`),
#'   - The number of jumps for both motifs,
#'   - a plot for the curve of the loglikelihood, the estimated theta logo and the estimated tilde{theta} logo,
#'   - a plot for the logos before and after a jump.
#'    @export
#' 
res_ana = function(motif_len_w,motif_len_g,burnin, N){
  # UG_loc = which(is.na(G_obs))
  # UW_loc = which(is.na(W_obs))
  
  ## read all samples ###
  dir_name = './result'
  if(dir.exists(dir_name)==0){dir.create(dir_name)}
  dir_name = paste0(dir_name,'/temp_data')
  if(dir.exists(dir_name)==0){dir.create(dir_name)}
  model_name = paste0(dir_name,'/motif_len_w=',motif_len_w,'_motif_len_g=',motif_len_g)
  if(dir.exists(dir_name)==0){dir.create(dir_name)}
  
  Theta_0 = matrix(NA, nrow=len_dict, ncol = N)
  Theta = rep(list(matrix(NA,nrow=len_dict,ncol = motif_len_w)),N)
  Theta_til = rep(list(matrix(NA,nrow=len_dict,ncol = motif_len_g)),N)
  Logllk = rep(NA,N)
  W_samp = matrix(NA, nrow=total_n, ncol = N)
  G_samp = matrix(NA, nrow=total_n, ncol = N)
  A_samp = rep(list(matrix(NA,nrow=total_n,ncol = motif_len_w)),N)
  B_samp = rep(list(matrix(NA,nrow=total_n,ncol = motif_len_g)),N)

  for (rep_N in 1:N){
    #rep_N=1
    # print(rep_N)
    Theta_0[,rep_N] = read.csv(paste0(model_name ,'/theta_0_rep_N=',rep_N,'.csv'),header = T)[,-1]
    Theta[[rep_N]] = as.matrix(read.csv(paste0(model_name,'/theta_rep_N=',rep_N,'.csv'),header = T)[,-1])
    Theta_til[[rep_N]] = as.matrix(read.csv(paste0(model_name,'/theta_til_rep_N=',rep_N,'.csv'),header = T)[,-1])
    Logllk[rep_N] = as.numeric(read.csv(paste0(model_name,'/logllk_rep_N=',rep_N,'.csv'),header = T)[,-1])
    
    if(length(UW_loc)>0){
    W_samp[,rep_N] = as.numeric(read.csv(paste0(model_name,'/W_rep_N=',rep_N,'.csv'),header = T)[,-1])
    }
    if(length(UG_loc)>0){
    G_samp[,rep_N] = as.numeric(read.csv(paste0(model_name,'/G_rep_N=',rep_N,'.csv'),header = T)[,-1])
    }
    A_samp[[rep_N]] = as.matrix(read.csv(paste0(model_name,'/A_rep_N=',rep_N,'.csv'),header = T)[,-1])
    B_samp[[rep_N]] = as.matrix(read.csv(paste0(model_name,'/B_rep_N=',rep_N,'.csv'),header = T)[,-1])
  }
  
  ### read MH jump theta theta_til ###
  num_change = read.csv(paste0(model_name,'/num_change.csv'),header = T)[,-1]
  if(num_change[1,1]!=0) {
    loc_jump_theta = as.numeric(read.csv(paste0(model_name,'/loc_jump_w.csv'),header = T)[,-1])[-1]
    theta_pre = rep(list(),num_change[1,1]);theta_post = rep(list(),num_change[1,1])
    for (rep_N in 1:num_change[1,1]) {
      theta_pre[[rep_N]] = as.matrix(read.csv(paste0(model_name,'/change_theta_pre_rep_N=',loc_jump_theta[rep_N],'.csv'),header = T)[,-1])
      theta_post[[rep_N]] = as.matrix(read.csv(paste0(model_name,'/change_theta_post_rep_N=',loc_jump_theta[rep_N],'.csv'),header = T)[,-1])
      row.names(theta_pre[[rep_N]]) = dict
      row.names(theta_post[[rep_N]]) = dict
    }
  }
  if(num_change[1,2]!=0) {
    loc_jump_theta_til = as.numeric(read.csv(paste0(model_name,'/loc_jump_g.csv'),header = T)[,-1])[-1]
    theta_til_pre = rep(list(),num_change[1,2]);theta_til_post = rep(list(),num_change[1,2])
    for (rep_N in 1:num_change[1,2]) {
      theta_til_pre[[rep_N]] = as.matrix(read.csv(paste0(model_name,'/change_theta_til_pre_rep_N=',loc_jump_theta_til[rep_N],'.csv'),header = T)[,-1])
      theta_til_post[[rep_N]] = as.matrix(read.csv(paste0(model_name,'/change_theta_til_post_rep_N=',loc_jump_theta_til[rep_N],'.csv'),header = T)[,-1])
      row.names(theta_til_pre[[rep_N]]) = dict
      row.names(theta_til_post[[rep_N]]) = dict
    }
  }
  
  #####################
  ### plot #############
  #######################
  plot_list = list()
  ##############################
  ############# logllk curve ###########
  ###############################
  data_logllk = data.frame(x = 1:length(Logllk),value = Logllk)
  
  logllk_plot = ggplot(data_logllk, aes(x = x, y = value)) +
    geom_line(color = "#3498DB") +
    labs(title = "Loglikelihood", x = " ", y = " ") +
    theme_bw()+theme_classic()+
    theme(panel.grid=element_blank(),plot.title = element_text(size = 10),
          axis.text=element_text(size=10),
          legend.key = element_blank(), 
          legend.title = element_blank(),
          legend.text = element_text(size = 10),  
          legend.position = "none")
  plot_list [[1]] = logllk_plot
  
  #########################
  ### parameters logo ####
  ##########################
  max_loc = burnin  + which.max(Logllk[(burnin+1):N])
  est_theta_0 = Theta_0[,max_loc]
  est_theta = Theta[[max_loc]]; est_theta = as.matrix(est_theta); row.names(est_theta) = dict
  plot_list[[2]] = ggseqlogo(est_theta) +
    theme(legend.position = "none", axis.text = element_text(size = 11),plot.title = element_text(size = 11))+
    labs(title = "First Motif Estimation") 

  est_theta_til = Theta_til[[max_loc]]; est_theta_til = as.matrix(est_theta_til); row.names(est_theta_til) = dict
  plot_list[[3]] = ggseqlogo(est_theta_til)+
    theme(legend.position = "none", axis.text = element_text(size = 11),plot.title = element_text(size = 11))+
    labs(title = "Second Motif Estimation") 

  
  all.plot=ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],
                     nrow=1,  widths = c(1, 1, 1),  
                     common.legend = FALSE)
  
  All_plot = annotate_figure(all.plot, left = text_grob(paste0('Lengths of the first and second motif are ',motif_len_w,' and ', motif_len_g), rot = 90, vjust = 1, size = 9))
  ###################################
  ###### Jump logo difference ####
  ############################
  jump_plot=list()
  num_change = read.csv(paste0(model_name ,'/num_change.csv'),header = T)[,-1]
  if(num_change[1,1]!=0) {
    change_logo = list()
    for (rep_N in 1:num_change[1,1]) {
      change_logo[[1]] = ggseqlogo(theta_pre[[rep_N]]) +theme(legend.position = "none")
      change_logo[[2]] = ggseqlogo(theta_post[[rep_N]]) +theme(legend.position = "none")
      change.plot=ggarrange(change_logo[[1]],change_logo[[2]],nrow=1, common.legend = FALSE)
      Change_plot = annotate_figure(change.plot,top = text_grob(paste0('Jump of the First Binding Motif When Time = ', loc_jump_theta[rep_N]), size = 11,hjust = 1)) 
      jump_plot[[rep_N]] = Change_plot
    }
  }
  if(num_change[1,2]!=0) {
    change_logo = list()
    for (rep_N in 1:num_change[1,2]) {
      change_logo[[1]] = ggseqlogo(theta_til_pre[[rep_N]]) +theme(legend.position = "none")
      change_logo[[2]] = ggseqlogo(theta_til_post[[rep_N]]) +theme(legend.position = "none")
      change.plot=ggarrange(change_logo[[1]],change_logo[[2]],nrow=1, common.legend = FALSE)
      Change_plot = annotate_figure(change.plot,top = text_grob(paste0('Jump of the Second Binding Motif When Time = ', loc_jump_theta[rep_N]), size = 11,hjust = 1)) 
      jump_plot[[num_change[1,1]+rep_N]] = Change_plot
    }
  }
  JUMP_plot = NA
  if(sum(num_change[1,])!=0) {
    jump.plot = do.call(ggarrange, c(jump_plot, ncol = 1, common.legend = FALSE))
    JUMP_plot = annotate_figure(jump.plot,top = text_grob(paste0('Case is (',motif_len_w,',',motif_len_g,')' ),face='bold',size = 11,hjust = 0.5))
  }
return(list(Theta_0,Theta,Theta_til,W_samp,G_samp,A_samp,B_samp,Logllk,est_theta_0, est_theta, est_theta_til,num_change,All_plot,JUMP_plot))
}
