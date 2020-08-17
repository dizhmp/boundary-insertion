insulation_across_insertions <- function (c1, c2, ins, nb, br, ws, file_name) {
  # c1, c2: insulation matrix for a control clone
  # ins: insulation matirx for a clone with insertion
  # nb: # of bins used for normalizing insulation score, must be smaller than matrix dimension
  # br: bin range, insulation score for # of bins on either side of the insertion
  # ws: max window size for calculating insulation score
  
  c1_col_no <- ncol(c1)
  c2_col_no <- ncol(c2)
  ins_col_no <- ncol(ins)
  
  c1_row_num_for_norm <- nrow(c1) - sum(is.na(c1[, ncol(c1)]))
  c2_row_num_for_norm <- nrow(c2) - sum(is.na(c1[, ncol(c2)]))
  ins_row_num_for_norm <- nrow(ins) - sum(is.na(c1[, ncol(ins)]))
  
  c1_row_num_from_center_for_norm <- (c1_row_num_for_norm - (c1_row_num_for_norm %% 2))/2 - 1
  c2_row_num_from_center_for_norm <- (c2_row_num_for_norm - (c2_row_num_for_norm %% 2))/2 - 1
  ins_row_num_from_center_for_norm <- (ins_row_num_for_norm - (ins_row_num_for_norm %% 2))/2 - 1
  
  # have to have same number of bins on either side of the insertion bin
  c1_for_norm <- c1[((nrow(c1)+1)/2 - c1_row_num_from_center_for_norm) : ((nrow(c1)+1)/2 + c1_row_num_from_center_for_norm) ,]
  c2_for_norm <- c2[((nrow(c2)+1)/2 - c2_row_num_from_center_for_norm) : ((nrow(c2)+1)/2 + c2_row_num_from_center_for_norm) ,]
  ins_for_norm <- ins[((nrow(ins)+1)/2 - ins_row_num_from_center_for_norm) : ((nrow(ins)+1)/2 + ins_row_num_from_center_for_norm) ,]
  
  c1_for_norm_colmean <- colMeans(c1_for_norm, na.rm=TRUE)
  c2_for_norm_colmean <- colMeans(c2_for_norm, na.rm=TRUE)
  ins_for_norm_colmean <- colMeans(ins_for_norm, na.rm=TRUE)
  
  # correct way of dividing each element in a column by the corresponding one in another vector
  log2_c1_normed <- log2(sweep(c1, MARGIN=2, FUN="/",STATS=c1_for_norm_colmean))
  log2_c2_normed <- log2(sweep(c2, MARGIN=2, FUN="/",STATS=c2_for_norm_colmean))
  log2_ins_normed <- log2(sweep(ins, MARGIN=2, FUN="/",STATS=ins_for_norm_colmean))
  
  x_ori <- ((nrow(c1)+1)/2 - br) : ((nrow(c1)+1)/2 + br)
  w_size_range <- 1 : ws
  
  log2_c1_normed_plot <- as.data.frame(log2_c1_normed[x_ori, w_size_range])
  log2_c2_normed_plot <- as.data.frame(log2_c2_normed[x_ori, w_size_range])
  log2_ins_normed_plot <- as.data.frame(log2_ins_normed[x_ori, w_size_range])
  
  library (ggplot2)
  library ('Rmisc')
  library (gridExtra)
  
  Bins_from_insertion <- -br: br
  
  # par(mfrow=c(p_row, p_col))
  # does not work with ggplot
  
pdf(file_name)
  for (i in 1:ws){
    plot_data_column(log2_c1_normed_plot, log2_c2_normed_plot, log2_ins_normed_plot, Bins_from_insertion, i)
  } 
dev.off()
}

  # myplots <- lapply(colnames(log2_c1_normed_plot), plot_data_column, data1 = log2_c1_normed_plot, data2= log2_c2_normed_plot , data3=log2_ins_normed_plot)
  
  # df <- rbind (WT, Clone25, Clone21)
  # p <- ggplot(df[1,], aes(Bins_from_insertion, y=value, color = samples)) +
  #   facet_wrap(~names, ncol=5) +
  #   geom_point(aes(y=df[1,], col='wt')) +
  #   geom_point(aes(y=df[2,], col='Clone 25')) +
  #   geom_point(aes(y=df[3,], col='Clone 21'))
  # p
  
  # melt command should be more elegant than using a for loop, but could not figure out the best way to do it!!
  # making a list of plots is super tricky, somehow multiplot only repeatedly plots the same plot in the iteration...
  
  # plots <- list()
  # for (i in w_size_range) {
  #   df <- data.frame(Bins_from_insertion, log2_c1_normed_plot[,i], log2_c2_normed_plot[,i], log2_ins_normed_plot[,i])
  #   p_i <- ggplot(df, aes(Bins_from_insertion, y=value, color = samples)) +
  #     geom_point(aes(y=log2_c1_normed_plot[,i], col='wt')) +
  #     geom_point(aes(y=log2_c2_normed_plot[,i], col='Clone 25')) +
  #     geom_point(aes(y=log2_ins_normed_plot[,i], col='Clone 21')) +
  #     ggtitle(paste("Insulation score for window size", i))
  #   plots[[i]] <- p_i
  # }
  # p