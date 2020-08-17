plot_data_column <- function (data1, data2, data3, bins, column) {
  p <- ggplot(data1, aes(bins, y=value, color = samples)) +
    geom_line(aes(y=data1[, column], col='WT'), size=3) + 
    geom_line(aes(y=data2[, column], col='Clone 25'), size=3) +
    geom_line(aes(y=data3[, column], col='Clone 21'), size=3) +
    ggtitle(paste("Insulation score for window size", column)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_blank(), panel.border = element_rect(colour = "grey35", fill= NA, size=1)) +
    scale_color_manual(values=c("darkorange", "deepskyblue", "grey31")) + 
    theme(legend.key=element_blank(), legend.key.width = unit(2, "cm"))
  print(p)
}