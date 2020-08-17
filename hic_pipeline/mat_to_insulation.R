# must be a dense square matrix
mat.to.insulation <- function (x, b) {
  # x the balanced full matrix
  # b is the upper limit of the window size
  insulation_score_matrix <- matrix(, nrow = nrow(x) - 2, ncol = b)
  rownames(insulation_score_matrix) <- c(seq((as.numeric(min(rownames(x))) + 1), (as.numeric(max(rownames(x))) - 1)))
  colnames(insulation_score_matrix) <- c(seq(1, b))
  
  for (i in 2:(nrow(x)-1)) {
    for (w in 1:b) {
      if (((i-w) < 0)||((i+w) > nrow(x)))  {
        insulation_score_matrix[i-1,w] <- NA
      } else {
          if (w==1) {
            insulation_score_matrix[i-1,w] <- x[(i+1):(i+w), (i-w):(i-1)]/w^2
          }
            else {
              windowed_mat <- (x[(i+1):(i+w), (i-w):(i-1)])
              insulation_score_matrix[i-1,w] <- sum(windowed_mat)/w^2
            }
      }
    }
  }
  insulation_score_matrix
}