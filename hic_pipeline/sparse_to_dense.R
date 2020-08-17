# need to use read.table to load the matrix first
spar.to.dense <- function(x) {
  if (max(x$V1) - min(x$V1) != max(x$V2) - min(x$V2)){
    print("Error: Input sparse matrix is not square!")
  } else {
    dense_matrix <- matrix(, nrow = (max(x$V1) - min(x$V1) + 1), ncol = (max(x$V1) - min(x$V1) + 1))
    rownames(dense_matrix) <- c(seq(min(x$V1), max(x$V1)))
    colnames(dense_matrix) <- c(seq(min(x$V1), max(x$V1)))
    for (i in 1:nrow(x)){
          dense_matrix[toString(x[i,]$V1),toString(x[i,]$V2)] = x[i,]$V3
          dense_matrix[toString(x[i,]$V2),toString(x[i,]$V1)] = x[i,]$V3
      }
    dense_matrix[is.na(dense_matrix)] <- 0.1
  }
  dense_matrix
}