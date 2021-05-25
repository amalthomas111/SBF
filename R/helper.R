#' Title Function to calulate Frobenium norm for decomposition
#'
#' @param matrix.list.start list with Di matrix
#' @param sig_i sigma values computed using SBF/A-SBF function
#' @param u_i  Ui values computed using SBF/A-SBF function
#' @param v_i V computed using SBF/A-SBF function
#'
#' @return
#' @export
#'
#' @examples
#' data(avg_counts,package="SBF")
#' sbf = SBF(avg_counts = avg_counts,colIndex = 2, approximate = FALSE,
#' transformData = FALSE,verbose = FALSE)
#' calcDecompError(avg_counts,sbf$sigma,sbf$u,sbf$v)
#' sbf.cor = SBF(avg_counts = avg_counts,colIndex = 2, approximate = FALSE,
#' transformData = TRUE,verbose = FALSE)
#' calcDecompError(avg_counts,sbf.cor$sigma,sbf.cor$u,sbf.cor$v)
calcDecompError <- function(matrix.list.start,sig_i,u_i,v_i){
  decomp.error <- 0
  matrix.list.new <- list()
  for(species in names(matrix.list.start)){
    if(length(sig_i[[species]]) == 1)
      phi <- as.matrix(diag(as.matrix(sig_i[[species]])))
    else
      phi <- as.matrix(diag(sig_i[[species]]))
    matrix.list.new[[species]] <- as.matrix(u_i[[species]]) %*% phi %*%
      as.matrix(t(v_i))
    delta <- as.matrix(matrix.list.start[[species]]) - matrix.list.new[[species]]
    decomp.error <-  decomp.error + norm(delta,type="F")^2
  }
  return(list(error=decomp.error))
}
