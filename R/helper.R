#' Calculate decomposition error (Forbenius norm)
#'
#' @param matrix_initial list with initial Di matrices
#' @param delta_i delta values computed using SBF/A-SBF function
#' @param u_i  Ui values computed using SBF/A-SBF function
#' @param v_i V computed using SBF/A-SBF function
#'
#' @return
#' @export
#'
#' @examples
#' # load dataset
#' avg_counts <- SBF::avg_counts
#' # call sbf
#' sbf <- SBF(matrix_list = avg_counts, col_index = 2, approximate = FALSE,
#'            transform_matrix = FALSE, verbose = FALSE)
#'
#' # calculate decomposition error
#' decomperror <- calcDecompError(avg_counts, sbf$delta, sbf$u, sbf$v)
#'
#' # e.g. 2
#' asbf_cor <- SBF(matrix_list = avg_counts, col_index = 2, approximate = FALSE,
#'                transform_matrix = TRUE, verbose = FALSE)
#' decomperror <- calcDecompError(avg_counts, asbf_cor$delta, asbf_cor$u_ortho,
#'                                asbf_cor$v)
calcDecompError <- function(matrix_initial, delta_i, u_i, v_i) {
    decomp_error <- 0
    matrix_new <- list()
    for (mat in names(matrix_initial)) {
        if (length(delta_i[[mat]]) == 1)
            phi <- as.matrix(diag(as.matrix(delta_i[[mat]])))
        else phi <- as.matrix(diag(delta_i[[mat]]))
        matrix_new[[mat]] <- as.matrix(u_i[[mat]]) %*% phi %*%
                                                as.matrix(t(v_i))
        delta <- as.matrix(matrix_initial[[mat]]) - matrix_new[[mat]]
        decomp_error <- decomp_error + norm(delta, type = "F")^2
    }
    return(decomp_error)
}
