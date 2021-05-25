#' Function to compute Shared Basis Factorization (SBF) and
#' Approximate Shared Basis Factorization (A-SBF)
#' @param avg_counts A list containing Di matrices for joint factorization.
#' Column name of each Di matrix will have information about tissue/celltype
#' separated by delimitter:species_tissue/celltype_otherdetails. For e.g.
#' "hsapiens_brain"  "hsapiens_heart"  "hsapiens_kidney" "hsapiens_liver"
#' "hsapiens_testis".
#' @param check_col_matching Check whether tissues/cell types are matching in
#' different Di. Default TRUE
#' @param colNameSep separator in column names to separate different fields.
#' For e.g "hsapiens_brain"  "hsapiens_heart" separator is "_". Default "_"
#' @param colIndex Index in column name corresponding to Tissue/Celltype.
#' For e.g. "hsapiens_brain" colIndex is 2. DEFAULT NULL
#' @param approximate TRUE will compute A-SBF. DEFAULT TRUE
#' @param transformData if TRUE then Di will be transfered to compute correlation
#' matrix and V is computed based on this. Otherwise V is computed using
#' Di^TDi. Default FALSE
#' @param verbose if TRUE print verbose lines. Default FALSE
#'
#' @return
#' @export
#'
#' @examples
#' data(avg_counts,package="SBF")
#' # SBF call
#' sbf = SBF(avg_counts = avg_counts,colIndex = 2, approximate = FALSE,
#'          transformData = FALSE)
#' # SBF call using correlation matrix
#' sbf.cor = SBF(avg_counts = avg_counts,colIndex = 2, approximate = FALSE,
#' transformData = TRUE)
#' # A-SBF call
#' asbf = SBF(avg_counts = avg_counts,colIndex = 2, approximate = TRUE,
#' transformData = FALSE)
#' # A-SBF call using correlation matrix
#' asbf.cor = SBF(avg_counts = avg_counts,colIndex = 2, approximate = TRUE,
#' transformData = TRUE)
SBF <- function(avg_counts = NULL, check_col_matching = TRUE,
                colNameSep="_", colIndex = NULL,
                approximate = TRUE, transformData = FALSE,verbose =FALSE){
  if(length(avg_counts) >= 2 & !is.null(avg_counts)){
    if(is.null(colIndex) | !is.numeric(colIndex))
      stop(paste0("\nInvalid Index to match Columns(Tissues/Celltypes).Exiting!"))
    if(check_col_matching){
      #require(data.table)
      c_tissues <- as.data.frame(sapply(avg_counts,
        function(x) data.table::tstrsplit(colnames(x),colNameSep)[[colIndex]]))
      if(!all(apply(c_tissues,1,function(x)all(x==x[1]))))
        stop(paste0("\nTissue names are not matching.Exiting!"))
    }
    if(!all(sapply(avg_counts,ncol) == ncol(avg_counts[[names(avg_counts)[1]]])))
      stop(paste0("\nAll data.frame should have same number of columns"))
    if(approximate)
      cat("\nA-SBF is computed")
    if(transformData)
      cat("\nV is computed using correlation matrix")
    species.names <- names(avg_counts)
    A <- list()
    S <- matrix(0L, nrow = ncol(avg_counts[[species.names[1]]]),
                ncol= ncol(avg_counts[[species.names[1]]]))
    if(transformData){
      for(species in species.names){
        A[[species]] <- stats::cor(avg_counts[[species]],method = "pearson")
        S <- S + A[[species]]
      }
    }else{
      for(species in species.names){
        A[[species]] <- t(as.matrix(avg_counts[[species]])) %*%
          as.matrix(avg_counts[[species]])
      }
    }
    S <- S/length(species.names)
    #eigen decomposition of S
    ev <- eigen(S)
    lambda <- ev$values
    V <- ev$vectors
    B = sigma = sigma.ortho = U = U.ortho = list()
    for(species in species.names){
      B[[species]] <- as.matrix(avg_counts[[species]]) %*% solve(t(V))
      sigma[[species]] <- sqrt(colSums(B[[species]]^2))
      U[[species]] <- sweep(B[[species]],2,sigma[[species]],"/")
      if(verbose){
        cat("\nspecies:",species)
        cat("\nDim of avg counts:",dim(avg_counts[[species]]))
        cat("\nDim of V:",dim(V))
        cat("\nDim of B:",dim(B[[species]]))
        cat("\nSigma vec:",sigma[[species]])
        if(length(sigma[[species]]) == 1)
          cat("\nDim of sigma:",
              dim(as.matrix(diag(as.matrix(sigma[[species]])))))
        else
          cat("\nDim of sigma:",
              dim(as.matrix(diag(sigma[[species]]))))
      }
      if(approximate){
        if(length(sigma[[species]]) == 1)
          phi <- as.matrix(diag(as.matrix(sigma[[species]])))
        else
          phi <- as.matrix(diag(sigma[[species]]))
        dvsig.svd <- svd(as.matrix(avg_counts[[species]]) %*%
                           V %*% phi)
        U.ortho[[species]] <- dvsig.svd$u %*% t(dvsig.svd$v)
        row.names(U.ortho[[species]]) <- row.names(U[[species]])
        sigma.ortho[[species]] <- diag(t(U.ortho[[species]]) %*%
                                         as.matrix(avg_counts[[species]]) %*%
                                         solve(t(V)))
      }
    }
    if(approximate){
      initial.error <- calcDecompError(avg_counts,sigma,U.ortho,V)
      initial.error.ortho <- calcDecompError(avg_counts,sigma.ortho,U.ortho,V)
      if(verbose){
        cat("\nInitial.Error =", round(initial.error$error,2),
            "after ortho sigma=", round(initial.error.ortho$error,2))
      }
    }
    if(approximate){
      output.list <- list(v = V, lambda = lambda, u = U,
                          u_ortho = U.ortho,
                          sigma = sigma,
                          sigma.ortho = sigma.ortho,
                          s = S,
                          error = initial.error$error,
                          error.ortho = initial.error.ortho$error)
    }else{
      output.list <- list(v = V, lambda = lambda, u = U,
                          sigma = sigma,
                          s = S)
    }
    return(output.list)
  }else{
    stop("\nInvalid counts file. counts file should be list \
              with matching column names.Exiting!")
  }
}
