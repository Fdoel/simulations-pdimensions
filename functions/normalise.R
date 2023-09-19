normalise <- function(M) {
 M * det(M)^-(1/(nrow(M)))
}

