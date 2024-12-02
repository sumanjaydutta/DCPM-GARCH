moore_penrose_pseudoinverse <- function(A, tol = .Machine$double.eps) {
  # Perform SVD on matrix A
  svd_decomp <- svd(A)
  
  # Extract components
  U <- svd_decomp$u
  D <- svd_decomp$d
  V <- svd_decomp$v
  
  # Invert the non-zero singular values
  D_plus <- diag(ifelse(D > tol, 1 / D, 0))
  
  # Calculate the pseudo-inverse
  A_plus <- V %*% D_plus %*% t(U)
  
  return(A_plus)
}

# Example usage
A <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
A_pseudo_inverse <- moore_penrose_pseudoinverse(A)
print(A_pseudo_inverse)
