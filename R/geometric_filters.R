
#' Estimate the probability density at each point in a sample.
#' @details Use KDE with a Gaussian kernal.
#' @param distances The distances among the sample.
#' @param bandwidth A characteristic distance between neighbors.
#' @value A vector with element's local density estimates.
#' @export

estimate_local_density <- function(distances, bandwidth){
  kernal          <- function(x) dnorm(x, sd = bandwidth)
  distance_matrix <- as.matrix(distances)
  image_matrix    <- kernal(distance_matrix)
  # don't count the point's contribution  to its own density
  diag(image_matrix) <- 0
  image_sums      <- apply(image_matrix, 1, sum)
  n_minus_one     <- attr(distances, 'Size') - 1
  image_sums / n_minus_one
}


#' Calculate the p-eccentricity for each element in a sample.
#' @details $$ e_p(x) \equiv (\mean_{y \neq x} dist(x, y) ^ p) ^ {1 \over p} $$
#' @param distances The distances among the sample.
#' @param p The exponenent to raise the distances to before taking the mean.
#' @value A vector of elements' eccentricities.
#' @export

estimate_local_eccentricity <- function(distances, p){
  distance_matrix <- distances %>% as.matrix()
  image_matrix    <- distance_matrix ** p
  # point's self-distances are zero
  image_sums      <- apply(image_matrix, 1, sum)
  n_minus_one     <- attr(distances, 'Size') - 1
  (image_sums / n_minus_one) ** (1 / p)
}


#' Find the distance to its n-th neighbor for each element in a sample.
#' @param distances The distances among the sample.
#' @param n Which neighbor, ordered by distance ascending.
#' @value A vector of elements' neighbor distances.
#' @export

find_neighbor_distance <- function(distances, n){
  distances %>% 
    as.matrix() %>%
    apply(1, function(x){
      sort(x, partial = n + 1)[n + 1]
    })
}


#' Calculate the mean distance of its first n neighbors for each element in
#' a sample.
#' @param distances The distances among the sample.
#' @param n How many neighbors, ordered by distance ascending.
#' @value A vector of elements' mean neighbor distances.
#' @export

estimate_neighbor_distance <- function(distances, n){
  distances %>% 
    as.matrix() %>%
    apply(1, function(x){
      sort(x, partial = n + 1)[2: n + 1] %>% mean()
    })
}

RBF_kernal <- function(bandwidth) {
  function(d) exp(-0.5 * (d / bandwidth) ** 2)
}

#' Construct the normalized Laplacian matrix for a sample.
#' @param distances The distances among the sample.
#' @param bandwidth A characteristic distance between neighbors.
#' @param kernal A kernal function for specifying weighted adjacency. Defaults
#' to an RBF kernal of the specified `bandwidth`.
#' @details If given distances and bandwidth, applies the RBF kernel to
#' calculate the normalized laplacian matrix.
#' @value The Laplacian matrix.
#' @export

normalized_laplacian_matrix <- function(
  distances, bandwidth, kernal = RBF_kernal(bandwidth)
){
  n                <- attr(distances, "Size")
  distance_matrix  <- as.matrix(distances)
  similarity       <- kernal(distance_matrix)
  total_similarity <- apply(similarity, 1, sum)
  rows             <- replicate(n, total_similarity) 
  cols             <- rows %>% t()
  diag(1, n) - similarity / sqrt(rows * cols)
}

# Find the index of the first non-approximately-zero entry in a sorted vector.
which_first_nonzero <- function(sorted_xs, descending = FALSE){
  if (descending) {
    xs <- rev(sorted_xs)
  } else {
    xs <- sorted_xs
  }
  # the absolute smallest nonzero entry, or 1e-12
  epsilon <- min(1e-12, min(abs(xs[xs != 0])))
  # peg the actual minimum to epsilon
  xs      <- xs - min(xs) + epsilon
  # take logs for geometric comparison
  # impute some new mimimum for -Inf
  # subtract (divide) and find minimum
  xs      <- log(xs)
  xs_lead <- c(xs[2:length(xs)], -Inf)
  xs_diff <- xs_lead - xs 
  retval  <- 1 + which.max(xs_diff)
  if (descending) {
    return(length(sorted_xs) + 1 - retval)
  } else {
    return(retval)
  }
}

#' Find the Fiedler vector, the eigenvector corresponding to the first nonzero
#' eigenvalue of the normalized Laplacian matrix.
#' @param distances The distances among the sample.
#' @param bandwidth A characteristic distance between neighbors.
#' @param n How many eigenvectors to return, defaults to 1.
#' @param eigen_solution Optional, pass the eigen solution for the
#' Laplacian matrix.
#' @details If given distances and bandwidth, applies the RBF kernel to
#' calculate the normalized Laplacian matrix, then solves the eigensystem.
#' @value The Fiedler vector(s).
#' @export

estimate_fiedler_vector <- function(distances, bandwidth,
                                    n = 1, eigen_solution = NULL){
  if (is.null(eigen_solution)) {
    laplacian <- normalized_laplacian_matrix(distances, bandwidth)
    eigen_solution <- eigen(laplacian)
  }
  vectors <- eigen_solution$vectors
  values  <- eigen_solution$values
  n       <- length(values)
  i       <- which_first_nonzero(values, descending = TRUE)
  vectors[, (i - n + 1):i]
}

