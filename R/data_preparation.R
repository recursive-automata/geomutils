
#' Calculate the mode of the elements of a vector.
#' @param x A numeric vector.
#' @return A numeric scalar.
#' @export

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Clip the values in a sample at upper and lower quantiles of the sample.
#' @param x A numeric vector.
#' @param p_min A numeric scalar between 0 and 1, the lower quantile.
#' @param p_max A numeric scalar between 0 and 1, the upper quantile.
#' @value A numeric vector of clipped values.
#' @export

clip_pct <- function(x, p_min, p_max){
  bounds <- quantile(x, probs = c(p_min, p_max))
  a      <- bounds[1]
  b      <- bounds[2]
  clip(x, a, b)
}


#' Find the cumulative percentile for each element in the sample.
#' @param x A numeric vector.
#' @details `value == rank(x, ties.method = "average") / length(x)`
#' @value A numeric vector of values between 0 and 1.
#' @export

percentilize <- function(x){
  r <- rank(x, ties.method = "average")
  r / length(x)
}


#' Make one-hot-encoded columns for a vector of characteristics.
#' @param chars A vector of characteristics, represented as `character`.
#' @param char_values Optional, the set of unique characteristics.
#' @param prefix How to prefix coluns in the output matrix.
#' @value A matrix of one-hot-encoded columns.
#' @export

make_ohe_columns <- function(chars, char_values = NULL, prefix = ''){
  if (typeof(char_values) != 'character') {
    char_values <- unique(chars)
  }
  ohe_columns <- char_values %>%
    lapply(function(cat) {as.numeric(chars == cat)}) %>%
    do.call(cbind, .)
  colnames(ohe_columns) <- paste0(prefix, char_values)
  ohe_columns
}
