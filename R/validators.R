validate_integer <- function(alphabet) {
  if(!has_integers_only(alphabet)) {
    stop("alphabet should contain integers")
  }
}

validate_string <- function(alphabet) {
  if(!is.character(alphabet)) {
    stop("alphabet should contain strings")
  }
}

validate_numeric <- function(alphabet) {
  if(!is.numeric(alphabet) || has_integers_only(alphabet)) {
    stop("alphabet should contain numerics")
  }
}

is_empty <- function(elem) {
  (length(elem) == 1 && is.na(elem)) || is.null(elem)
}

has_integers_only <- function(v, tol=1e-9) {
  is.numeric(v) && all(abs(v - as.integer(v)) < tol)
}


is_positive_integer <- function(batch_size) {
  length(batch_size) == 1 && has_integers_only(batch_size) && batch_size > 0
}