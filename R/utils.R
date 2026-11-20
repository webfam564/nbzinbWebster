# Small internal helpers shared by nbreg() and zinbreg().

# Stop unless x is a single positive finite number.
check_pos <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("'%s' must be a single positive number.", name),
         call. = FALSE)
  }
}

# Stop unless x is TRUE or FALSE.
check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
  }
}

# Validate the response, prior weights and offset shared by both fitters.
check_count_response <- function(y) {
  if (is.null(y)) stop("The formula must have a response.", call. = FALSE)
  if (is.matrix(y)) stop("The response must be a vector.", call. = FALSE)
  if (any(!is.finite(y)) || any(y < 0) || any(abs(y - round(y)) > 1e-8)) {
    stop("The response must contain non-negative integer counts.",
         call. = FALSE)
  }
  invisible(as.vector(y))
}

check_weights <- function(a, n) {
  if (is.null(a)) return(rep(1, n))
  if (!is.numeric(a) || any(!is.finite(a)) || any(a < 0)) {
    stop("'weights' must be non-negative and finite.", call. = FALSE)
  }
  as.vector(a)
}

check_offset <- function(o, n) {
  if (is.null(o)) return(rep(0, n))
  o <- as.vector(o)
  if (length(o) != n || any(!is.finite(o))) {
    stop("'offset' must be finite and have one value per observation.",
         call. = FALSE)
  }
  o
}
