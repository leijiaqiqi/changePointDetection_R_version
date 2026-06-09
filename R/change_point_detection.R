#' Fixed-Segment Change-Point Detection with Gurobi
#'
#' Solves a mixed-integer optimization formulation for segmented empirical
#' cumulative distribution functions. The number of segments is fixed by the
#' user. The returned vector gives the cumulative ending positions of the
#' detected segments.
#'
#' @param Length Integer. Number of segments. Must be at least 2.
#' @param X Numeric vector containing the sequential observations.
#' @param Delta Integer. Minimum number of observations assigned to each
#'   segment. Must be positive.
#' @param output_flag Integer. Use 1 to display Gurobi output and 0 to suppress
#'   it. Defaults to 1.
#'
#' @return A numeric vector of cumulative segment ending positions.
#'
#' @examples
#' \dontrun{
#' x <- c(1, 1, 2, 2, 10, 10)
#' change_point_detection_fixed_num(Length = 2, X = x, Delta = 2)
#' }
#'
#' @export
change_point_detection_fixed_num <- function(Length, X, Delta, output_flag = 1L) {
  L <- as.integer(Length)
  X <- as.numeric(X)
  delta <- as.integer(Delta)
  n <- length(X)

  if (length(L) != 1L || is.na(L) || L < 2L) {
    stop("Length must be a single integer greater than or equal to 2.")
  }

  if (n == 0L || anyNA(X) || any(!is.finite(X))) {
    stop("X must be a non-empty numeric vector with only finite values.")
  }

  if (length(delta) != 1L || is.na(delta) || delta < 1L) {
    stop("Delta must be a single positive integer.")
  }

  if (L * delta > n) {
    stop("The model is infeasible because Length * Delta is greater than length(X).")
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required.")
  }

  if (!requireNamespace("gurobi", quietly = TRUE)) {
    stop(
      paste(
        "The gurobi R package is required.",
        "Install it from the R directory of your Gurobi installation."
      )
    )
  }

  u_vals <- sort(X)
  indicators <- outer(
    u_vals,
    X,
    FUN = function(u, x) as.numeric(x <= u)
  )

  n_z <- n * L
  n_s <- n * L
  n_cdf <- n * L
  n_t <- n * L
  n_k <- L
  n_b <- n * L * n
  n_diff <- n * L

  z_start <- 1L
  s_start <- z_start + n_z
  cdf_start <- s_start + n_s
  t_start <- cdf_start + n_cdf
  k_start <- t_start + n_t
  b_start <- k_start + n_k
  diff_start <- b_start + n_b
  num_vars <- diff_start + n_diff - 1L

  z_idx <- function(i, l) z_start + (l - 1L) * n + i - 1L
  s_idx <- function(u, l) s_start + (l - 1L) * n + u - 1L
  cdf_idx <- function(u, l) cdf_start + (l - 1L) * n + u - 1L
  t_idx <- function(i, l) t_start + (l - 1L) * n + i - 1L
  k_idx <- function(l) k_start + l - 1L
  b_idx <- function(i, l, u) {
    b_start + (u - 1L) * n * L + (l - 1L) * n + i - 1L
  }
  diff_idx <- function(u, l) diff_start + (l - 1L) * n + u - 1L

  row_ids <- integer(0)
  col_ids <- integer(0)
  coefficients <- numeric(0)
  senses <- character(0)
  rhs_values <- numeric(0)
  row_count <- 0L

  add_constraint <- function(indices, values, sense, rhs) {
    row_count <<- row_count + 1L

    if (length(indices) != length(values)) {
      stop("Internal error: constraint indices and values have different lengths.")
    }

    keep <- values != 0

    if (any(keep)) {
      row_ids <<- c(row_ids, rep.int(row_count, sum(keep)))
      col_ids <<- c(col_ids, indices[keep])
      coefficients <<- c(coefficients, values[keep])
    }

    senses <<- c(senses, sense)
    rhs_values <<- c(rhs_values, rhs)
  }

  # Each observation must be assigned to exactly one segment.
  for (i in seq_len(n)) {
    indices <- vapply(seq_len(L), function(l) z_idx(i, l), integer(1))
    add_constraint(indices, rep(1, L), "=", 1)
  }

  # Each segment must contain at least Delta observations.
  for (l in seq_len(L)) {
    indices <- vapply(seq_len(n), function(i) z_idx(i, l), integer(1))
    add_constraint(indices, rep(1, n), ">", delta)
  }

  # Segment assignments must be monotone over the sequence.
  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      for (l in seq_len(L)) {
        later_segments <- l:L
        indices <- c(
          z_idx(i, l),
          vapply(later_segments, function(lp) z_idx(i + 1L, lp), integer(1))
        )
        values <- c(1, rep(-1, length(later_segments)))
        add_constraint(indices, values, "<", 0)
      }
    }
  }

  # Constraints from the original Python implementation for the first two
  # segments. These are why Length must be at least 2.
  for (i in seq_len(n)) {
    for (j in i:n) {
      add_constraint(c(z_idx(i, 1L), z_idx(j, 1L)), c(1, -1), ">", 0)
      add_constraint(c(z_idx(i, 2L), z_idx(j, 2L)), c(1, -1), "<", 0)
    }
  }

  for (l in seq_len(L)) {
    # The t values for each segment sum to one.
    t_indices <- vapply(seq_len(n), function(i) t_idx(i, l), integer(1))
    add_constraint(t_indices, rep(1, n), "=", 1)

    for (u in seq_len(n)) {
      active <- which(indicators[u, ] != 0)

      # cdf_l[u, l] = sum_i indicator[u, i] * t[i, l]
      add_constraint(
        c(cdf_idx(u, l), vapply(active, function(i) t_idx(i, l), integer(1))),
        c(1, rep(-1, length(active))),
        "=",
        0
      )

      # t[u, l] <= z[u, l]
      add_constraint(
        c(t_idx(u, l), z_idx(u, l)),
        c(1, -1),
        "<",
        0
      )

      # t[u, l] <= k[l] + z[u, l] / n - 1 / n
      add_constraint(
        c(t_idx(u, l), k_idx(l), z_idx(u, l)),
        c(1, -1, -1 / n),
        "<",
        -1 / n
      )

      # t[u, l] >= z[u, l] / n
      add_constraint(
        c(t_idx(u, l), z_idx(u, l)),
        c(1, -1 / n),
        ">",
        0
      )

      # t[u, l] >= k[l] + z[u, l] - 1
      add_constraint(
        c(t_idx(u, l), k_idx(l), z_idx(u, l)),
        c(1, -1, -1),
        ">",
        -1
      )
    }
  }

  for (l in seq_len(L)) {
    for (u in seq_len(n)) {
      for (i in seq_len(n)) {
        # b[i, l, u] <= diff[u, l]
        add_constraint(
          c(b_idx(i, l, u), diff_idx(u, l)),
          c(1, -1),
          "<",
          0
        )

        # b[i, l, u] <= z[i, l]
        add_constraint(
          c(b_idx(i, l, u), z_idx(i, l)),
          c(1, -1),
          "<",
          0
        )

        # b[i, l, u] >= z[i, l] + diff[u, l] - 1
        add_constraint(
          c(b_idx(i, l, u), z_idx(i, l), diff_idx(u, l)),
          c(1, -1, -1),
          ">",
          -1
        )
      }
    }
  }

  for (l in seq_len(L)) {
    for (u in seq_len(n)) {
      # diff[u, l] + cdf_l[u, l] = 1
      add_constraint(
        c(diff_idx(u, l), cdf_idx(u, l)),
        c(1, 1),
        "=",
        1
      )

      active <- which(indicators[u, ] != 0)

      # s[u, l] = sum_i indicator[u, i] * b[i, l, u]
      add_constraint(
        c(s_idx(u, l), vapply(active, function(i) b_idx(i, l, u), integer(1))),
        c(1, rep(-1, length(active))),
        "=",
        0
      )
    }
  }

  model <- list()
  model$modelname <- "Segmented_CDF_Diff"
  model$A <- Matrix::sparseMatrix(
    i = row_ids,
    j = col_ids,
    x = coefficients,
    dims = c(row_count, num_vars)
  )
  s_columns <- unlist(
    lapply(seq_len(L), function(l) {
      vapply(seq_len(n), function(u) s_idx(u, l), integer(1))
    }),
    use.names = FALSE
  )

  model$obj <- rep(0, num_vars)
  model$obj[s_columns] <- 1
  model$modelsense <- "min"
  model$rhs <- rhs_values
  model$sense <- senses
  model$lb <- rep(0, num_vars)
  model$ub <- rep(Inf, num_vars)
  model$vtype <- rep("C", num_vars)

  z_columns <- unlist(
    lapply(seq_len(L), function(l) {
      vapply(seq_len(n), function(i) z_idx(i, l), integer(1))
    }),
    use.names = FALSE
  )

  b_columns <- unlist(
    lapply(seq_len(n), function(u) {
      unlist(
        lapply(seq_len(L), function(l) {
          vapply(seq_len(n), function(i) b_idx(i, l, u), integer(1))
        }),
        use.names = FALSE
      )
    }),
    use.names = FALSE
  )

  cdf_columns <- unlist(
    lapply(seq_len(L), function(l) {
      vapply(seq_len(n), function(u) cdf_idx(u, l), integer(1))
    }),
    use.names = FALSE
  )

  diff_columns <- unlist(
    lapply(seq_len(L), function(l) {
      vapply(seq_len(n), function(u) diff_idx(u, l), integer(1))
    }),
    use.names = FALSE
  )

  k_columns <- vapply(seq_len(L), k_idx, integer(1))

  model$vtype[z_columns] <- "B"
  model$ub[z_columns] <- 1
  model$ub[b_columns] <- 1
  model$ub[cdf_columns] <- 1
  model$ub[diff_columns] <- 1
  model$ub[k_columns] <- 1 / delta

  params <- list(OutputFlag = as.integer(output_flag))
  result <- gurobi::gurobi(model, params)

  if (!identical(result$status, "OPTIMAL")) {
    stop(
      paste0(
        "Gurobi did not return an optimal solution. Status: ",
        result$status
      )
    )
  }

  z_sol <- matrix(result$x[z_columns], nrow = n, ncol = L)
  segment_sizes <- colSums(z_sol)

  cumsum(segment_sizes)
}
