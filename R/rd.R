# %% fuzzy regression discontinuity: build the instrument and the controls

# In a fuzzy RD the instrument is the indicator for crossing the cutoff. The
# running variable enters as controls, with a separate polynomial on each side,
# so that only the jump at the cutoff identifies the effect. A treatment
# "predicted" from the running variable is not a valid instrument by itself:
# it also moves smoothly with the running variable, and that smooth variation
# inflates the first stage.
#
# Returns the data restricted to the bandwidth with the new columns, and the
# names to use as Z, controls, FE and weights.
rd_prep <- function(data, running, cutoff = 0, Z = NULL, controls = NULL,
                    FE = NULL, weights = NULL, poly = 1, bandwidth = NULL,
                    kernel = "triangular") {
  if (length(running) != 1 || !running %in% names(data)) {
    stop("\"running\" must be one variable in the dataset.")
  }
  if (!is.numeric(cutoff) || length(cutoff) < 1 || any(!is.finite(cutoff))) {
    stop("\"cutoff\" must be one or more finite numbers.")
  }
  if (!poly %in% 0:4) stop("\"rd.poly\" must be 0, 1, 2, 3 or 4.")
  kernel <- match.arg(kernel, c("triangular", "uniform"))
  data <- as.data.frame(data)
  data <- data[!is.na(data[, running]), , drop = FALSE]
  r <- data[, running]
  # with several cutoffs, each unit is compared with its nearest cutoff
  cutoff <- sort(unique(cutoff))
  near <- sapply(r, function(v) which.min(abs(v - cutoff)))
  x <- r - cutoff[near]
  # bandwidth and kernel
  w <- rep(1, length(x))
  if (!is.null(bandwidth)) {
    if (!is.numeric(bandwidth) || length(bandwidth) != 1 || bandwidth <= 0) {
      stop("\"bandwidth\" must be one positive number.")
    }
    if (kernel == "triangular") {
      keep <- abs(x) < bandwidth
      w <- 1 - abs(x) / bandwidth
    } else {
      keep <- abs(x) <= bandwidth
    }
  } else {
    keep <- rep(TRUE, length(x))
  }
  data <- data[keep, , drop = FALSE]
  x <- x[keep]; w <- w[keep]; near <- near[keep]
  above <- as.numeric(x >= 0)
  if (length(unique(above)) < 2) stop("No observations on one side of the cutoff.")
  # the instrument
  if (is.null(Z)) {
    data$rd_above <- above
    Z <- "rd_above"
  } else {
    z <- data[, Z]
    step <- length(Z) == 1 && all(z %in% c(0, 1)) &&
      length(unique(z[x > 0])) == 1 && length(unique(z[x < 0])) == 1 &&
      z[x > 0][1] != z[x < 0][1]
    if (!isTRUE(step)) {
      message("Z is not the indicator for crossing the cutoff. Only the variation in Z ",
        "left after the controls, the running-variable polynomial (by side) and the ",
        "cutoff fixed effects identifies the effect; make sure that variation is the ",
        "jump at the cutoff, not smooth variation in the running variable.")
    }
  }
  # polynomial in the running variable, separately on each side
  rd.controls <- c()
  if (poly >= 1) {
    for (k in 1:poly) {
      v1 <- paste0("rd_x", k)
      v2 <- paste0("rd_x", k, "_above")
      data[, v1] <- x^k
      data[, v2] <- x^k * above
      rd.controls <- c(rd.controls, v1, v2)
    }
  }
  controls <- unique(c(controls, rd.controls))
  # several cutoffs: compare within cutoff
  if (length(cutoff) > 1) {
    data$rd_cutoff <- factor(cutoff[near])
    FE <- unique(c(FE, "rd_cutoff"))
  }
  # kernel weights (times the user's weights)
  if (!is.null(bandwidth) && kernel == "triangular") {
    if (is.null(weights)) {
      data$rd_w <- w
    } else {
      data$rd_w <- w * data[, weights]
    }
    weights <- "rd_w"
  }
  info <- list(running = running, cutoff = cutoff, poly = poly,
    bandwidth = bandwidth, kernel = if (is.null(bandwidth)) "uniform" else kernel,
    N_left = sum(above == 0), N_right = sum(above == 1))
  list(data = data, Z = Z, controls = controls, FE = FE, weights = weights, info = info)
}
