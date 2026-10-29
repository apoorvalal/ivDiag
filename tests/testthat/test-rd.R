# fuzzy RD option of ivDiag(): instrument = crossing the cutoff, running
# variable partialled out on each side

sim_rd <- function(n = 2000, seed = 1, jump = 0.4) {
  set.seed(seed)
  x <- runif(n, -1, 1)
  above <- as.numeric(x >= 0)
  u <- rnorm(n)
  # take-up rises smoothly with x and jumps at the cutoff
  d <- as.numeric(0.2 + 0.3 * x + jump * above + 0.3 * u + rnorm(n, sd = 0.3) > 0.5)
  y <- 2 * d + 1.5 * x - x^2 + u
  # a "predicted" treatment that also moves smoothly with x
  p_hat <- 0.2 + 0.3 * x + jump * above
  data.frame(y = y, d = d, x = x, p_hat = p_hat, g = rep(1:(n / 10), each = 10))
}

test_that("rd_prep builds the indicator, the polynomial by side and the kernel", {
  df <- sim_rd()
  rd <- ivDiag:::rd_prep(df, running = "x", cutoff = 0, poly = 2, bandwidth = 0.5)
  expect_equal(rd$Z, "rd_above")
  expect_setequal(rd$controls, c("rd_x1", "rd_x1_above", "rd_x2", "rd_x2_above"))
  expect_true(all(abs(rd$data$x) < 0.5))
  expect_equal(rd$weights, "rd_w")
  expect_equal(rd$data$rd_w, 1 - abs(rd$data$x) / 0.5)
  expect_equal(rd$data$rd_above, as.numeric(rd$data$x >= 0))
  expect_equal(rd$info$N_left + rd$info$N_right, nrow(rd$data))
})

test_that("several cutoffs: nearest cutoff and cutoff fixed effects", {
  df <- sim_rd()
  rd <- ivDiag:::rd_prep(df, running = "x", cutoff = c(-0.5, 0.5))
  expect_true("rd_cutoff" %in% rd$FE)
  expect_equal(rd$data$rd_x1, df$x - ifelse(df$x < 0, -0.5, 0.5))
})

test_that("errors and messages", {
  df <- sim_rd()
  expect_error(ivDiag:::rd_prep(df, running = "x", cutoff = 2), "one side")
  expect_error(ivDiag:::rd_prep(df, running = "nope"), "running")
  expect_message(ivDiag:::rd_prep(df, running = "x", Z = "p_hat"), "not the indicator")
  df$z_step <- as.numeric(df$x >= 0)
  expect_silent(ivDiag:::rd_prep(df, running = "x", Z = "z_step"))
})

test_that("ivDiag with a running variable recovers the effect", {
  df <- sim_rd(n = 4000, seed = 2)
  g <- suppressMessages(ivDiag(df, Y = "y", D = "d", running = "x", cutoff = 0,
    rd.poly = 2, cl = "g", bootstrap = FALSE, parallel = FALSE))
  expect_equal(g$est_2sls[1, "Coef"], 2, tolerance = 0.25)
  expect_false(is.null(g$rd))
  expect_equal(rownames(g$est_fs), "rd_above")
  # same numbers as building the columns by hand
  df$a <- as.numeric(df$x >= 0); df$x1 <- df$x; df$x1a <- df$x * df$a
  df$x2 <- df$x^2; df$x2a <- df$x^2 * df$a
  h <- suppressMessages(ivDiag(df, Y = "y", D = "d", Z = "a",
    controls = c("x1", "x1a", "x2", "x2a"), cl = "g", bootstrap = FALSE, parallel = FALSE))
  expect_equal(unname(g$est_2sls), unname(h$est_2sls))
  expect_equal(g$AR$ci, h$AR$ci)
})

test_that("a predicted treatment without the running variable overstates the first stage", {
  df <- sim_rd(n = 4000, seed = 3)
  naive <- suppressMessages(ivDiag(df, Y = "y", D = "d", Z = "p_hat",
    bootstrap = FALSE, run.AR = FALSE, parallel = FALSE))
  rd <- suppressMessages(ivDiag(df, Y = "y", D = "d", running = "x", cutoff = 0,
    bootstrap = FALSE, run.AR = FALSE, parallel = FALSE))
  expect_gt(naive$F_stat["F.effective"], 2 * rd$F_stat["F.effective"])
})
