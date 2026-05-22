# Monte Carlo simulation for the ordinal VRP estimator
#
# Specification: z = x  (mirrors the ANES application)
#
# Factors varied:
#   miss_p in {0.5, 0.65, 0.8}
#   rho    in {0.0, 0.3, 0.6}
#   N      fixed at N_sample (default 2000)
#
# Estimators compared:
#   naive   - ordered probit on respondents only (MAR, rho=0)
#   binary  - Peress (2010) binary VRP on median-split outcome
#   ordinal - proposed ordinal VRP (vrpoprob)
#
# True DGP parameters calibrated from life-satisfaction estimates.


source("code/vrpoprob.R")

tic("Simulation total")

# config ------------------------------------------------------------------
# All of these can be overridden by setting them before source()-ing this file.

if (!exists("N_sample"))    N_sample    <- 2000
if (!exists("N_rep"))       N_rep       <- 500
if (!exists("miss_grid"))   miss_grid   <- c(0.5, 0.65, 0.8)
if (!exists("rho_grid"))    rho_grid    <- c(0.0, 0.3, 0.6)
if (!exists("N_cores"))     N_cores     <- max(1, parallel::detectCores() - 1)
if (!exists("seed_base"))   seed_base   <- 20260422
if (!exists("save_output")) save_output <- TRUE

if (!exists("res_array")) load("results/results.RData")


# calibration from life-satisfaction, miss_p = 0.65 -----------------------

calib       <- res_array[["0.65"]][[1]]
alpha_core  <- calib$alpha
beta_core   <- calib$beta
lambda_true <- calib$lambda
theta_true  <- calib$theta

M <- length(lambda_true) + 1
R <- length(theta_true)


# population grid ---------------------------------------------------------

Xpop_s  <- data.matrix(Xpop)
Wpop_s <- Wpop


# helpers -----------------------------------------------------------------

calibrate_miss_shift <- function(miss_target) {
  f <- function(b0) sum(Wpop_s * pnorm(Xpop_s %*% beta_core + b0)) - miss_target
  uniroot(f, c(-15, 15))$root
}

true_pphat <- function(alpha_use, lambda_use) {
  ystarhat <- as.numeric(Xpop_s %*% alpha_use)
  P <- matrix(0, nrow(Xpop_s), M)
  P[, 1] <- pnorm(lambda_use[1] - ystarhat)
  if (M > 2) {
    for (m in 2:(M - 1))
      P[, m] <- pnorm(lambda_use[m] - ystarhat) - pnorm(lambda_use[m - 1] - ystarhat)
  }
  P[, M] <- 1 - pnorm(lambda_use[M - 1] - ystarhat)
  colSums(P * Wpop_s)
}


# one simulated dataset ---------------------------------------------------

sim_one_dataset <- function(N, rho, miss_p) {

  b0      <- calibrate_miss_shift(miss_p)
  theta_t <- theta_true - b0

  idx_x    <- sample(seq_len(nrow(Xpop_s)), N, replace = TRUE, prob = Wpop_s)
  X_sample <- Xpop_s[idx_x, , drop = FALSE]

  Sigma  <- matrix(c(1, rho, rho, 1), 2, 2)
  shocks <- mnorm::rmnorm(N, mean = c(0, 0), sigma = Sigma)
  eps <- shocks[, 1]; eta <- shocks[, 2]

  ystar <- as.numeric(X_sample %*% alpha_core) + eps
  rstar <- as.numeric(X_sample %*% beta_core)  + eta

  y          <- findInterval(ystar, c(-Inf, lambda_true))
  is_nonresp <- rstar >= theta_t[R]
  r          <- findInterval(rstar, c(-Inf, theta_t[seq_len(R - 1)]))

  keep  <- !is_nonresp
  Nmiss <- sum(is_nonresp)

  list(
    ydata   = y[keep],
    rdata   = pmin(pmax(r[keep], 1), R),
    xdata   = X_sample[keep, , drop = FALSE],
    Nmiss   = Nmiss,
    pp_true = true_pphat(alpha_core, lambda_true),
    rho_true = rho
  )
}


# naive estimator ---------------------------------------------------------

estimate_naive <- function(D) {
  x_cols <- colnames(D$xdata)
  if (is.null(x_cols)) x_cols <- paste0("x", seq_len(ncol(D$xdata)))
  df_fit <- data.frame(y = factor(D$ydata, levels = 1:M), D$xdata)
  colnames(df_fit)[-1] <- x_cols
  fml <- as.formula(paste0("y ~ ", paste(x_cols, collapse = " + ")))
  fit <- try(MASS::polr(fml, data = df_fit, method = "probit"), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, M))
  new_df <- as.data.frame(Xpop_s)
  colnames(new_df) <- x_cols
  preds <- predict(fit, newdata = new_df, type = "probs")
  if (is.null(dim(preds))) preds <- matrix(preds, nrow = 1)
  full <- matrix(0, nrow(preds), M)
  cn <- as.integer(colnames(preds))
  full[, cn] <- preds
  as.numeric(colSums(full * Wpop_s))
}


# binary VRP estimator ----------------------------------------------------

estimate_binary <- function(D) {
  med <- median(D$ydata)
  if (med >= M) med <- M - 1
  y_bin <- ifelse(D$ydata <= med, 1L, 2L)
  if (length(unique(y_bin)) < 2)
    return(list(pphat_bin = c(NA, NA), rho = NA, converged = FALSE, split = med))
  fit <- try(suppressWarnings(vrpoprob_estim(
    ydata = y_bin,  rdata = D$rdata,
    xdata = D$xdata, zdata = D$xdata,
    Nmiss = D$Nmiss,
    Wpop = Wpop_s, Xpop = Xpop_s,
    Zpop = Xpop_s
  )), silent = TRUE)
  if (inherits(fit, "try-error") || !isTRUE(fit$converged))
    return(list(pphat_bin = c(NA, NA), rho = NA, converged = FALSE, split = med))
  list(pphat_bin = fit$pphat, rho = fit$rho, converged = TRUE, split = med)
}


# ordinal VRP estimator ---------------------------------------------------

estimate_ordinal <- function(D) {
  fit <- try(suppressWarnings(vrpoprob_estim(
    ydata = D$ydata,  rdata = D$rdata,
    xdata = D$xdata,  zdata = D$xdata,
    Nmiss = D$Nmiss,
    Wpop = Wpop_s, Xpop = Xpop_s,
    Zpop = Xpop_s
  )), silent = TRUE)
  if (inherits(fit, "try-error") || !isTRUE(fit$converged))
    return(list(pphat = rep(NA_real_, M), pphat_se = rep(NA_real_, M),
                rho = NA_real_, converged = FALSE))
  list(pphat = fit$pphat, pphat_se = fit$pphat_se,
       rho = fit$rho, converged = TRUE)
}


# one replication ---------------------------------------------------------

run_replication <- function(rep_id, rho, miss_p, N) {
  seed <- seed_base + 1000L * rep_id +
    7L * which(rho_grid  == rho) +
         which(miss_grid == miss_p)
  set.seed(seed)

  D         <- sim_one_dataset(N, rho, miss_p)
  out_naive <- estimate_naive(D)
  out_bin   <- estimate_binary(D)
  out_ord   <- estimate_ordinal(D)

  list(
    rep = rep_id, rho = rho, miss_p = miss_p, N = N,
    realized_miss  = D$Nmiss / (D$Nmiss + length(D$ydata)),
    pp_true        = D$pp_true,
    naive_pp       = out_naive,
    bin_pp         = out_bin$pphat_bin,
    bin_rho        = out_bin$rho,
    bin_split      = out_bin$split,
    bin_converged  = out_bin$converged,
    ord_pp         = out_ord$pphat,
    ord_pp_se      = out_ord$pphat_se,
    ord_rho        = out_ord$rho,
    ord_converged  = out_ord$converged
  )
}


# main loop ---------------------------------------------------------------

grid <- expand.grid(rep_id = seq_len(N_rep),
                    rho    = rho_grid,
                    miss_p = miss_grid,
                    stringsAsFactors = FALSE,
                    KEEP.OUT.ATTRS   = FALSE)

message("Replications: ", nrow(grid), "  Cores: ", N_cores)

sim_results <- mclapply(seq_len(nrow(grid)), function(i) {
  row <- grid[i, ]
  tryCatch(
    run_replication(row$rep_id, row$rho, row$miss_p, N_sample),
    error = function(e) list(error = conditionMessage(e),
                             rep = row$rep_id, rho = row$rho, miss_p = row$miss_p)
  )
}, mc.cores = N_cores, mc.preschedule = FALSE)

if (save_output) {
  save(sim_results, grid,
       N_sample, N_rep, miss_grid, rho_grid,
       alpha_core, beta_core, lambda_true, theta_true,
       file = "results/simulation.RData")
  message("Saved to results/simulation.RData")
}


# ------------------------------------------------------------------
# Summarize results -> simulation_table1.csv, simulation_table2.csv
# ------------------------------------------------------------------

summarize_cell <- function(sub) {

  n_total   <- length(sub)
  n_ord_ok  <- sum(vapply(sub, function(r) isTRUE(r$ord_converged), logical(1)))
  n_bin_ok  <- sum(vapply(sub, function(r) isTRUE(r$bin_converged), logical(1)))

  pp_true  <- sub[[1]]$pp_true
  rho_true <- sub[[1]]$rho

  med_split <- round(stats::median(vapply(sub, function(r)
    if (is.null(r$bin_split)) NA_real_ else r$bin_split, numeric(1)), na.rm = TRUE))
  bin_truth <- c(sum(pp_true[seq_len(med_split)]),
                 sum(pp_true[(med_split + 1):M]))

  pp_naive <- do.call(rbind, lapply(sub, `[[`, "naive_pp"))
  pp_ord   <- do.call(rbind, lapply(sub, `[[`, "ord_pp"))
  pp_bin   <- do.call(rbind, lapply(sub, `[[`, "bin_pp"))

  cover_ord <- mean(vapply(seq_len(M), function(m) {
    pt <- vapply(sub, function(r) r$ord_pp[m],    numeric(1))
    se <- vapply(sub, function(r) r$ord_pp_se[m], numeric(1))
    ok <- is.finite(pt) & is.finite(se) & se > 0
    lo <- pt - 1.96 * se;  hi <- pt + 1.96 * se
    mean((pp_true[m] >= lo & pp_true[m] <= hi)[ok])
  }, numeric(1)))

  L1_naive   <- sum(abs(colMeans(pp_naive, na.rm = TRUE) - pp_true))
  L1_ord     <- sum(abs(colMeans(pp_ord,   na.rm = TRUE) - pp_true))
  RMSE_naive <- mean(sqrt(colMeans((pp_naive - matrix(pp_true, nrow(pp_naive), M, byrow = TRUE))^2, na.rm = TRUE)))
  RMSE_ord   <- mean(sqrt(colMeans((pp_ord   - matrix(pp_true, nrow(pp_ord),   M, byrow = TRUE))^2, na.rm = TRUE)))

  ord_bin   <- cbind(rowSums(pp_ord[,   seq_len(med_split), drop = FALSE]),
                     rowSums(pp_ord[,   (med_split + 1):M,  drop = FALSE]))
  naive_bin <- cbind(rowSums(pp_naive[, seq_len(med_split), drop = FALSE]),
                     rowSums(pp_naive[, (med_split + 1):M,  drop = FALSE]))
  L1_naive_bin <- sum(abs(colMeans(naive_bin, na.rm = TRUE) - bin_truth))
  L1_bin       <- sum(abs(colMeans(pp_bin,    na.rm = TRUE) - bin_truth))
  L1_ord_bin   <- sum(abs(colMeans(ord_bin,   na.rm = TRUE) - bin_truth))

  rho_ord <- Filter(is.finite, vapply(sub, `[[`, numeric(1), "ord_rho"))
  rho_bin <- Filter(is.finite, vapply(sub, `[[`, numeric(1), "bin_rho"))

  list(
    n_total = n_total, n_ord_ok = n_ord_ok, n_bin_ok = n_bin_ok,
    L1_naive = L1_naive,   L1_ord = L1_ord,
    RMSE_naive = RMSE_naive, RMSE_ord = RMSE_ord,
    L1_naive_bin = L1_naive_bin, L1_bin = L1_bin, L1_ord_bin = L1_ord_bin,
    cover_ord    = cover_ord,
    rho_true     = rho_true,
    rho_ord_bias = mean(rho_ord) - rho_true,
    rho_ord_rmse = sqrt(mean((rho_ord - rho_true)^2)),
    rho_bin_bias = if (length(rho_bin) > 0) mean(rho_bin) - rho_true else NA,
    rho_bin_rmse = if (length(rho_bin) > 0) sqrt(mean((rho_bin - rho_true)^2)) else NA
  )
}

summary_grid <- list()
for (mp in miss_grid) for (rr in rho_grid) {
  key <- paste(mp, rr, sep = "_")
  sub <- sim_results[vapply(sim_results, function(r)
    !is.null(r) && is.null(r$error) && r$miss_p == mp && r$rho == rr,
    logical(1))]
  if (length(sub) == 0) next
  summary_grid[[key]] <- summarize_cell(sub)
}

save(summary_grid, file = "results/simulation_summary.RData")

# Table 1: Naive vs Ordinal VRP -- all 9 cells
rows1 <- list()
for (mp in miss_grid) for (rr in rho_grid) {
  cs <- summary_grid[[paste(mp, rr, sep = "_")]]
  if (is.null(cs)) next
  rows1[[length(rows1) + 1L]] <- data.frame(
    rho = rr, miss_p = mp,
    L1_bias_naive  = cs$L1_naive,    L1_bias_ord   = cs$L1_ord,
    RMSE_naive     = cs$RMSE_naive,   RMSE_ord      = cs$RMSE_ord,
    rho_ord_bias   = cs$rho_ord_bias, rho_ord_rmse  = cs$rho_ord_rmse,
    coverage       = cs$cover_ord,
    conv_rate      = cs$n_ord_ok / cs$n_total,
    stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, rows1), "results/simulation_table1.csv", row.names = FALSE)

# Table 2: Binary vs Ordinal VRP -- nonignorable cells (rho > 0)
rows2 <- list()
for (mp in miss_grid) for (rr in rho_grid[rho_grid > 0]) {
  cs <- summary_grid[[paste(mp, rr, sep = "_")]]
  if (is.null(cs)) next
  rows2[[length(rows2) + 1L]] <- data.frame(
    rho = rr, miss_p = mp,
    L1_naive_bin = cs$L1_naive_bin,
    L1_bin       = cs$L1_bin,
    L1_ord_bin   = cs$L1_ord_bin,
    rho_bin_bias = cs$rho_bin_bias, rho_bin_rmse = cs$rho_bin_rmse,
    rho_ord_bias = cs$rho_ord_bias, rho_ord_rmse = cs$rho_ord_rmse,
    stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, rows2), "results/simulation_table2.csv", row.names = FALSE)

message("Wrote simulation_table1.csv and simulation_table2.csv")


toc()
