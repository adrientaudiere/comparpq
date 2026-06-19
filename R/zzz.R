.onLoad <- function(libname, pkgname) {
  # dabestr 2025.3.15 incompatibility with boot::boot.ci's BCA interval path:
  # boot.ci's BCA computation calls empinf -> inf.jack -> statistic(data, weights, ...)
  # but dabestr's internal bootboot statistic expects row indices, not frequency
  # weights, so it produces garbage that makes the acceleration 'a' = NA.
  # Fix: replace get_boot_row to use type = "perc" only, setting bca_ci_* = NA.
  if (requireNamespace("dabestr", quietly = TRUE)) {
    tryCatch(
      {
        ns <- getNamespace("dabestr")
        effsize_boot_fn <- get("effsize_boot", envir = ns, inherits = FALSE)
        var_w_df_fn <- get("var_w_df", envir = ns, inherits = FALSE)
        calc_grp_var_fn <- get(
          "calculate_group_variance",
          envir = ns,
          inherits = FALSE
        )

        new_get_boot_row <- function(
          ctrl_measurement,
          test_measurement,
          effect_size_func,
          seed,
          reps,
          is_paired,
          control_group,
          test_group,
          ci
        ) {
          control_test <- list(
            control = ctrl_measurement,
            test = test_measurement
          )
          ctrl_size <- length(ctrl_measurement)
          ctrl_var <- var_w_df_fn(ctrl_measurement, ctrl_size)
          test_size <- length(test_measurement)
          test_var <- var_w_df_fn(test_measurement, test_size)
          grp_var <- calc_grp_var_fn(
            ctrl_var = ctrl_var,
            ctrl_N = ctrl_size,
            test_var = test_var,
            test_N = test_size
          )
          weight <- 1 / grp_var
          set.seed(seed)
          boots <- effsize_boot_fn(
            data = control_test,
            effect_size_func = effect_size_func,
            reps = reps,
            paired = is_paired
          )
          bootci <- boot::boot.ci(boots, conf = ci / 100, type = "perc")
          list(
            control_group = control_group,
            test_group = test_group,
            bootstraps = list(as.vector(boots$t)),
            nboots = length(boots$t),
            bca_ci_low = NA_real_,
            bca_ci_high = NA_real_,
            pct_ci_low = bootci$percent[4],
            pct_ci_high = bootci$percent[5],
            ci = ci,
            difference = boots$t0,
            weight = weight
          )
        }

        utils::assignInNamespace(
          "get_boot_row",
          new_get_boot_row,
          ns = "dabestr"
        )
      },
      error = function(e) NULL
    )
  }
}
