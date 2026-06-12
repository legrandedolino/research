##################################################################
#
# Original version by Le Dolino
# Current version by Justin Tuazon
#
##################################################################

# 1. Preliminaries
# 1.1 Dependencies
install_and_load <- function(...) {
  packages <- unlist(list(...))
  present <-
    unlist(lapply(packages, require, character.only = TRUE))
  needed <- packages[present == FALSE]
  if (length(needed) > 0) {
    install.packages(needed, repos = "https://cran.r-project.org/")
    lapply(needed, require, character.only = TRUE)
  }
  for (package in packages) {
    library(package, character.only = TRUE)
  }
}

install_and_load("lme4", "lmerTest", "optimx", "rockchalk", "openxlsx",
                 "mvtnorm", "performance", "pracma")

# 1.2 Working directory
# Change this
setwd("~/Documents/Longitudinal/New")

# 2. Defining functions
# 2.1 Data generator
generate_data <- function(N = 10, T = 3, ICC = 0.1, var_epsilon = 1, b_00 = 0, 
                          b_10 = 0, b_11_prop = 0.1, tau_01_sq = 0, pE = 0.5, 
                          dist = "normal", attrition = 0, nonsamp_p = 0,
                          nonsamp_q = 0) {
  var_int <- (var_epsilon * ICC) / ((1 - ICC) * (1 + ((T - 1) * (2 * T - 1)) / 6))
  # b_11 is in terms of sd not var to be on same scale
  # no variance for slope fixed effect because at time = 0
  b_11 <- b_11_prop * sqrt(var_int + var_epsilon) 
  if (dist == "normal") {
    r_ijs <- rnorm(N * T, sd = sqrt(var_epsilon))
    if (tau_01_sq == 0) {
      u_0is <- rep(rnorm(N, sd = sqrt(var_int)), each = T) 
      u_1is <- rep(rnorm(N, sd = sqrt(var_int)), each = T)
    } else {
      if (abs(tau_01_sq) > var_int) {
        # from the Cauchy-Schwarz inequality, |Cov(X, Y)| <= sqrt(Var(X)Var(Y))
        stop("tau_01_sq must be at most var_int")
      }
      sigma <- matrix(c(var_int, tau_01_sq, tau_01_sq, var_int), nrow = 2, 
                      byrow = TRUE)
      # rockchalk::mvnrnorm is a wrapper for MASS:mvrnorm, the rockchalk
      # version follows random seeding correctly while MASS version does not
      us <- rockchalk::mvrnorm(N, rep(0, 2), sigma)
      u_0is <- rep(us[, 1], each = T)
      u_1is <- rep(us[, 2], each = T)
    }
  } else if (dist == "unif") {
    # The variance of U(a, b) is (1/12) * (b - a)^2 and its mean is 
    # (1 / 2) * (a + b). Since we set the mean to 0, it follows that a = -b so
    # that the variance is (1 / 12) * (4b^2) = (1 / 3) * b^2, where b > 0.
    # Thus, if the variance is v, it follows that b = sqrt(3v).
    r_ijs <- runif(N * T, -sqrt(3 * var_epsilon), sqrt(3 * var_epsilon))
    if (tau_01_sq == 0) {
      u_0is <- rep(runif(N, -sqrt(3 * var_int), sqrt(3 * var_int)), each = T)
      u_1is <- rep(runif(N, -sqrt(3 * var_int), sqrt(3 * var_int)), each = T)
    } else {
      if (abs(tau_01_sq) > var_int) {
        # from the Cauchy-Schwarz inequality, |Cov(X, Y)| <= sqrt(Var(X)Var(Y))
        stop("tau_01_sq must be at most var_int")
      }
      # See https://ieeexplore.ieee.org/abstract/document/5408380
      resid_corr <- tau_01_sq / var_int
      z_corr <- 2 * (sin((pi / 6) * resid_corr))
      sigma <- matrix(c(1, z_corr, z_corr, 1), nrow = 2, byrow = TRUE)
      us <- rockchalk::mvrnorm(N, rep(0, 2), sigma)
      u_0is <- rep(us[, 1], each = T)
      u_1is <- rep(us[, 2], each = T)
      # The desired distribution is U(-sqrt(3 * var_int), sqrt(3 * var_int)) 
      # because this has mean 0 and variance var_int.
      # Then, the inverse CDF of this distribution is 
      # -sqrt(3 * var_int) + 2 * x * sqrt(3 * var_int).
      u_inv <- function(x) {
        return(-sqrt(3 * var_int) + 2 * x * sqrt(3 * var_int))
      }
      u_0is <- u_inv(pnorm(u_0is))
      u_1is <- u_inv(pnorm(u_1is))
      # The correlation between the residuals is resid_corr.
    }
    
  }
  data <- data.frame(
    subject = as.factor(rep(1:N, each = T)),
    time = rep(0:(T - 1), N),
    r_ij = r_ijs,
    u_0i = u_0is,
    u_1i = u_1is,
    treatment = 0
  )
  treatment_group <- sample(N, N * pE)
  data[data$subject %in% treatment_group, ]$treatment <- 1
  data$Y <- b_00 + b_10 * data$time + b_11 * data$treatment * data$time + 
    data$r_ij + data$u_0i + data$u_1i * data$time
  if (attrition > 0) {
    # The math is in the supplementary material, Simulating Attrition.
    last_periods <- rep(NULL, N)
    for (i in 1:N) {
      bernoullis <- rbinom(T - 1, 1, 1 - attrition)
      min_bernoulli <- suppressWarnings(min(which(bernoullis == 1)))
      last_periods[i] <- if (min_bernoulli == Inf) 0 else (T - min_bernoulli)
    }
    data <- data[data$time <= last_periods[data$subject], ]
    rownames(data) <- 1:nrow(data)
  }
  if (nonsamp_p > 0 & nonsamp_q > 0) {
    # The math is in the supplementary material, Simulating Non-sampling Errors.
    obs <- 1:nrow(data)
    contaminated_group <- sample(obs, nrow(data) * nonsamp_p)
    bound <- sqrt(3 * nonsamp_q * (var_int + var_epsilon))
    contaminants <- runif(length(contaminated_group), -bound, bound)
    for (i in contaminated_group) {
      data$Y[i] <- data$Y[i] + contaminants[i]
    }
  }
  return(data)
}

# 2.2 Model simulator
simulate_model <- function(N = 10, T = 3, ICC = 0.1, var_epsilon = 1, b_00 = 0, 
                           b_10 = 0, b_11_prop = 0.1, tau_01_sq = 0, pE = 0.5,
                           alpha = 0.05, dist = "normal", attrition = 0, 
                           nonsamp_p = 0, nonsamp_q = 0, niter = 1000, 
                           seed = NA, input_num = NA) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  param_metrics <- list(
    b_11 = list(bias = rep(NA, niter), reject_null = rep(NA, niter), noncov = rep(NA, niter)),
    b_00 = list(bias = rep(NA, niter), reject_null = rep(NA, niter), noncov = rep(NA, niter)),
    b_10 = list(bias = rep(NA, niter), reject_null = rep(NA, niter), noncov = rep(NA, niter)),
    tau_00_sq = list(bias = rep(NA, niter), reject_null = rep(NA, niter), noncov = rep(NA, niter)),
    tau_11_sq = list(bias = rep(NA, niter), reject_null = rep(NA, niter), noncov = rep(NA, niter)),
    sigma_sq = list(bias = rep(NA, niter), reject_null = rep(NA, niter), noncov = rep(NA, niter))
  )
  opt_problem <- rep(NA, niter)
  conv_problem <- rep(NA, niter)
  rand_ci_problem <- rep(NA, niter)
  soln_problem <- rep(NA, niter)
  conv_boundary_problem <- rep(NA, niter)
  conv_gradient_evaluation_problem <- rep(NA, niter)
  conv_gradient_converge_fail_problem <- rep(NA, niter)
  conv_hessian_check_problem <- rep(NA, niter)
  conv_hessian_degenerate_problem <- rep(NA, niter)
  conv_hessian_singular_problem <- rep(NA, niter)
  conv_eigenvalue_large_problem <- rep(NA, niter)
  inadequacy <- rep(1, niter)
  models <- list()
  for (i in 1:niter) {
    print(sprintf("Input %s - Iteration %s", input_num, i))
    models[[paste("model", i, sep = "_")]] <- NA
    data <- generate_data(N, T, ICC, var_epsilon, b_00, b_10, b_11_prop, tau_01_sq, pE, 
                          dist, attrition, nonsamp_p, nonsamp_q)
    var_int <- (var_epsilon * ICC) / ((1 - ICC) * (1 + ((T - 1) * (2 * T - 1)) / 6))
    b_11 <- b_11_prop * sqrt(var_int + var_epsilon)
    tryCatch(
      {
        # Optimizer minqa::bobyqa appears to be the fastest
        # https://cran.r-project.org/web/packages/lme4/vignettes/lmerperf.html
        # By default, max iterations is 10,000 from source code of
        # bobyqa from the minqa package.
        model <- lmerTest::lmer(Y ~ time + treatment:time + 
                                  (1 + time || subject), data = data, 
                                REML = TRUE, control = lmerControl(
                                  optimizer = "bobyqa"
                                ))
        if (!is.null(attr(getME(model, "X"), "col.dropped"))) {
          # If there is not enough data to fit the model or estimate it,
          # then the data is inadequate (i.e., inadequacy == 1)
          stop()
        }
        models[[paste("model", i, sep = "_")]] <- model
        inadequacy[i] <- 0
        summ <- summary(model, ddf = "Kenward-Roger")[["coefficients"]]
        rands <- as.data.frame(VarCorr(model))
        with_rand_ci <- TRUE
        tryCatch({
          intervals <- confint(model, method = "profile", level = 1 - alpha,
                               parm = "theta_", oldNames = FALSE)
          rand_ci_problem[i] <- 0
        }, error = function(e) {
          print("- - - - -")
          print(e)
          print("This error was caught! The error is at the first CI level.")
          print("- - - - -")
        })
        if (is.na(rand_ci_problem[i])) {
          tryCatch({
            pp <- profile(model, which = "theta_", devtol = Inf, 
                          signames = FALSE)
            intervals <- confint(pp, method = "profile", level = 1 - alpha, 
                                 parm = "theta_", oldNames = FALSE)
            rand_ci_problem[i] <- 0
          }, error = function(e) {
            print("- - - - -")
            print(e)
            print("This error was caught! The error is at the second CI level.")
            print("- - - - -")
          })
          with_rand_ci <- !is.na(rand_ci_problem[i])
          rand_ci_problem[i] <- 1
        }
        
        # b_11
        vals <- summ["time:treatment", ]
        est <- vals[["Estimate"]]
        me <- qt(p = alpha / 2, df = vals[["df"]], 
                 lower.tail = FALSE) * vals[["Std. Error"]]
        p_val <- vals[["Pr(>|t|)"]]
        param_metrics$b_11$bias[i] <- est - b_11
        param_metrics$b_11$reject_null[i] <- as.integer(p_val <= alpha)
        param_metrics$b_11$noncov[i] <- as.integer((b_11 < est - me) || (b_11 > est + me))
        
        # b_00
        vals <- summ["(Intercept)", ]
        est <- vals[["Estimate"]]
        me <- qt(p = alpha / 2, df = vals[["df"]], 
                 lower.tail = FALSE) * vals[["Std. Error"]]
        p_val <- vals[["Pr(>|t|)"]]
        param_metrics$b_00$bias[i] <- est - b_00
        param_metrics$b_00$reject_null[i] <- as.integer(p_val <= alpha)
        param_metrics$b_00$noncov[i] <- as.integer((b_00 < est - me) || (b_00 > est + me))
        
        # b_10
        vals <- summ["time", ]
        est <- vals[["Estimate"]]
        me <- qt(p = alpha / 2, df = vals[["df"]], 
                 lower.tail = FALSE) * vals[["Std. Error"]]
        p_val <- vals[["Pr(>|t|)"]]
        param_metrics$b_10$bias[i] <- est - b_10
        param_metrics$b_10$reject_null[i] <- as.integer(p_val <= alpha)
        param_metrics$b_10$noncov[i] <- as.integer((b_10 < est - me) || (b_10 > est + me))
        
        # tau_00_sq
        est <- as.numeric(rands[rands$grp == "subject", "vcov"])
        true_val <- var_int
        param_metrics$tau_00_sq$bias[i] <- est - true_val
        if (with_rand_ci) {
          ci <- intervals["sd_(Intercept)|subject", ]
          if (!(is.na(ci[1]) || is.na(ci[2]))) {
            param_metrics$tau_00_sq$reject_null[i] <- as.integer(0 < ci[1])
            param_metrics$tau_00_sq$noncov[i] <- as.integer(
              sqrt(true_val) < ci[1] || sqrt(true_val) > ci[2]
            )
          } else {
            rand_ci_problem[i] <- 1
          }
        }

        # tau_11_sq
        est <- as.numeric(rands[rands$grp == "subject.1", "vcov"])
        true_val <- var_int
        param_metrics$tau_11_sq$bias[i] <- est - true_val
        if (with_rand_ci) {
          ci <- intervals["sd_time|subject", ]
          if (!(is.na(ci[1]) || is.na(ci[2]))) {
            param_metrics$tau_11_sq$reject_null[i] <- as.integer(0 < ci[1])
            param_metrics$tau_11_sq$noncov[i] <- as.integer(
              sqrt(true_val) < ci[1] || sqrt(true_val) > ci[2]
            )
          } else {
            rand_ci_problem[i] <- 1
          }
        }

        # sigma_sq
        est <- as.numeric(rands[rands$grp == "Residual", "vcov"])
        true_val <- var_epsilon
        param_metrics$sigma_sq$bias[i] <- est - true_val
        if (with_rand_ci) {
          ci <- intervals["sigma", ]
          if (!(is.na(ci[1]) || is.na(ci[2]))) {
            param_metrics$sigma_sq$reject_null[i] <- as.integer(0 < ci[1])
            param_metrics$sigma_sq$noncov[i] <- as.integer(
              sqrt(true_val) < ci[1] || sqrt(true_val) > ci[2]
            )
          } else {
            rand_ci_problem[i] <- 1
          }
        }

        # Checks if there is an error or warning returned by the 
        # optimization algorithm
        opt_problem[i] <- as.integer(model@optinfo$conv$opt != 0)
        
        # Checks if there is a problem with the gradient, Hessian, or
        # eigenvalues based on the Karush-Kuhn-Tucker conditions and also 
        # whether the fit is singular (e.g., est. var. is exactly 0).
        conv_boundary_problem[i] <- as.integer(isSingular(model))
        conv_problem[i] <- as.integer(
          (!is.null(model@optinfo$conv$lme4$code)) | 
            conv_boundary_problem[i]
        )
        
        # Checks whether there is a problem with the solution in general
        soln_problem[i] <- as.integer(opt_problem[i] | conv_problem[i])
        
        if (conv_problem[i] == 1) {
          # The conditions here are based on the source code for lme4
          # particularly, checkConv.R
          # https://github.com/lme4/lme4/blob/master/R/checkConv.R
          msgs <- model@optinfo$conv$lme4$messages
          conv_gradient_evaluation_problem[i] <- as.integer(any(grepl(
            "unable to evaluate scaled gradient",
            msgs,
            ignore.case = TRUE
          )))
          conv_gradient_converge_fail_problem[i] <- as.integer(any(
            grepl(
              "Model failed to converge with max",
              msgs,
              ignore.case = TRUE
            )
          ))
          conv_hessian_check_problem[i] <- as.integer(any(grepl(
            "Problem with Hessian check",
            msgs,
            ignore.case = TRUE
          )))
          conv_hessian_degenerate_problem[i] <- as.integer(any(grepl(
            "degenerate",
            msgs,
            ignore.case = TRUE
          )))
          conv_hessian_singular_problem[i] <- as.integer(any(grepl(
            "Hessian is numerically singular",
            msgs,
            ignore.case = TRUE
          )))
          conv_eigenvalue_large_problem[i] <- as.integer(any(grepl(
            "Model is nearly unidentifiable",
            msgs,
            ignore.case = TRUE
          )))
        } else {
          conv_gradient_evaluation_problem[i] <- 0
          conv_gradient_converge_fail_problem[i] <- 0
          conv_hessian_check_problem[i] <- 0
          conv_hessian_degenerate_problem[i] <- 0
          conv_hessian_singular_problem[i] <- 0
          conv_eigenvalue_large_problem[i] <- 0
        }
      },
      error = function(e) {
        print("- - - - -")
        print(e)
        print("This error was caught! The error is at the model-fitting level.")
        print("- - - - -")
      }
    )
  }
  results <- list(
    "iterations" = data.frame(
      iteration = 1:niter,
      
      b_11_bias = param_metrics$b_11$bias,
      b_11_reject_null = param_metrics$b_11$reject_null,
      b_11_noncov = param_metrics$b_11$noncov,
      
      b_00_bias = param_metrics$b_00$bias,
      b_00_reject_null = param_metrics$b_00$reject_null,
      b_00_noncov = param_metrics$b_00$noncov,
      
      b_10_bias = param_metrics$b_10$bias,
      b_10_reject_null = param_metrics$b_10$reject_null,
      b_10_noncov = param_metrics$b_10$noncov,
      
      tau_00_sq_bias = param_metrics$tau_00_sq$bias,
      tau_00_sq_reject_null = param_metrics$tau_00_sq$reject_null,
      tau_00_sq_noncov = param_metrics$tau_00_sq$noncov,
      
      tau_11_sq_bias = param_metrics$tau_11_sq$bias,
      tau_11_sq_reject_null = param_metrics$tau_11_sq$reject_null,
      tau_11_sq_noncov = param_metrics$tau_11_sq$noncov,
      
      sigma_sq_bias = param_metrics$sigma_sq$bias,
      sigma_sq_reject_null = param_metrics$sigma_sq$reject_null,
      sigma_sq_noncov = param_metrics$sigma_sq$noncov,
      
      soln_problem = soln_problem,
      opt_problem = opt_problem,
      conv_problem = conv_problem,
      conv_boundary_problem = conv_boundary_problem,
      conv_gradient_evaluation_problem = conv_gradient_evaluation_problem,
      conv_gradient_converge_fail_problem = conv_gradient_converge_fail_problem, 
      conv_hessian_check_problem = conv_hessian_check_problem,
      conv_hessian_degenerate_problem = conv_hessian_degenerate_problem,
      conv_hessian_singular_problem = conv_hessian_singular_problem,
      conv_eigenvalue_large_problem = conv_eigenvalue_large_problem,
      inadequacy = inadequacy,
      rand_ci_problem = rand_ci_problem
    ),
    "models" = models
  )
  # iterations contains the results for each iteration
  # models contains the model objects (for troubleshooting or inspection)
  return(results)
}

# 2.3 Simulation summarizer
summarize_simulation <- function(results) {
  # This gets the proportions for the metrics from the return value of 
  # simulate_model
  iterations <- results$iterations
  summary <- list(
    "b_11_bias" = mean(iterations$b_11_bias, na.rm = TRUE),
    "b_11_rmse" = sqrt(mean((iterations$b_11_bias)^2, na.rm = TRUE)),
    "b_11_power" = mean(iterations$b_11_reject_null, na.rm = TRUE),
    "b_11_noncov" = mean(iterations$b_11_noncov, na.rm = TRUE),
    
    "b_00_bias" = mean(iterations$b_00_bias, na.rm = TRUE),
    "b_00_rmse" = sqrt(mean((iterations$b_00_bias)^2, na.rm = TRUE)),
    "b_00_power" = mean(iterations$b_00_reject_null, na.rm = TRUE),
    "b_00_noncov" = mean(iterations$b_00_noncov, na.rm = TRUE),
    
    "b_10_bias" = mean(iterations$b_10_bias, na.rm = TRUE),
    "b_10_rmse" = sqrt(mean((iterations$b_10_bias)^2, na.rm = TRUE)),
    "b_10_power" = mean(iterations$b_10_reject_null, na.rm = TRUE),
    "b_10_noncov" = mean(iterations$b_10_noncov, na.rm = TRUE),
    
    "tau_00_sq_bias" = mean(iterations$tau_00_sq_bias, na.rm = TRUE),
    "tau_00_sq_rmse" = sqrt(mean((iterations$tau_00_sq_bias)^2, na.rm = TRUE)),
    "tau_00_sq_power" = mean(iterations$tau_00_sq_reject_null, na.rm = TRUE),
    "tau_00_sq_noncov" = mean(iterations$tau_00_sq_noncov, na.rm = TRUE),
    
    "tau_11_sq_bias" = mean(iterations$tau_11_sq_bias, na.rm = TRUE),
    "tau_11_sq_rmse" = sqrt(mean((iterations$tau_11_sq_bias)^2, na.rm = TRUE)),
    "tau_11_sq_power" = mean(iterations$tau_11_sq_reject_null, na.rm = TRUE),
    "tau_11_sq_noncov" = mean(iterations$tau_11_sq_noncov, na.rm = TRUE),
    
    "sigma_sq_bias" = mean(iterations$sigma_sq_bias, na.rm = TRUE),
    "sigma_sq_rmse" = sqrt(mean((iterations$sigma_sq_bias)^2, na.rm = TRUE)),
    "sigma_sq_power" = mean(iterations$sigma_sq_reject_null, na.rm = TRUE),
    "sigma_sq_noncov" = mean(iterations$sigma_sq_noncov, na.rm = TRUE),
    
    "soln_problem" = mean(iterations$soln_problem, na.rm = TRUE),
    "opt_problem" = mean(iterations$opt_problem, na.rm = TRUE),
    "conv_problem" = mean(iterations$conv_problem, na.rm = TRUE),
    "conv_boundary_problem" = mean(iterations$conv_boundary_problem, na.rm = TRUE),
    "conv_gradient_evaluation_problem" = mean(
      iterations$conv_gradient_evaluation_problem, na.rm = TRUE
    ),
    "conv_gradient_converge_fail_problem" = mean(
      iterations$conv_gradient_converge_fail_problem, na.rm = TRUE
    ),
    "conv_hessian_check_problem" = mean(
      iterations$conv_hessian_check_problem, na.rm = TRUE
    ),
    "conv_hessian_degenerate_problem" = mean(
      iterations$conv_hessian_degenerate_problem, na.rm = TRUE
    ),
    "conv_hessian_singular_problem" = mean(
      iterations$conv_hessian_singular_problem, na.rm = TRUE
    ),
    "conv_eigenvalue_large_problem" = mean(
      iterations$conv_eigenvalue_large_problem, na.rm = TRUE
    ),
    "inadequacy" = mean(iterations$inadequacy),
    "rand_ci_problem" = mean(iterations$rand_ci_problem)
  )
  return(summary)
}

# 3. Running simulations
# 3.1 Inputs
inputs <- read.xlsx("inputs_full.xlsx", "inputs")

# 3.2 Simulations
start_time <- Sys.time()
name <- rep("input", nrow(inputs))
power <- rep(NA, nrow(inputs))
param_metrics <- list(
  b_11 = list(bias = rep(NA, nrow(inputs)), rmse = rep(NA, nrow(inputs)), 
              power = rep(NA, nrow(inputs)), noncov = rep(NA, nrow(inputs))),
  b_00 = list(bias = rep(NA, nrow(inputs)), rmse = rep(NA, nrow(inputs)), 
              power = rep(NA, nrow(inputs)), noncov = rep(NA, nrow(inputs))),
  b_10 = list(bias = rep(NA, nrow(inputs)), rmse = rep(NA, nrow(inputs)), 
              power = rep(NA, nrow(inputs)), noncov = rep(NA, nrow(inputs))),
  tau_00_sq = list(bias = rep(NA, nrow(inputs)), rmse = rep(NA, nrow(inputs)), 
                   power = rep(NA, nrow(inputs)), noncov = rep(NA, nrow(inputs))),
  tau_11_sq = list(bias = rep(NA, nrow(inputs)), rmse = rep(NA, nrow(inputs)), 
                   power = rep(NA, nrow(inputs)), noncov = rep(NA, nrow(inputs))),
  sigma_sq = list(bias = rep(NA, nrow(inputs)), rmse = rep(NA, nrow(inputs)), 
                  power = rep(NA, nrow(inputs)), noncov = rep(NA, nrow(inputs)))
)
soln_problem <- rep(NA, nrow(inputs))
opt_problem <- rep(NA, nrow(inputs))
conv_problem <- rep(NA, nrow(inputs))
conv_boundary_problem <- rep(NA, nrow(inputs))
conv_gradient_evaluation_problem <- rep(NA, nrow(inputs))
conv_gradient_converge_fail_problem <- rep(NA, nrow(inputs))
conv_hessian_check_problem <- rep(NA, nrow(inputs))
conv_hessian_degenerate_problem <- rep(NA, nrow(inputs))
conv_hessian_singular_problem <- rep(NA, nrow(inputs))
conv_eigenvalue_large_problem <- rep(NA, nrow(inputs))
inadequacy <- rep(NA, nrow(inputs))
rand_ci_problem <- rep(NA, nrow(inputs))
for (i in 1:nrow(inputs)) {
  name[i] <- paste(name[i], i, sep = "_")
  inp <- inputs[i, ]

  tmp <- summarize_simulation(simulate_model(
    N = inp[["N"]], T = inp[["T"]], ICC = inp[["ICC"]],
    var_epsilon = inp[["var_epsilon"]], b_00 = inp[["b_00"]],
    b_10 = inp[["b_10"]], b_11_prop = inp[["b_11_prop"]], tau_01_sq = inp[["tau_01_sq"]],
    pE = inp[["pE"]], alpha = inp[["alpha"]], dist = inp[["dist"]],
    attrition = inp[["attrition"]], nonsamp_p = inp[["nonsamp_p"]],
    nonsamp_q = inp[["nonsamp_q"]], niter = inp[["niter"]],
    seed = inp[["seed"]], input_num = i
  ))

  param_metrics$b_11$bias[i] <- tmp$b_11_bias
  param_metrics$b_11$rmse[i] <- tmp$b_11_rmse
  param_metrics$b_11$power[i] <- tmp$b_11_power
  param_metrics$b_11$noncov[i] <- tmp$b_11_noncov

  param_metrics$b_00$bias[i] <- tmp$b_00_bias
  param_metrics$b_00$rmse[i] <- tmp$b_00_rmse
  param_metrics$b_00$power[i] <- tmp$b_00_power
  param_metrics$b_00$noncov[i] <- tmp$b_00_noncov

  param_metrics$b_10$bias[i] <- tmp$b_10_bias
  param_metrics$b_10$rmse[i] <- tmp$b_10_rmse
  param_metrics$b_10$power[i] <- tmp$b_10_power
  param_metrics$b_10$noncov[i] <- tmp$b_10_noncov

  param_metrics$tau_00_sq$bias[i] <- tmp$tau_00_sq_bias
  param_metrics$tau_00_sq$rmse[i] <- tmp$tau_00_sq_rmse
  param_metrics$tau_00_sq$power[i] <- tmp$tau_00_sq_power
  param_metrics$tau_00_sq$noncov[i] <- tmp$tau_00_sq_noncov

  param_metrics$tau_11_sq$bias[i] <- tmp$tau_11_sq_bias
  param_metrics$tau_11_sq$rmse[i] <- tmp$tau_11_sq_rmse
  param_metrics$tau_11_sq$power[i] <- tmp$tau_11_sq_power
  param_metrics$tau_11_sq$noncov[i] <- tmp$tau_11_sq_noncov

  param_metrics$sigma_sq$bias[i] <- tmp$sigma_sq_bias
  param_metrics$sigma_sq$rmse[i] <- tmp$sigma_sq_rmse
  param_metrics$sigma_sq$power[i] <- tmp$sigma_sq_power
  param_metrics$sigma_sq$noncov[i] <- tmp$sigma_sq_noncov

  soln_problem[i] <- tmp$soln_problem
  opt_problem[i] <- tmp$opt_problem
  conv_problem[i] <- tmp$conv_problem
  conv_boundary_problem[i] <- tmp$conv_boundary_problem
  conv_gradient_evaluation_problem[i] <- tmp$conv_gradient_evaluation_problem
  conv_gradient_converge_fail_problem[i] <- tmp$conv_gradient_converge_fail_problem
  conv_hessian_check_problem[i] <- tmp$conv_hessian_check_problem
  conv_hessian_degenerate_problem[i] <- tmp$conv_hessian_degenerate_problem
  conv_hessian_singular_problem[i] <- tmp$conv_hessian_singular_problem
  conv_eigenvalue_large_problem[i] <- tmp$conv_eigenvalue_large_problem
  inadequacy[i] <- tmp$inadequacy
  rand_ci_problem[i] <- tmp$rand_ci_problem
}
names <- data.frame(
  name = name
)
outputs <- data.frame(
  b_11_bias = param_metrics$b_11$bias,
  b_11_rmse = param_metrics$b_11$rmse,
  b_11_power = param_metrics$b_11$power,
  b_11_noncov = param_metrics$b_11$noncov,
  
  b_00_bias = param_metrics$b_00$bias,
  b_00_rmse = param_metrics$b_00$rmse,
  b_00_power = param_metrics$b_00$power,
  b_00_noncov = param_metrics$b_00$noncov,
  
  b_10_bias = param_metrics$b_10$bias,
  b_10_rmse = param_metrics$b_10$rmse,
  b_10_power = param_metrics$b_10$power,
  b_10_noncov = param_metrics$b_10$noncov,
  
  tau_00_sq_bias = param_metrics$tau_00_sq$bias,
  tau_00_sq_rmse = param_metrics$tau_00_sq$rmse,
  tau_00_sq_power = param_metrics$tau_00_sq$power,
  tau_00_sq_noncov = param_metrics$tau_00_sq$noncov,
  
  tau_11_sq_bias = param_metrics$tau_11_sq$bias,
  tau_11_sq_rmse = param_metrics$tau_11_sq$rmse,
  tau_11_sq_power = param_metrics$tau_11_sq$power,
  tau_11_sq_noncov = param_metrics$tau_11_sq$noncov,
  
  sigma_sq_bias = param_metrics$sigma_sq$bias,
  sigma_sq_rmse = param_metrics$sigma_sq$rmse,
  sigma_sq_power = param_metrics$sigma_sq$power,
  sigma_sq_noncov = param_metrics$sigma_sq$noncov,
  
  soln_problem = soln_problem,
  opt_problem = opt_problem,
  conv_problem = conv_problem,
  conv_boundary_problem = conv_boundary_problem,
  conv_gradient_evaluation_problem = conv_gradient_evaluation_problem,
  conv_gradient_converge_fail_problem = conv_gradient_converge_fail_problem,
  conv_hessian_check_problem = conv_hessian_check_problem,
  conv_hessian_degenerate_problem = conv_hessian_degenerate_problem,
  conv_hessian_singular_problem = conv_hessian_singular_problem,
  conv_eigenvalue_large_problem = conv_eigenvalue_large_problem,
  inadequacy = inadequacy,
  rand_ci_problem = rand_ci_problem
)
outputs <- cbind(names, inputs, outputs)
write.xlsx(outputs, "outputs.xlsx")
end_time <- Sys.time()
print(end_time - start_time)
