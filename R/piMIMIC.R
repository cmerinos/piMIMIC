#' @title DIF analysis using PI-MIMIC with Score Test (Oort adjustment optional)
#'
#' @description
#' Implements the product indicator (PI) approach for MIMIC models to detect
#' uniform and non‑uniform DIF using the score test (modification indices).
#' Optionally applies Oort's critical value adjustment to control Type I error.
#'
#' @param data Data frame containing items and the covariate.
#' @param items Character vector of item names.
#' @param cov Name of the covariate (numeric or factor).
#' @param lvname Name of the latent variable (default "LatFact").
#' @param est Estimator for lavaan (default "MLM").
#' @param Oort.adj Logical; if `TRUE`, applies Oort's adjustment to the critical value.
#' @param p.crit Numeric; significance level for the Oort adjustment.
#'
#' @return A list with:
#' \item{DIF.Global}{Data frame with global DIF test (2 df): Chi², df, p-value, and Oort critical value if requested.}
#' \item{DIF.Uniforme}{Data frame with uniform DIF test (1 df).}
#' \item{DIF.NoUniforme}{Data frame with non‑uniform DIF test (1 df).}
#' \item{SEPC.uDIF}{Standardized expected parameter change for uniform DIF.}
#' \item{SEPC.nuDIF}{Standardized expected parameter change for non‑uniform DIF.}
#' \item{fit}{The fitted lavaan object.}
#'
#' @details
#' This function implements Differential Item Functioning (DIF) analysis using the Product of Indicators approach
#' (PI; Kolbe & Jorgensen, 2018) within the Restricted Factor Analysis (RFA; Oort, 1998) framework. This method operates
#' under a MIMIC scheme (Finch, 2005), incorporating latent variable interactions using PI; see Kolbe et al. (2018, 2019),
#' Kolbe, Jorgensen, & Molenaar (2020). It allows for the evaluation of uniform (uDIF) and non-uniform DIF (nuDIF)
#' with covariates that can be categorical (e.g., sex) or continuous (e.g., self-esteem, conscientiousness).
#'
#' Estimation is performed via lavaan::cfa, and DIF statistical tests are based on the Score test (Lagrange Multiplier
#' test used to evaluate fixed or constrained parameters). By default, chi-square tests compare an unrestricted model
#' (with covariate effect and interaction parameters freely estimated) and a restricted model (with these parameters
#' fixed to zero), using the standard chi-square distribution.
#'
#' Simulation studies (Oort, 1992, 1998; Kim, Yoon & Lee, 2011) have shown that
#' the LR test under MIMIC may suffer from inflated Type I error rates. To address this, Oort proposed
#' a correction of the chi-square critical value:
#'
#' \deqn{K^\prime = (\chi^2_0 / (K + df_0 - 1)) * K}
#' where \eqn{\chi^2_0} and \eqn{df_0} are from the baseline (full invariance) model, and \eqn{K} is
#' the original critical value. This adjustment is recommended when the baseline
#' model shows evidence of misfit (\eqn{\chi^2_0 / df_0 > 1}), as it helps control Type I error.
#'
#' where:
#' \itemize{
#'   \item \eqn{K^\prime} is the adjusted critical value,
#'   \item \eqn{K} is the original critical value from the chi-square distribution at significance level \eqn{p.crit},
#'   \item \eqn{\chi^2_0} is the chi-square statistic of the baseline model,
#'   \item \eqn{df_0} is the corresponding degrees of freedom of the baseline model.
#' }
#'
#' Interpretation:
#' \itemize{
#'   \item The reported \code{p.value} corresponds to the standard chi-square test.
#'   \item When \code{Oort.adj = TRUE}, the function also reports the adjusted critical value (\code{crit.Oort}).
#'   \item Users can compare the observed chi-square statistic against both thresholds:
#'         the conventional critical value (via \code{p.value}) and the Oort-adjusted critical value.
#'   \item This comparison allows assessment of how conclusions may differ when controlling for
#'         potential inflation of Type I error in MIMIC-PI models.
#'         }
#'
#' @examples
#' \donttest{
#' ### Example 1: simulated data -------------
#' set.seed(123)
#' Exmp1.data <- data.frame(
#'   grp = sample(0:1, 100, replace = TRUE),  # Group variable
#'   item1 = sample(1:5, 100, replace = TRUE),
#'   item2 = sample(1:5, 100, replace = TRUE),
#'   item3 = sample(1:5, 100, replace = TRUE),
#'   item4 = sample(1:5, 100, replace = TRUE))
#'
#' res1 <- piMIMIC(data = Exmp1.data , items = c("item1","item2","item3"), cov = "grp")
#' res1$DIF.Global
#'
#' ### Example 2: Using the 'bfi' dataset from the 'psych' package -------------
#' library(psych)
#' data("bfi")
#'
#' data.bfi <- bfi[, c("N1","N2","N3","N4","N5","gender")]
#'
#' data.bfi <- data.bfi[complete.cases(data.bfi), ]
#'
#' data.bfi$gender <- as.factor(data.bfi$gender)
#'
#' neuro.items <- c("N1","N2","N3","N4","N5")
#'
#' # Run DIF analysis with Oort adjustment
#' res.bfi <- piMIMIC(data = data.bfi, items = neuro.items, cov = "gender",
#'                  lvname = "Neuroticism", est = "MLM",
#'                  Oort.adj = TRUE, p.crit = 0.05)
#' res.bfi$DIF.Global
#' }
#'
#' @references
#' Kim, E. S., Yoon, M., & Lee, T. (2011). Testing Measurement Invariance Using MIMIC:
#' Likelihood Ratio Test With a Critical Value Adjustment. *Educational and Psychological Measurement, 72*(3), 469–492.
#' https://doi.org/10.1177/0013164411427395
#'
#' Oort, F. J. (1992). Using restricted factor analysis to detect item bias. *Psychological Methods*, 37, 547–567.
#'
#' Oort, F. J. (1998). Simulation study of item bias detection with restricted factor analysis.
#' *Structural Equation Modeling*, 5, 107–124.
#'
#' French, B. F., & Finch, W. H. (2008). Multigroup confirmatory factor analysis: Locating the invariant referent variables.
#' *Structural Equation Modeling*, 15(1), 96–113.
#'
#' Stark, S., Chernyshenko, O. S., & Drasgow, F. (2006). Detecting differential item functioning with confirmatory factor analysis and item response theory: Toward a unified strategy.
#' *Journal of Applied Psychology*, 91(6), 1292–1306.
#'
#' Kolbe, L., & Jorgensen, T. D. (2018). Using product indicators in restricted factor analysis models to detect nonuniform measurement bias.
#' In *Quantitative Psychology* (pp. 235–245). Springer.
#'
#' Kolbe, L., & Jorgensen, T. D. (2019). Using restricted factor analysis to select anchor items and detect differential item functioning.
#' *Behavior Research Methods*, 51, 138–151.
#'
#' Kolbe, L., Jorgensen, T. D., & Molenaar, D. (2020). The Impact of Unmodeled Heteroskedasticity on Assessing Measurement Invariance in Single-group Models.
#' *Structural Equation Modeling*, 28(1), 82–98.
#'
#' Whittaker, T. A. (2012). Estimation of Standardized Expected Parameter Change for DIF Detection.
#' *Educational and Psychological Measurement*, 72(3), 342-357.
#'
#' Garnier-Villarreal, M., & Jorgensen, T. D. (2024). Evaluating Local Model Misspecification with Modification Indices in Bayesian Structural Equation Modeling.
#' *Structural Equation Modeling*, 1–15.
#'
#' @importFrom lavaan cfa lavTestScore parameterestimates
#' @importFrom semTools indProd
#' @importFrom stats qchisq
#'
#' @export
piMIMIC <- function(data, items, cov, lvname = "LatFact", est = "MLM",
                    Oort.adj = FALSE, p.crit = 0.05) {

  # ---- Basic checks ----
  if (!is.character(est) || nchar(est) == 0) {
    stop("Error: Estimator must be a non-empty string (see lavaan documentation).")
  }
  if (!all(c(items, cov) %in% colnames(data))) {
    stop("Error: Some item names or the covariate are not present in the data frame.")
  }
  if (!is.character(lvname) || lvname == "") {
    stop("Error: The latent variable name ('lvname') must be a non-empty string.")
  }

  # ---- Prepare covariate: convert to numeric and center ----
  if (is.factor(data[[cov]])) {
    cov_num <- as.numeric(data[[cov]]) - 1
  } else {
    cov_num <- as.numeric(data[[cov]])
  }
  cov_centered <- cov_num - mean(cov_num, na.rm = TRUE)
  data[[paste0(cov, "_cent")]] <- cov_centered

  # ---- Create product indicators using semTools::indProd ----
  prod_data <- semTools::indProd(
    data = data,
    var1 = items,
    var2 = paste0(cov, "_cent"),
    doubleMC = TRUE,
    match = FALSE
  )
  prod_names <- paste0(items, ".", cov, "_cent")
  if (!all(prod_names %in% names(prod_data))) {
    stop("Product indicator columns were not created correctly.")
  }

  # ---- Build MIMIC model syntax ----
  cov_lat <- paste0(cov, "lat")
  int_fac <- paste0("LFacX", cov)

  syntax_lv <- paste0(lvname, " =~ ", paste(items, collapse = " + "))
  syntax_cov <- paste0(cov_lat, " =~ 1*", paste0(cov, "_cent"), "\n",
                       paste0(cov, "_cent"), " ~~ 0*", paste0(cov, "_cent"))
  syntax_int <- paste0(int_fac, " =~ ", paste(prod_names, collapse = " + "))
  residual_cov <- paste(paste0(items, " ~~ ", items, ".", cov, "_cent"), collapse = "\n")

  model_mimic <- paste(syntax_lv, syntax_cov, syntax_int, residual_cov, sep = "\n")

  # ---- Fit the MIMIC model ----
  fit <- lavaan::cfa(model_mimic,
                     data = prod_data,
                     estimator = est,
                     meanstructure = TRUE)

  # ---- Internal function: generate parameters for score tests ----
  generate_mimic_params <- function(items, cov_lat, int_fac) {
    flat_params <- character(length(items) * 2)
    grouped_params <- vector("list", length(items))
    names(grouped_params) <- items

    for (i in seq_along(items)) {
      flat_params[2*i - 1] <- paste0(items[i], " ~ ", cov_lat)
      flat_params[2*i]     <- paste0(items[i], " ~ ", int_fac)
      grouped_params[[i]] <- c(paste0(items[i], " ~ ", cov_lat),
                               paste0(items[i], " ~ ", int_fac))
    }
    return(list(flat = flat_params, grouped = grouped_params))
  }

  params <- generate_mimic_params(items, cov_lat, int_fac)

  # ---- Internal function: extract DIF results and SEPC ----
  mimicout_modificado <- function(fit.mimic, params, cov, Oort.adj, p.crit) {
    ests <- as.data.frame(lavaan::parameterestimates(fit.mimic))
    uniqnames <- unique(ests$lhs)
    lvname <- uniqnames[1]

    # ---- Global test (2 df per item) using grouped parameters ----
    # IMPORTANT: lavTestScore returns multiple rows if add has multiple parameters.
    # We only want the first row (the joint 2-df test) for each item.
    global_list <- list()
    for (i in seq_along(params$grouped)) {
      test_result <- lavaan::lavTestScore(fit.mimic, add = params$grouped[[i]])$test
      # Keep only the first row (joint test with 2 df)
      global_list[[i]] <- test_result[1, , drop = FALSE]
    }
    global_test <- do.call(rbind, global_list)
    rownames(global_test) <- NULL

    # ---- Univariate tests (1 df per parameter) using flat list ----
    uni_test <- lavaan::lavTestScore(fit.mimic, add = as.character(params$flat))

    # ---- Baseline chi2 and df for Oort adjustment ----
    baseline <- fit.mimic@test$standard
    chi0 <- as.numeric(baseline[["stat"]])
    df0  <- as.numeric(baseline[["df"]])

    if (Oort.adj) {
      K_global <- qchisq(1 - p.crit, 2)
      K_uniform <- qchisq(1 - p.crit, 1)
      crit.global <- (chi0 / (K_global + df0 - 1)) * K_global
      crit.uniform <- (chi0 / (K_uniform + df0 - 1)) * K_uniform
    }

    # ---- DIF.Global (now one row per item) ----
    df_dif_global <- data.frame(
      Item = items,
      Chi2 = round(global_test$X2, 3),
      df = global_test$df,
      p.value = round(global_test$p, 4)
    )
    if (Oort.adj) df_dif_global$crit.Oort <- round(crit.global, 3)

    # ---- DIF.Uniforme (odd indices: direct effects) ----
    oddnum <- seq(1, length(uni_test$uni$lhs), 2)
    df_dif_uniforme <- data.frame(
      Item = gsub("~.*$", "", uni_test$uni$lhs[oddnum]),
      Chi2 = round(uni_test$uni$X2[oddnum], 3),
      df = uni_test$uni$df[oddnum],
      p.value = round(uni_test$uni$p.value[oddnum], 4)
    )
    if (Oort.adj) df_dif_uniforme$crit.Oort <- round(crit.uniform, 3)

    # ---- DIF.NoUniforme (even indices: interaction effects) ----
    evennum <- seq(2, length(uni_test$uni$lhs), 2)
    df_dif_nouniforme <- data.frame(
      Item = gsub("~.*$", "", uni_test$uni$lhs[evennum]),
      Chi2 = round(uni_test$uni$X2[evennum], 3),
      df = uni_test$uni$df[evennum],
      p.value = round(uni_test$uni$p.value[evennum], 4)
    )
    if (Oort.adj) df_dif_nouniforme$crit.Oort <- round(crit.uniform, 3)

    # ---- SEPC ----
    sepc_values <- lavaan::lavTestScore(fit.mimic,
                                        add = as.character(params$flat),
                                        univariate = TRUE,
                                        standardized = TRUE,
                                        cov.std = TRUE,
                                        epc = TRUE)$epc

    df_sepc <- sepc_values[sepc_values$lhs %in% items &
                             sepc_values$rhs %in% c(cov_lat, int_fac),
                           c("lhs", "op", "rhs", "epc", "sepc.all")]

    colnames(df_sepc) <- c("Item", "Operator", "Effect", "EPC", "SEPC.ALL")
    df_sepc_u <- df_sepc[grepl(paste0(cov_lat, "$"), df_sepc$Effect), ]
    df_sepc_nu <- df_sepc[grepl(paste0("^", int_fac), df_sepc$Effect), ]

    return(list(
      DIF.Global = df_dif_global,
      DIF.Uniforme = df_dif_uniforme,
      DIF.NoUniforme = df_dif_nouniforme,
      SEPC.uDIF = df_sepc_u,
      SEPC.nuDIF = df_sepc_nu
    ))
  }

  # ---- Extract results ----
  resultados_DIF <- mimicout_modificado(fit, params, cov, Oort.adj, p.crit)
  resultados_DIF$fit <- fit
  return(resultados_DIF)
}
