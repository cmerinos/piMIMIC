#' @title Mahalanobis Distance for DIF Detection with PI-MIMIC
#'
#' @description
#' Calculates the Mahalanobis distance for each item based on selected DIF metrics
#' standardized expected parameter change (SEPC), expected parameter change (EPC) or delta R^2). Items with large
#' distances are potential outliers and may indicate substantial DIF. A plot can be generated to visualize
#' the distances.
#'
#' @param output A list returned by \code{\link{piMIMIC}} or \code{\link{piMIMIClrt}}.
#' @param method Character. Either `"score"` (for output from `piMIMIC`) or `"lrt"` (for output from `piMIMIClrt`).
#' @param target Character. Which DIF dimension to evaluate:
#'   - `"global"`: uses both uniform and non-uniform effects (two dimensions, only available for `method = "score"`).
#'   - `"udif"`: uses only uniform DIF effects.
#'   - `"nudif"`: uses only non-uniform DIF effects.
#' @param metric Character. Only for `method = "score"`. Either `"epc"` (unstandardized) or `"stdepc"` (standardized, default).
#' @param alpha Numeric. Significance level for the chi-square cutoff (default = 0.95). The cutoff is `qchisq(alpha, df)`.
#' @param plot Logical. If `TRUE`, generates a bar plot of the distances with a cutoff line.
#' @param cut.line Numeric. Optional custom cutoff value. If `NULL` (default), it is calculated from `alpha` and the number of dimensions.
#' @param ... Additional arguments passed to `ggplot2` (e.g., `theme`).
#'
#' @return A list with components:
#'   \item{distances}{A data frame with items and their Mahalanobis distances (sorted descending).}
#'   \item{cutoff}{The cutoff value used.}
#'   \item{plot}{A ggplot object (if `plot = TRUE`).}
#'
#' @details
#' **Mahalanobis distance**
#' It measures how far each item is from the center of the distribution of all items,
#' taking into account the correlations between the DIF metrics. A large distance
#' indicates that comparatively the item behaves unusually show strong DIF.
#'
#' **When you have only one metric** (e.g., `target = "udif"`), the Mahalanobis
#' distance reduces to the squared z-score: `(x - mean)^2 / var`. This follows a
#' chi-square distribution with 1 degree of freedom, so the cutoff is
#' `qchisq(alpha, df = 1)`.
#'
#' **Two dimensions (`target = "global"`):**
#' The distance is calculated in a bivariate space. The covariance matrix is estimated
#' from the data. For the distance to be reliable, you need at least 3 items
#' (preferably more). Still with fewer than 5 items, the covariance matrix may be unstable,
#' and the distances should be interpreted with caution.
#'
#' **Choosing the cutoff:**
#' The default cutoff is the `alpha` quantile of the chi-square distribution with
#' `df = number of dimensions` (1 or 2). Items with distance above the cutoff
#' are flagged as potential outliers. You can also supply a custom value via
#' `cut.line`.
#'
#' @examples
#' \dontrun{
#' ### Loading database ####
#' library(psych)
#' #Loading data
#' data("bfi")
#'
#' # Choosing variables
#' data.bfi <- bfi[, c("N1", "N2", "N3", "N4", "N5", "gender", "age")]
#'
#' # Clean for missing values
#' data.bfi <- data.bfi[complete.cases(data.bfi), ]
#'
#' data.bfi$gender <- as.factor(data.bfi$gender)
#'
#' neuro.items <- c("N1", "N2", "N3", "N4", "N5")
#'
#' # Using piMIMIC output (score test)
#' result <- piMIMIC(data = data.bfi, items = neuro.items, cov = "gender")
#'
#' # Mahalanobis distance based on both uniform and non-uniform SEPC
#' dist_global <- distancepiMIMIC(output = result, method = "score",
#'                                target = "global", metric = "stdepc")
#' print(dist_global$distances)
#'
#' # Focus only on non-uniform DIF (one dimension)
#' dist_nu <- distancepiMIMIC(output = result, method = "score",
#'                            target = "nudif", metric = "stdepc",
#'                            plot = TRUE)
#'
#' # Using piMIMIClrt output (LRT) with delta R^2
#' result_lrt <- piMIMIClrt(data = bfi, items = neuro.items, cov = "gender",
#'                          anchor = "rest")
#' dist_lrt <- distancepiMIMIC(output = result_lrt, method = "lrt",
#'                             target = "global", plot = TRUE)
#' }
#'
#' @importFrom stats var cov mahalanobis qchisq complete.cases
#' @importFrom ggplot2 ggplot aes geom_col geom_hline coord_flip labs theme_minimal
#' @importFrom rlang sym .data
#'
#' @export
distancepiMIMIC <- function(output,
                            method = c("score", "lrt"),
                            target = c("global", "udif", "nudif"),
                            metric = c("stdepc", "epc"),
                            alpha = 0.95,
                            plot = FALSE,
                            cut.line = NULL,
                            ...) {

  # ---- Initial Validations ----
  method <- match.arg(method)
  target <- match.arg(target)
  if (method == "score") metric <- match.arg(metric)

  # ---- Extract metrics by method and target ----
  df_metrics <- NULL
  n_dim <- 0

  if (method == "score") {
    # Extract SEPC.uDIF and SEPC.nuDIF
    if (is.null(output$SEPC.uDIF) || is.null(output$SEPC.nuDIF)) {
      stop("Output does not contain SEPC.uDIF or SEPC.nuDIF. Did you run piMIMIC()?")
    }
    # Select a metric column
    col_metric <- if (metric == "stdepc") "SEPC.ALL" else "EPC"

    if (target %in% c("udif", "global")) {
      df_u <- output$SEPC.uDIF[, c("Item", col_metric)]
    }
    if (target %in% c("nudif", "global")) {
      df_nu <- output$SEPC.nuDIF[, c("Item", col_metric)]
    }

    if (target == "udif") {
      df_metrics <- df_u
      colnames(df_metrics)[2] <- "value"
      n_dim <- 1
    } else if (target == "nudif") {
      df_metrics <- df_nu
      colnames(df_metrics)[2] <- "value"
      n_dim <- 1
    } else if (target == "global") {
      # Combining both Dimensions
      df_metrics <- merge(df_u, df_nu, by = "Item", suffixes = c("_u", "_nu"))
      colnames(df_metrics)[2:3] <- c("udif", "nudif")
      n_dim <- 2
    }

  } else if (method == "lrt") {
    # Extract DeltaR^2
    if (is.null(output$DeltaR2.Global) && target == "global") {
      stop("Output does not contain DeltaR2.Global. Did you run piMIMIClrt()?")
    }
    if (is.null(output$DeltaR2.uDIF) && target == "udif") {
      stop("Output does not contain DeltaR2.uDIF. Did you run piMIMIClrt()?")
    }
    if (is.null(output$DeltaR2.nuDIF) && target == "nudif") {
      stop("Output does not contain DeltaR2.nuDIF. Did you run piMIMIClrt()?")
    }

    if (target == "udif") {
      df_metrics <- output$DeltaR2.uDIF
      colnames(df_metrics)[2] <- "value"
      n_dim <- 1
    } else if (target == "nudif") {
      df_metrics <- output$DeltaR2.nuDIF
      colnames(df_metrics)[2] <- "value"
      n_dim <- 1
    } else if (target == "global") {
      df_metrics <- output$DeltaR2.Global
      colnames(df_metrics)[2] <- "value"
      n_dim <- 1  # Just one column for the overall delta R^2
    }
  }

  # ---- Check to see if there is any data ----
  if (is.null(df_metrics) || nrow(df_metrics) == 0) {
    stop("No metrics extracted. Check your 'method' and 'target' arguments.")
  }

  # ---- Build a feature matrix ----
  if (n_dim == 1) {
    # A single column: a vector of values
    vals <- df_metrics$value
    # Remover NAs
    if (anyNA(vals)) {
      warning("NA values found in metrics. They were removed for distance calculation.")
      df_metrics <- df_metrics[!is.na(vals), ]
      vals <- vals[!is.na(vals)]
    }
    # Calculate the distance as z^2
    mu <- mean(vals)
    var_vals <- var(vals)
    if (var_vals == 0) stop("All values are identical. Cannot compute distance.")
    dist_vals <- (vals - mu)^2 / var_vals
    dist_df <- data.frame(Item = df_metrics$Item, Distance = dist_vals)
    # Sort in descending order
    dist_df <- dist_df[order(dist_df$Distance, decreasing = TRUE), ]
    # Cutoff
    if (is.null(cut.line)) {
      cut.line <- qchisq(alpha, df = 1)
    }

  } else if (n_dim == 2) {
    # Two columns: matrix 2D
    mat <- as.matrix(df_metrics[, c("udif", "nudif")])
    # Remove rows containing NA
    if (anyNA(mat)) {
      warning("NA values found in metrics. Rows with NA were removed.")
      complete <- complete.cases(mat)
      mat <- mat[complete, ]
      df_metrics <- df_metrics[complete, ]
    }
    if (nrow(mat) < 3) {
      warning("Less than 3 complete observations. Distance may be unstable.")
    }
    if (nrow(mat) < 2) {
      stop("At least 2 complete observations are needed for 2D Mahalanobis distance.")
    }
    # Calculate the mean and covariance
    mu <- colMeans(mat)
    Sigma <- cov(mat)
    # If the covariance matrix is singular, add a small constant to make it invertible
    # (this is a simple and transparent form of regularization)
    if (rcond(Sigma) < 1e-10) {
      message("Covariance matrix is near-singular. A small constant was added to the diagonal for invertibility.")
      eps <- 0.001 * mean(diag(Sigma))
      Sigma <- Sigma + diag(eps, nrow = ncol(Sigma))
    }
    # Distance
    dist_vals <- mahalanobis(mat, center = mu, cov = Sigma)
    dist_df <- data.frame(Item = df_metrics$Item, Distance = dist_vals)
    dist_df <- dist_df[order(dist_df$Distance, decreasing = TRUE), ]
    # Cutoff
    if (is.null(cut.line)) {
      cut.line <- qchisq(alpha, df = n_dim)
    }
  }

  # ---- Generate a chart if requested ----
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' needed for plotting. Returning distances only.")
    } else {
      # Sort items by distance (so the graph is in descending order)
      dist_df$Item <- factor(dist_df$Item, levels = rev(dist_df$Item))
      p <- ggplot2::ggplot(dist_df, ggplot2::aes(x = .data$Item, y = .data$Distance)) +
        ggplot2::geom_col(fill = "steelblue") +
        ggplot2::geom_hline(yintercept = cut.line, linetype = "dashed", color = "red") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = paste("Mahalanobis Distance -", target, "DIF"),
          x = "Item",
          y = "Mahalanobis Distance",
          caption = paste("Cutoff (", alpha*100, "%): ", round(cut.line, 3), sep = "")
        ) +
        ggplot2::theme_minimal()
      print(p)
    }
  }

  # ---- Output ----
  result <- list(
    distances = dist_df,
    cutoff = cut.line
  )
  if (plot && exists("p")) result$plot <- p

  return(result)
}
