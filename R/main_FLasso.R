#' GP-Based Pareto Front for Fused Lasso
#'
#' This function performs a Pareto optimality search on the Fused Lasso problem
#' using GP-based optimization (`easyGParetoptim`) and creates a 3D scatter plot
#' of the solution space. The plot displays the number of non-zero coefficients,
#' residual sum of squares (RSS), and roughness (mean absolute difference of consecutive beta coefficients) for different Pareto-optimal combinations
#'  of the regularization parameters \eqn{\lambda_1} and \eqn{\lambda_2}.
#'
#' @param X A numeric matrix of predictor variables (n x p), where n is the number
#' of samples and p is the number of predictors.
#' @param y A numeric vector of length n representing the response variable.
#' @param lb A numeric vector specifying the lower bounds (log10-scale) for \eqn{\lambda_1} and \eqn{\lambda_2}. Default is c(-3, -3).
#' @param ub A numeric vector specifying the upper bounds (log10-scale) for \eqn{\lambda_1} and \eqn{\lambda_2}. Default is c(1, 1).
#' @param Pareto_budget Number of evaluations in the GP-based Pareto optimization. Default is 100.
#' @param plot_title A character string for the title of the plot. Default is "Fused LASSO".
#' @param plot_marker_color Fill color of markers. Default is red with transparency.
#' @param plot_marker_border_color Border color of markers. Default is black.
#' @param plot_marker_size Size of Pareto front markers. Default value is 10.
#' @param plot_marker_border_width Width of marker borders. Default value is 2.
#' @param plot_marker_symbol Symbol shape for plot markers (e.g., 'circle', 'square', 'diamond'). Default shape is 'circle'.
#' @param fontsize_x Font size of the x-axis title. Default value is 20.
#' @param fontsize_y Font size of the y-axis title. Default value is 20.
#' @param fontsize_z Font size of the z-axis title. Default value is 20.
#' @param plot_title_size Font size for the plot title. Default value is 40.
#' @param title_pos_x Horizontal position of the title (0 to 1 scale). Default is 0.5.
#' @param title_pos_y Vertical position of the title (0 to 1 scale). Default is 0.95.
#' @param draw_projection Logical flag (0 or 1) to add projections from each point to the XY plane. Default is 1.
#' @param proj_line_color Color for projection lines. Default is black.
#' @param proj_line_width Width of projection lines. Default value is 2.
#' @param proj_marker_size Size of projection plane markers. Default value is 2.
#' @param proj_marker_color Fill color of projection markers. Default is red with transparency.
#' @param proj_marker_border_width Border width of projection markers. Default value is 2.
#' @param proj_marker_border_color Border color of projection markers. Default is black.
#'
#' @return A list with two componenets:
#' \itemize{
#'   \item \code{figure}: A \code{plotly} 3D scatter plot object showing the Pareto front.
#'   \item \code{summary_table}: A data frame containing the summary of the Pareto
#'   optimal solutions, including the values for \eqn{\lambda_1}, \eqn{\lambda_2},
#'   number of non-zero coefficients, RSS, and roughness.
#' }
#'
#' @details
#' The Fused Lasso objective function used is:
#' \deqn{
#' \mathcal{L}(\beta) = \frac{1}{2n} \|y - X\beta\|_2^2 + \lambda_1 \|\beta\|_1 + \lambda_2 \sum_{j=2}^{p} |\beta_j - \beta_{j-1}|
#' }
#' where \eqn{X} is the design matrix, \eqn{y} is the response vector,
#' \eqn{\beta} is the coefficient vector, \eqn{\lambda_1 > 0} controls the sparsity of \eqn{\beta},
#' and \eqn{\lambda_2 > 0} encourages smoothness between adjacent coefficients.
#' The function fits Fused Lasso and uses `easyGParetoptim` for Pareto optimization.
#' Each point in the plot corresponds to a model, where:
#' \itemize{
#'   \item x-axis: Number of non-zero coefficients (sparsity).
#'   \item y-axis: Residual sum of squares (RSS).
#'   \item z-axis: Roughness, i.e., mean absolute difference of consecutive beta coefficients.
#' }
#' The lower these three quantities, the better the model according to their respective criteria.
#'
#' @examples
#' \dontrun{
#' # Example data
#' set.seed(123)
#' n <- 100
#' p <- 10
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' beta_true <- c(0, 0, 1, 1, 1, 0, 0, 0, 2, 2)
#' y <- X %*% beta_true + rnorm(n)
#'
#' # Run the function and plot
#' result <- pared_FLasso(X, y, Pareto_budget = 80, plot_marker_symbol = 'diamond', plot_marker_size = 7)
#' result$figure
#' result$summary_table
#' }
#'
#' @importFrom plotly plot_ly add_trace layout
#' @importFrom dplyr %>%
#' @importFrom RMPSH RMPSolveH
#' @importFrom htmlwidgets saveWidget
#' @importFrom psych tr
#' @importFrom GPareto easyGParetoptim
#' @export

pared_FLasso <- function(X, y, lb = c(-3, -3), ub = c(1, 1),
                                   Pareto_budget = 100, plot_title = "Fused LASSO",
                                   plot_marker_color = 'rgba(255, 0, 0, 0.7)', plot_marker_border_color = 'black',
                                   plot_marker_size = 10, plot_marker_border_width = 2, plot_marker_symbol = 'circle',
                                   fontsize_x = 20, fontsize_y = 20, fontsize_z = 20, plot_title_size = 40,
                                   title_pos_x = 0.5, title_pos_y = 0.95, draw_projection = 1,
                                   proj_line_color = 'black', proj_line_width = 2, proj_marker_size = 2,
                                   proj_marker_color = 'rgba(255, 0, 0, 0.7)', proj_marker_border_width = 2,
                                   proj_marker_border_color = 'black') {


  ##############################################################################
  ## Required internal functions ###############################################
  ##############################################################################
  Moop_FusedLasso <- function(log_lambdas, X, y) {
    if (is.null(dim(log_lambdas))) {
      log_lambdas <- matrix(log_lambdas, nrow = 1)
    }

    log_lam1 <- log_lambdas[1]
    log_lam2 <- log_lambdas[2]

    lambda1 <- 10 ^ log_lam1
    lambda2 <- 10 ^ log_lam2
    n <- nrow(X)
    p <- ncol(X)


    g <- function(beta) FLasso_objective(X, y, beta, log_lam1, log_lam2)
    beta <- RMPSolveH(rep(0, p), g, rep(-10 ^ 3, p), rep(10 ^ 3, p), no_runs = 1, max_time = 10, print = 0)

    # Number of non-zero betas (absolute value greater than 10^(-4))
    non_zero_betas <- sum(abs(beta) > 10 ^ (- 3))

    # Residual sum of squares (RSS)
    residuals <- y - X %*% beta
    rss <- sum(residuals^2) / length(y)

    # Roughness term: sum of absolute differences between adjacent coefficients
    roughness <- sum(abs(diff(beta)))

    mult.obj <- cbind(non_zero_betas, rss, roughness)
    dimnames(mult.obj) <- NULL
    return(mult.obj)
  }

  ##############################################################################

  Moop_FusedLasso_wrapper <- function(log_lambdas) Moop_FusedLasso(log_lambdas, X, y)

  ##############################################################################
  ##############################################################################
  ##############################################################################

  # Run GP-based Pareto optimization with appropriate wrapper
  FusedLasso.Pareto.Optimals <- easyGParetoptim(fn = Moop_FusedLasso_wrapper, budget = Pareto_budget,
                                                lower = lb, upper = ub)

  # Extract parameter values and objective values
  ParetoOpt.params <- FusedLasso.Pareto.Optimals$par
  ParetoOpt.fitVal.nzVal.raw <- FusedLasso.Pareto.Optimals$value

  # Keep only unique solutions (non-duplicate rows)
  unique.rows <- !duplicated(ParetoOpt.fitVal.nzVal.raw)
  unique.row.numbers <- which(unique.rows)
  ParetoOpt.fitVal.nzVal.unique <- ParetoOpt.fitVal.nzVal.raw[unique.row.numbers, ]
  unique.opt.params <- ParetoOpt.params[unique.row.numbers, ]
  Solution <- ParetoOpt.fitVal.nzVal.unique

  # Summarize and format results
  Summary.Table.temp <- cbind(round(10 ^ unique.opt.params, 3), round(Solution, 3))
  Summary.Table <- Summary.Table.temp
  colnames(Summary.Table) <- c("lambda1", "lambda2", "Num. non-zero coeffs.", "RSS", "Roughness")

  # Prepare plot data
  plot_data <- data.frame(
    x = as.array(Solution[, 1]),
    y = as.array(round(Solution[, 2], 3)),
    z = as.array(round(Solution[, 3], 3)),
    lam1 = round(10 ^ unique.opt.params[, 1], 3),
    lam2 = round(10 ^ unique.opt.params[, 2], 3)
  )


  fig <- plot_ly(plot_data, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
                 name = 'Pareto-optimal',
                 marker = list(
                   size = plot_marker_size,
                   color = plot_marker_color,
                   line = list(
                     color = plot_marker_border_color,
                     width = plot_marker_border_width
                   ),
                   symbol = plot_marker_symbol
                 ),
                 text = ~paste("λ₁: ", t(lam1), "<br>","λ₂: ", t(lam2)),
                 hovertemplate = paste(
                   "No. of non-zero coeffs.: %{x}<br>",
                   "Residual sum of squares: %{y}<br>",
                   "Roughness: %{z}<br>",
                   "%{text}<extra></extra>"
                 ))%>%
    layout(scene = list(
      xaxis =list(title = list(text = "No. of non-zero coeffs.", font = list(size = fontsize_x))),
      yaxis =list(title = list(text = "RSS", font = list(size = fontsize_y))),
      zaxis =list(title = list(text = "Roughness", font = list(size = fontsize_z)))
    ),
    title = list(
      text = plot_title,
      font = list(size = plot_title_size),     # Plot title
      x = title_pos_x,                         # Centering the title horizontally
      y = title_pos_y,                         # Position slightly below the top edge
      xanchor = "center", yanchor = "top"
    )
    )

  # Optional: Add projections to XY plane
  if(draw_projection == 1) {
    # Add projection lines to the XY-plane (z = 0)
    for (i in 1:length(plot_data$x)) {
      fig <- fig %>% add_trace(
        x = c(plot_data$x[i], plot_data$x[i]),
        y = c(plot_data$y[i], plot_data$y[i]),
        z = c(min(plot_data$z), plot_data$z[i]),
        text = paste("λ₁: ", plot_data$lam1[i], "<br>","λ₂: ", plot_data$lam2[i]),
        type = 'scatter3d',
        mode = 'lines + markers',
        line = list(color = proj_line_color, width = proj_line_width, dash = 'dash'),
        marker = list(
          size = proj_marker_size,
          color = proj_marker_color,
          line = list(
            color = proj_marker_border_color,
            width = proj_marker_border_width
          )
        )
      )
    }
  }
  return(list(figure = fig, summary_table = Summary.Table))
}


#' Fused Lasso Objective Function value
#'
#' Computes the objective function value for the Fused Lasso model.
#'
#' @param X A numeric matrix of predictor variables (n x p), where n is the number
#' of samples and p is the number of predictors.
#' @param y A numeric vector of length n representing the response variable.
#' @param beta A numeric vector of length p representing the regression coefficients.
#' @param log_lam1 A numeric value representing the logarithm (log10-scale) of the regularization
#' parameter \eqn{\lambda_1} for the Lasso penalty (L1 regularization).
#' @param log_lam2 A numeric value representing the logarithm (log10-scale) of the regularization
#' parameter \eqn{\lambda_2} for the Fused Lasso penalty.
#'
#' @return A numeric value representing the total objective function value.
#'
#' @details The function calculates the total objective function for the Fused Lasso model.
#' The Fused Lasso objective function used is:
#' \deqn{
#' \mathcal{L}(\beta) = \frac{1}{2n} \|y - X\beta\|_2^2 + \lambda_1 \|\beta\|_1 + \lambda_2 \sum_{j=2}^{p} |\beta_j - \beta_{j-1}|
#' }
#' where \eqn{X} is the design matrix, \eqn{y} is the response vector,
#' \eqn{\beta} is the coefficient vector, \eqn{\lambda_1 > 0} controls the sparsity of \eqn{\beta},
#' and \eqn{\lambda_2 > 0} encourages smoothness between adjacent coefficients.
#'
#' @examples
#' \dontrun{
#' # Example data
#' set.seed(123)
#' X <- matrix(rnorm(100), 20, 5)
#' y <- rnorm(20)
#' beta <- rep(0, 5)
#' log_lam1 <- -2
#' log_lam2 <- -3
#'
#' # Compute the fused lasso objective function value
#' objective_value <- FLasso_objective(X, y, beta, log_lam1, log_lam2)
#' print(objective_value)
#' }
#'
#' @export


FLasso_objective <- function(X, y, beta, log_lam1, log_lam2) {
  lambda1 <- 10 ^ log_lam1
  lambda2 <- 10 ^ log_lam2
  # Residual sum of squares (RSS)
  n <- length(y)
  residuals <- y - X %*% beta
  rss <- sum(residuals ^ 2) / (2 * n)

  # L1 penalty on individual coefficients (lasso)
  l1_penalty <- lambda1 * sum(abs(beta))

  # Fused lasso penalty: sum of absolute differences between adjacent coefficients
  fused_penalty <- lambda2 * sum(abs(diff(beta)))  # L1 penalty on differences

  # Total objective (RSS + L1 penalty + fused lasso penalty)
  return(rss + l1_penalty + fused_penalty)
}
