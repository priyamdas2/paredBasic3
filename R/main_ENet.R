#' GP-Based Pareto Front for Elastic-Net
#'
#' This function performs a Pareto optimality search on the Elastic-Net problem
#' using GP-based optimization (`easyGParetoptim`) and creates a 3D scatter plot
#' of the solution space. The plot displays the number of non-zero coefficients,
#' L2 norm of coefficients, and deviance for different Pareto-optimal combinations
#' of the regularization parameters \eqn{\alpha} (0 < \eqn{\alpha} < 1) and \eqn{\lambda}.
#'
#' @param X A numeric matrix of predictor variables (n x p), where n is the number
#' of samples and p is the number of predictors.
#' @param y A numeric vector of length n representing the response variable.
#' @param lb A numeric vector specifying the lower bounds for the parameters: α (`alpha`) and log₁₀(λ) (`lambda`). Default is c(0, -6).
#' @param ub A numeric vector specifying the upper bounds for the parameters: α (`alpha`) and log₁₀(λ) (`lambda`). Default is c(1, 6).
#' @param Pareto_budget Number of evaluations in the GP-based Pareto optimization. Default is 100.
#' @param plot_title A character string for the title of the plot. Default is "Elastic-Net".
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
#' @return A list with two components:
#' \itemize{
#'   \item \code{figure}: A \code{plotly} 3D scatter plot object showing the Pareto front.
#'   \item \code{summary_table}: A data frame containing the summary of the Pareto
#'   optimal solutions, including the values for \eqn{\alpha}, \eqn{\log10} of \eqn{\lambda},
#'   number of non-zero coefficients, L2 norm, and deviance.
#' }
#'
#' @details
#' The Elastic Net objective function optimized in this process is:
#' \deqn{
#' \mathcal{L}(\beta) = \frac{1}{2n} \|y - X\beta\|_2^2 + \lambda \left[ \alpha \|\beta\|_1 + (1 - \alpha) \|\beta\|_2^2 \right]
#' }
#' where \eqn{X} is the design matrix, \eqn{y} is the response vector,
#' \eqn{\beta} is the coefficient vector, \eqn{\lambda > 0} is the regularization strength,
#' and \eqn{\alpha \in [0, 1]} balances Lasso and Ridge penalties.
#'
#' This function performs Gaussian Process (GP)-based Pareto optimization to find and visualize optimal Elastic Net models with respect to multiple objectives: sparsity (number of non-zero coefficients), L2 norm of coefficients, and deviance.
#'
#' The function uses `glmnet` for fitting Elastic Net models and `easyGParetoptim` for Pareto optimization.
#' Each point in the plot corresponds to a model, where:
#' \itemize{
#'   \item x-axis: Number of non-zero coefficients (sparsity).
#'   \item y-axis: L2 norm of the coefficients (shrinkage).
#'   \item z-axis: Deviance (model fit).
#' }
#' The lower these three quantities, the better the model according to their respective criteria.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' p <- 5
#' X <- matrix(rnorm(100), ncol = p)
#' n <- dim(X)[1]
#' beta.true <- matrix(c(1,2), ncol = 1)  # only first few coordinates are non-zero
#' y <- X[, 1:2] %*% beta.true + rnorm(n)
#'
#' A <- pared_ENet(X, y, Pareto_budget = 50)
#' A$fig
#' A$summary_table
#' }
#'
#' @importFrom plotly plot_ly add_trace layout
#' @importFrom dplyr %>%
#' @importFrom htmlwidgets saveWidget
#' @importFrom psych tr
#' @importFrom GPareto easyGParetoptim
#' @importFrom MASS mvrnorm
#' @importFrom glmnet glmnet
#' @export


pared_ENet <- function(X, y, lb = c(0, -6), ub = c(1, 6),
                                        Pareto_budget = 100, plot_title = "Elastic-Net",
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
  Moop_ElasticNet <- function(alpha.loglambda, X, y) {

    # Converting into 1x2 dimensional object

    if (is.null(dim(alpha.loglambda))) {
      alpha.loglambda <- matrix(alpha.loglambda, nrow = 1)
    }

    # Fitting elastic net

    ElasticNet.fit <- glmnet(X, y, alpha = alpha.loglambda[1],
                             lambda = 10 ^ (alpha.loglambda[2]))


    ## Extracting number of non-zero coefficients ##############################
    # Smaller g1 <=> less non-zero coeffs <=> better model

    g1 <- sum(abs(ElasticNet.fit$beta) > 1e-6)


    ## Extracting L2 norm of coefficients ######################################
    # Smaller g2 <=> less L2 norm of coeffs <=> better model

    g2 <- sqrt(sum(ElasticNet.fit$beta ^ 2))


    ## Extracting deviance #####################################################
    # Smaller dev <=> higher log-likelihood <=> better model

    Dev <- deviance(ElasticNet.fit)
    g3 <- Dev[1]

    if(length(Dev) > 1) {
      message("Deviance is of length more than 1")
    }
    mult.obj <- cbind(g1, g2, g3)
    dimnames(mult.obj) <- NULL
    return(mult.obj)
  }

  ##############################################################################

  Moop_ElasticNet_wrapper <- function(alpha.loglambda) Moop_ElasticNet(alpha.loglambda, X, y)

  ##############################################################################
  ##############################################################################
  ##############################################################################


  # Run GP-based Pareto optimization with appropriate wrapper
  ElasticNet.Pareto.Optimals <- easyGParetoptim(fn = Moop_ElasticNet_wrapper, budget = Pareto_budget,
                                                lower = lb, upper = ub)

  # Extract parameter values and objective values
  ParetoOpt.params <- ElasticNet.Pareto.Optimals$par
  ParetoOpt.fitVal.nzVal.raw <- ElasticNet.Pareto.Optimals$value


  # Keep only unique solutions (non-duplicate rows)
  unique.rows <- !duplicated(ParetoOpt.fitVal.nzVal.raw)
  unique.row.numbers <- which(unique.rows)
  ParetoOpt.fitVal.nzVal.unique <- ParetoOpt.fitVal.nzVal.raw[unique.row.numbers, ]
  unique.opt.params <- ParetoOpt.params[unique.row.numbers, ]
  Solution <- ParetoOpt.fitVal.nzVal.unique

  # Summarize and format results
  Summary.Table.temp <- cbind(round(unique.opt.params, 3), round(Solution, 3))
  Summary.Table <- Summary.Table.temp
  colnames(Summary.Table) <- c("alpha", "log10_lambda", "Num. non-zero coeffs.", "L2 norm", "Deviance")

  # Prepare plot data
  plot_data <- data.frame(
    x = as.array(Solution[, 1]),
    y = as.array(round(Solution[, 2], 3)),
    z = as.array(round(Solution[, 3], 3)),
    alpha = round(unique.opt.params[, 1], 3),
    loglambda = round(unique.opt.params[, 2], 3)
  )
  ##############################################
  ############################################


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
                 text = ~paste("\u03B1: ", alpha, "<br>","log10(λ): ", loglambda),
                 hovertemplate = paste(
                   "No. of non-zero coeffs.: %{x}<br>",
                   "L2 (\u03B2) : %{y}<br>",
                   "Deviance: %{z}<br>",
                   "%{text}<extra></extra>"
                 ))%>%
    layout(scene = list(
      xaxis =list(title = list(text = "No. of non-zero coeffs.", font = list(size = fontsize_x))),
      yaxis =list(title = list(text = "L2(\u03B2)", font = list(size = fontsize_y))),
      zaxis =list(title = list(text = "Deviance", font = list(size = fontsize_z)))
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
        text = paste("\u03B1: ", plot_data$alpha[i], "<br>","log10(λ): ", plot_data$loglambda[i]),
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
