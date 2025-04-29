#' Generate Synthetic Multivariate Gaussian Samples with Evolving Graph Structures
#'
#' This function generates synthetic data from 4 (fixed) multivariate Gaussian distributions with 
#' evolving precision matrices (graph structures) of dimension 20 x 20. Starting from a base AR(2) graph, edges 
#' are progressively added and removed to simulate structural changes across four groups.
#'
#' @param sample_sizes A numeric vector of length 4 (fixed) specifying the number of observations 
#'   to generate for each of the four groups. Default is \code{c(20, 40, 60, 80)}.
#' @param rand_seed An integer specifying the random seed for reproducibility. Default is \code{1}.
#'
#' @return A list of length 4, each element being a matrix of synthetic samples (observations in rows)
#'   generated from a multivariate normal distribution with distinct precision matrices.
#'
#' @details 
#' The function constructs four precision matrices (A1 to A4), where A1 is a banded Toeplitz matrix
#' resembling an AR(2) structure. The subsequent matrices, A2, A3, and A4, are created by perturbing A1 
#' through random edge additions and deletions. The resulting matrices are adjusted to ensure 
#' positive definiteness, and samples are drawn using \code{rmvnorm()} 
#' from the \pkg{mvtnorm} package. 
#'
#'
#' @importFrom mvtnorm rmvnorm
#' @export

generate_sample <- function(sample_sizes = c(20, 40, 60, 80), rand_seed = 1) {
  set.seed(rand_seed)
  # Parameters (fixed, don't change)
  K <- 4
  p <- 20
  
  # A1 is an AR(2) graph
  A1 <- matrix(0, p, p)
  A1[1, ] <- c(1, 0.5, 0.4, rep(0, p - 3))
  A1 <- toeplitz(c(1, 0.5, 0.4, rep(0, p - 3)))
  
  # Number of edges
  n_edges <- (sum(A1 != 0) - p) / 2
  n_possible <- p * (p - 1) / 2
  
  # Locations of all nonzero entries of A1 above the diagonal
  rposA1 <- which(upper.tri(A1) & A1 != 0, arr.ind = TRUE)[, 1]
  cposA1 <- which(upper.tri(A1) & A1 != 0, arr.ind = TRUE)[, 2]
  
  # Locations of all zero entries of A1 above the diagonal
  rzeroA1 <- which(upper.tri(A1) & A1 == 0, arr.ind = TRUE)[, 1]
  czeroA1 <- which(upper.tri(A1) & A1 == 0, arr.ind = TRUE)[, 2]
  
  # Sample locations without replacement to perturb
  pos_inds <- sample(length(rposA1), n_edges)
  zero_inds <- sample(length(rzeroA1), n_possible - n_edges)
  
  A2 <- A1
  
  # Add 5 new and remove 5 edges to get A2 from A1
  for (j in 1:5) {
    sign <- ifelse(runif(1) > 0.5, 1, -1)
    nonzero_val <- runif(1, 0.4, 0.6) * sign
    A2[rzeroA1[zero_inds[j]], czeroA1[zero_inds[j]]] <- nonzero_val
    A2[czeroA1[zero_inds[j]], rzeroA1[zero_inds[j]]] <- nonzero_val
    
    # Remove a current nonzero value
    A2[rposA1[pos_inds[j]], cposA1[pos_inds[j]]] <- 0
    A2[cposA1[pos_inds[j]], rposA1[pos_inds[j]]] <- 0
  }
  
  # Now add 10 and remove 10 more to get A3
  A3 <- A2
  for (j in 6:15) {
    sign <- ifelse(runif(1) > 0.5, 1, -1)
    nonzero_val <- runif(1, 0.4, 0.6) * sign
    A3[rzeroA1[zero_inds[j]], czeroA1[zero_inds[j]]] <- nonzero_val
    A3[czeroA1[zero_inds[j]], rzeroA1[zero_inds[j]]] <- nonzero_val
    
    # Remove a current nonzero value
    A3[rposA1[pos_inds[j]], cposA1[pos_inds[j]]] <- 0
    A3[cposA1[pos_inds[j]], rposA1[pos_inds[j]]] <- 0
  }
  
  # Now change remaining to get A4
  A4 <- A3
  for (j in 16:20) {
    sign <- ifelse(runif(1) > 0.5, 1, -1)
    nonzero_val <- runif(1, 0.4, 0.6) * sign
    A4[rzeroA1[zero_inds[j]], czeroA1[zero_inds[j]]] <- nonzero_val
    A4[czeroA1[zero_inds[j]], rzeroA1[zero_inds[j]]] <- nonzero_val
    
    # Remove a current nonzero value
    A4[rposA1[pos_inds[j]], cposA1[pos_inds[j]]] <- 0
    A4[cposA1[pos_inds[j]], rposA1[pos_inds[j]]] <- 0
  }
  
  # Adjust revised matrices to ensure positive definiteness 
  A2 <- fix_matrix(A2, 1)
  A3 <- fix_matrix(A3, 1)
  A4 <- fix_matrix(A4, 1)
  
  # Generate precision matrices and sample
  C <- 4
  Precisions_true <- array(1, c(p, p, C))
  
  Precisions_true[, , 1] <- A1
  Precisions_true[, , 2] <- A2
  Precisions_true[, , 3] <- A3
  Precisions_true[, , 4] <- A4
  
  # Generating sample_list 
  sample_list <- vector("list", C)
  S_mat <- vector("list", C)
  
  for (c in 1:C) {
    sample_list[[c]] <- matrix(0, sample_sizes[c], p)
    for (kk in 1:sample_sizes[c]) {
      sample_list[[c]] <- rmvnorm(sample_sizes[c], mean = rep(0, p), 
                                  sigma = solve(Precisions_true[, , c]))
    }
  }
  
  return(sample_list)
}


#' GP-Based Pareto Front for Joint Graphical Lasso
#'
#' Generates an interactive 3D Pareto front plot for the Joint Graphical Lasso (JGL) using
#' GP-based optimization (`easyGParetoptim`). Supports both group and fused penalties.
#' Additionally, it provides the list of Pareto-optimal values of the tuning parameters \eqn{\lambda_1} 
#' and \eqn{\lambda_2}.
#'
#' @param sample_list A list of class-specific data samples used within the JGL wrapper functions.
#' @param lb A numeric vector specifying the lower bounds (log10-scale) for \eqn{\lambda_1} and \eqn{\lambda_2}. Default is c(-2, -2).
#' @param ub A numeric vector specifying the upper bounds (log10-scale) for \eqn{\lambda_1} and \eqn{\lambda_2}. Default is c(0, 0).
#' @param Pareto_budget Number of evaluations in the GP-based Pareto optimization. Default is 100.
#' @param method Type of JGL penalty. Either \code{"group"} or \code{"fused"}. Default is \code{"group"}.
#' @param plot_title Optional custom title for the plot. Automatically assigned if NULL.
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
#' @return A list with two elements:
#' \describe{
#'   \item{figure}{A \code{plotly} 3D scatter plot object showing the Pareto front.}
#'   \item{summary_table}{A data frame with optimal \code{lambda1}, \code{lambda2},
#'    and objectives, namely, total number of edges (for both group and fused), shared number of edges (group)
#'    or MAD of Precision matrices (fused) computed pair-wise, and AIC (for both group and fused).}
#' }
#' 
#' @details 
#' The function fits Joint Graphical Lasso using `JGL` and uses `easyGParetoptim` for Pareto optimization. 
#' Each point in the plot corresponds to a model, where:
#' \itemize{
#'   \item x-axis: Total number of edges (for both group and fused).
#'   \item y-axis: Shared number of edges (group)/Mean absolute deviation of Precision matrices computed pair-wise (fused).
#'   \item z-axis: Akaike Information Criteria (AIC; for both group and fused).
#' }
#' Lower values of all quantities but Shared number of edges (for that higher value is favored),
#' are favored while computing the Pareto-optimal points. For details on the objective functions of the Joint 
#' Graphical Lasso with 'group' and 'fused' penalties, see Danaher et al. (2014), linked below.
#'
#' @examples
#' \dontrun{
#' sample_data <- generate_sample(sample_sizes = c(30, 50, 40, 70), rand_seed = 123)
#' result <- pared_JGL(sample_list = sample_data, method = "group", Pareto_budget = 50)
#' result$figure
#' result$summary_table
#' 
#' sample_data2 <- generate_sample()
#' resultFused <- pared_JGL(sample_list = sample_data2, method = "fused", 
#' plot_marker_symbol = 'diamond', plot_marker_color = 'blue', 
#' plot_marker_size = 7, Pareto_budget = 40)
#' 
#' resultFused$figure
#' resultFused$summary_table
#' }
#'
#' @references
#' \itemize{
#'
#'   \item Danaher P, Wang P, Witten DM. \emph{The joint graphical lasso for inverse covariance estimation across multiple classes.}\cr
#'    J R Stat Soc Series B Stat Methodol. 2014 Mar;76(2):373-397. Doi: 10.1111/rssb.12033. \cr
#'          (available at \url{https://pubmed.ncbi.nlm.nih.gov/24817823/}).
#' }
#' 
#' @importFrom GPareto easyGParetoptim
#' @importFrom mvtnorm rmvnorm
#' @importFrom plotly plot_ly add_trace layout
#' @importFrom htmlwidgets saveWidget
#' @importFrom psych tr
#' @importFrom JGL JGL
#' @importFrom dplyr %>%
#' @export

pared_JGL <- function(sample_list, lb = c(-2, -2), ub = c(0, 0),
                            Pareto_budget = 100, method = "group", plot_title = NULL,
                            plot_marker_color = 'rgba(255, 0, 0, 0.7)', plot_marker_border_color = 'black',
                            plot_marker_size = 10, plot_marker_border_width = 2, plot_marker_symbol = 'circle',
                            fontsize_x = 20, fontsize_y = 20, fontsize_z = 20, plot_title_size = 40,
                            title_pos_x = 0.5, title_pos_y = 0.95, draw_projection = 1,
                            proj_line_color = 'black', proj_line_width = 2, proj_marker_size = 2,
                            proj_marker_color = 'rgba(255, 0, 0, 0.7)', proj_marker_border_width = 2,
                            proj_marker_border_color = 'black') {
  
  # Assign default plot title based on method
  if (is.null(plot_title)) {
    if (method == "group") {
      plot_title <- "JGL: Group LASSO"
    } else if (method == "fused") {
      plot_title <- "JGL: Fused LASSO"
    } 
  }
  C <- length(sample_list)
  
  ##############################################################################
  ## Required internal functions ###############################################
  ##############################################################################
  
  Moop_JGL_group <- function(log_lambdas, sample_list) {             
    # JGL has two tuning parameters: 
    # (i) lambda.1
    # (ii) lambda.2
    # Our MOOP function for 'Moop_JGL' takes input a given set of (scalar) 
    # values for lambda.1 and lambda.2, and it calculates the 
    # (i) AIC 
    # (ii) number of non-zero elements in precision matrix.
    # LL (for 1 matrix) = (n/2) * [ tr(X'X * W) - log(det(W))] + constant
    # (iii) AIC = 2K - 2LL
    
    
    C <- length(sample_list)
    p <- dim(sample_list[[1]])[2]
    
    # Converting into 1x2 dimensional object
    
    if (is.null(dim(log_lambdas))) {
      log_lambdas <- matrix(log_lambdas, nrow = 1)
    }
    
    # Fitting JGL
    
    JGL.fit <- JGL(sample_list, penalty = "group", 
                   lambda1 = 10 ^ log_lambdas[1],
                   lambda2 = 10 ^ log_lambdas[2])
    
    
    
    Precision.est <- JGL.fit$theta
    num.nodes.now <- sqrt(length(Precision.est[[1]]))
    
    
    if(num.nodes.now  < p) {
      cat(sprintf(">> Warning: At (lambda_1, lambda_2) = (%.3f,%.3f) the number of estimated nodes is observed to be less than p.\n", 10 ^ log_lambdas[1], 10 ^ log_lambdas[2]));
      g1 <- 10 ^ (10)
      g2 <- 10 ^ (10)
      g2.alt <- 10 ^ (10)
      g3 <- 10 ^ (10)
    } else {
      
      num.non.zero.edges <- rep(0, C)
      num.common.edges <- rep(0,C)
      sample.sizes <- rep(0, C)
      AIC.sum <- 0
      non.zero.edge.serials <- list()
      for(c in 1:C) {
        sample.sizes[c] <- dim(sample_list[[c]])[1]
        temp <- which(abs(Precision.est[[c]]) > 0.0001)
        non.zero.edge.serials[[c]] <- temp
        num.non.zero.edges[c] <- (length(temp) - p) / 2
        S.mat <- t(sample_list[[c]]) %*% sample_list[[c]]
        LL <-  - (sample.sizes[c] / 2) * (log(det(Precision.est[[c]])) - 
                                            tr(S.mat %*% Precision.est[[c]]))
        AIC <- 2 * num.non.zero.edges[c] - 2 * LL
        
        AIC.sum <- AIC.sum + AIC
        
        if (c == 1) {
          common.elements <- non.zero.edge.serials[[c]]
        } else {
          common.elements <- intersect(common.elements, non.zero.edge.serials[[c]])
        }
      }
      
      num.common.edges <- (length(common.elements) - p) / 2
      
      
      ## Extracting number of non-zero edges 
      # Smaller g2 <=> less non-zero coeffs <=> better model
      
      g1 <- sum(num.non.zero.edges) 
      
      ## Extracting number of common edges 
      # more common edges <=> better model (for JGL GROUP) 
      
      g2 <- -num.common.edges
      
      
      ## Extracting AIC 
      # Lower AIC <=> higher log-likelihood <=> better model 
      
      g3 <- AIC.sum
    } 
    
    mult.obj <- cbind(g1, g2, g3)
    
    dimnames(mult.obj) <- NULL
    return(mult.obj)  
  }
  
  
  
  
  Moop_JGL_fused <- function(log_lambdas, sample_list) {             
    # JGL has two tuning parameters: 
    # (i) lambda.1
    # (ii) lambda.2
    # Our MOOP function for 'Moop_JGL' takes input a given set of (scalar) 
    # values for lambda.1 and lambda.2, and it calculates the 
    # (i) AIC 
    # (ii) Precision MAD
    # LL (for 1 matrix) = (n/2) * [ tr(X'X * W) - log(det(W))] + constant
    # (iii) AIC = 2K - 2LL
    
    
    C <- length(sample_list)
    p <- dim(sample_list[[1]])[2]
    
    # Converting into 1x2 dimensional object
    
    if (is.null(dim(log_lambdas))) {
      log_lambdas <- matrix(log_lambdas, nrow = 1)
    }
    
    # Fitting JGL
    
    
    JGL.fit <- JGL(sample_list, penalty = "fused",
                   lambda1 = 10 ^ log_lambdas[1],
                   lambda2 = 10 ^ log_lambdas[2])
    
    
    
    Precision.est <- JGL.fit$theta
    num.nodes.now <- sqrt(length(Precision.est[[1]]))
    
    
    if(num.nodes.now  < p) {
      cat(sprintf(">> Warning: At (lambda_1, lambda_2) = (%.3f,%.3f) the number of estimated nodes is observed to be less than p.\n", 10 ^ log_lambdas[1], 10 ^ log_lambdas[2]));
      g1 <- 10^(10)
      g2 <- 10^(10)
      g2.alt <- 10^(10)
      g3 <- 10^(10)
    } else {
      
      num.non.zero.edges <- rep(0, C)
      num.common.edges <- rep(0,C)
      sample.sizes <- rep(0, C)
      AIC.sum <- 0
      non.zero.edge.serials <- list()
      for(c in 1:C) {
        sample.sizes[c] <- dim(sample_list[[c]])[1]
        temp <- which(abs(Precision.est[[c]]) > 0.0001)
        non.zero.edge.serials[[c]] <- temp
        num.non.zero.edges[c] <- (length(temp) - p) / 2
        S.mat <- t(sample_list[[c]]) %*% sample_list[[c]]
        LL <-  - (sample.sizes[c] / 2) * (log(det(Precision.est[[c]])) - 
                                            tr(S.mat %*% Precision.est[[c]]))
        AIC <- 2 * num.non.zero.edges[c] - 2 * LL
        
        AIC.sum <- AIC.sum + AIC
        
        if (c == 1) {
          common.elements <- non.zero.edge.serials[[c]]
        } else {
          common.elements <- intersect(common.elements, non.zero.edge.serials[[c]])
        }
      }
      
      num.common.edges <- (length(common.elements) - p) / 2
      
      
      ## Extracting number of non-zero edges 
      # Smaller g2 <=> less non-zero coeffs <=> better model
      
      g1 <- sum(num.non.zero.edges) 
      
      
      ## Extracting L1 norm of pairwise difference of precision matrices
      # less difference <=> better model (for JGL FUSED) 
      
      g2.alt <- mean_abs_distance_between_matrices(Precision.est)
      
      ## Extracting AIC 
      # Lower AIC <=> higher log-likelihood <=> better model 
      
      g3 <- AIC.sum
      
      
    } 
    
    mult.obj <- cbind(g1, g2.alt, g3)
    dimnames(mult.obj) <- NULL
    return(mult.obj)  
  }
  
  ##############################################################################
  
  Moop_JGL_group_wrapper <- function(log_lambdas) Moop_JGL_group(log_lambdas, sample_list)
  
  ##############################################################################
  
  Moop_JGL_fused_wrapper <- function(log_lambdas) Moop_JGL_fused(log_lambdas, sample_list)
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  
  
  # Run GP-based Pareto optimization with appropriate wrapper
  if(method == "group") {
    JGL.Pareto.Optimals <- easyGParetoptim(fn = Moop_JGL_group_wrapper, budget = Pareto_budget, 
                                           lower = lb, upper = ub)
  } else if (method == "fused") {
    JGL.Pareto.Optimals <- easyGParetoptim(fn = Moop_JGL_fused_wrapper, budget = Pareto_budget, 
                                           lower = lb, upper = ub)
  }
  
  
  # Extract parameter values and objective values
  ParetoOpt.lambda1.lambda2 <- JGL.Pareto.Optimals$par
  ParetoOpt.fitVal.nzVal.raw <- JGL.Pareto.Optimals$value
  
  # Keep only unique solutions (non-duplicate rows)
  unique.rows <- !duplicated(ParetoOpt.fitVal.nzVal.raw)
  unique.row.numbers <- which(unique.rows)
  ParetoOpt.fitVal.nzVal.unique <- ParetoOpt.fitVal.nzVal.raw[unique.row.numbers, ]
  unique.opt.params <- ParetoOpt.lambda1.lambda2[unique.row.numbers, ]
  Solution <- ParetoOpt.fitVal.nzVal.unique
  
  # Summarize and format results
  Summary.Table.temp <- cbind(round(unique.opt.params, 3), Solution)
  Summary.Table <- Summary.Table.temp
  Summary.Table[, 4] <- - Summary.Table.temp[, 4]
  
  if(method == "group") {
    colnames(Summary.Table) <- c("lambda1", "lambda2", "total edges", "shared edges", "AIC")
  } else if (method == "fused") {
    colnames(Summary.Table) <- c("lambda1", "lambda2", "total edges", "precision MAD", "AIC")
  }
  
  
  # Prepare plot data
  if(method == "group") {
    plot_data <- data.frame(
      x = Solution[, 1],
      y = -Solution[, 2],
      z = Solution[, 3],
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
                     "Total no. of edges: %{x}<br>",
                     "Shared no. of edges: %{y}<br>",
                     "AIC: %{z}<br>",
                     "%{text}<extra></extra>"
                   )
    )%>%
      layout(scene = list(
        xaxis =list(title = list(text = "Total no. of edges", font = list(size = fontsize_x))),
        yaxis =list(title = list(text = "Shared no. of edges", font = list(size = fontsize_y))),
        zaxis =list(title = list(text = "AIC", font = list(size = fontsize_z)))
      ),
      title = list(
        text = as.character(plot_title), 
        font = list(size = plot_title_size),
        x = title_pos_x,                         # Centering the title horizontally
        y = title_pos_y,                         # Position slightly below the top edge
        xanchor = "center", yanchor = "top"
      )
      )
  } else if (method == "fused") {
    plot_data <- data.frame(
      x = Solution[, 1],
      y = Solution[, 2],
      z = Solution[, 3],
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
                     "Total no. of edges: %{x}<br>",
                     "Precision MAD: %{y}<br>",
                     "AIC: %{z}<br>",
                     "%{text}<extra></extra>"
                   )
    )%>%
      layout(scene = list(
        xaxis =list(title = list(text = "Total no. of edges", font = list(size = fontsize_x))),
        yaxis =list(title = list(text = "Precision MAD", font = list(size = fontsize_y))),
        zaxis =list(title = list(text = "AIC", font = list(size = fontsize_z)))
      ),
      title = list(
        text = as.character(plot_title), 
        font = list(size = plot_title_size),
        x = title_pos_x,                         # Centering the title horizontally
        y = title_pos_y,                         # Position slightly below the top edge
        xanchor = "center", yanchor = "top"
      )
      )
  }
  
  
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



#' @keywords internal
mean_abs_distance_between_matrices <- function(matrix_list) {
  # Check that all elements are matrices of the same dimension
  dims <- sapply(matrix_list, dim)
  if (!all(apply(dims, 1, function(x) length(unique(x)) == 1))) {
    stop("All matrices must have the same dimensions.")
  }
  
  K <- length(matrix_list)
  if (K < 2) return(0)  # No pairs if less than 2 matrices
  
  # Generate all unique index pairs
  pair_indices <- combn(K, 2)
  
  # Compute total distances using vectorized apply
  total_dists <- vapply(
    seq_len(ncol(pair_indices)),
    function(k) {
      i <- pair_indices[1, k]
      j <- pair_indices[2, k]
      sum(abs(matrix_list[[i]] - matrix_list[[j]]))
    },
    numeric(1)
  )
  
  # Return the mean distance
  return(mean(total_dists))
}

#' @keywords internal
fix_matrix <- function(A, denom_factor) {
  # Ensure positive definiteness by normalizing off-diagonal elements
  p <- nrow(A)
  
  for (cur_row in 1:p) {
    cur_sum <- sum(abs(A[cur_row, ])) - 1
    if (cur_sum != 0) {
      A[cur_row, ] <- A[cur_row, ] / (denom_factor * cur_sum)
    }
    
    # Make sure diagonal entries are still 1
    A[cur_row, cur_row] <- 1
  }
  
  # Final matrix is the average of the matrix with its transpose
  A <- (A + t(A)) / 2
  
  return(A)
}