# List of cancers to work with: OV, UCEC, UCS
rm(list=ls())
setwd("U:/MOOP_Christine/Github_vignette/case study")
library(paredBasic3)
library(ggplot2)
library(reshape2)
library(magick)
################################################################################
### Reading dataset ############################################################
################################################################################

names_vector <- read.csv("NExUS data/cancer_serial_no.csv", header = FALSE)[[1]]
target_terms <- c("OV", "UCEC", "UCS")
match_indices <- which(names_vector %in% target_terms)
matched_names <- names_vector[match_indices]
matrix_list <- lapply(match_indices, function(i) {
  file_path <- file.path("NExUS data", paste0(i, ".csv"))
  as.matrix(read.csv(file_path, header = FALSE))
})
names(matrix_list) <- matched_names


OV <- matrix_list$OV
UCEC <- matrix_list$UCEC
UCS <- matrix_list$UCS

################################################################################
### Selecting variables ########################################################
################################################################################

### Full pathway information ###################################################
# load("U:/MOOP_Christine/Github_vignette/Case study/NExUS data/RPPA_12_pathway.rda")
name_list <- vector("list", 12)
array_names_FULL <- c("APOPTOSIS", "CELL CYCLE", "DNA DMG RSPNS", "EMT", 
                      "HORMONE RECPTR", "HORMONE SIG BRST", "PI3K/AKT",
                      "RAS/MAPK", "RTK", "TSC/mTOR", "BREAST REACTIVE", 
                      "CORE REACTIVE") 
Apoptosis <- c("BAK", "BAX", "BID", "BIM", "CASPASE7CLEAVEDD198", "BAD_pS112", 
               "BCL2", "BCLXL", "CIAP")
Cell_cycle <- c("CDK1", "CYCLINB1", "CYCLINE2", "P27_pT157", "P27_pT198", "PCNA",
                "FOXM1")
DNA_damage_response <- c("X53BP1", "ATM", "CHK1_pS345", "CHK2_pT68", "KU80", 
                         "MRE11", "P53", "RAD50", "RAD51", "XRCC1")
EMT <- c("FIBRONECTIN", "NCADHERIN", "COLLAGENVI", "CLAUDIN7", "ECADHERIN", 
         "BETACATENIN", "PAI1")
Hormone_receptor <- c("ERALPHA", "ERALPHA_pS118", "PR", "AR")
Hormone_signaling_Breast <- c("BCL2", "INPP4B", "GATA3")
PI3K_AKT <- c("P27_pT157", "P27_pT198", "INPP4B", "AKT_pS473", "AKT_pT308", 
              "GSK3ALPHABETA_pS21S9",
              "GSK3_pS9", "PRAS40_pT246", "TUBERIN_pT1462", "PTEN")
RAS_MAPK <- c("ARAF_pS299", "CJUN_pS73", "CRAF_pS338", "JNK_pT183Y185", 
              "MAPK_pT202Y204", "MEK1_pS217S221", "P38_pT180Y182", 
              "P90RSK_pT359S363", "YB1_pS102")
RTK <- c("EGFR_pY1068", "EGFR_pY1173", "HER2_pY1248", "HER3_pY1289", "SHC_pY317",
         "SRC_pY416", "SRC_pY527")
TSC_mTOR <- c("X4EBP1_pS65", "X4EBP1_pT37T46", "X4EBP1_pT70", "P70S6K_pT389", 
              "MTOR_pS2448", "S6_pS235S236", "S6_pS240S244", "RB_pS807S811")
Breast_reactive <- c("BETACATENIN", "CAVEOLIN1", "MYH11", "RAB11", "GAPDH", "RBM15")
Core_reactive <- c("CLAUDIN7", "ECADHERIN", "BETACATENIN", "CAVEOLIN1", "RBM15")
################################################################################

proteins_here <- unique(c(Breast_reactive, Cell_cycle, Hormone_receptor, Hormone_signaling_Breast))
selected_variables <- read.csv("NExUS data/selected_variables.csv", header = FALSE)[[1]]
match_indices <- match(proteins_here,selected_variables)


ALL_samples <- list()
ALL_samples[[1]] <- OV[, match_indices]
ALL_samples[[2]] <- UCEC[, match_indices]
ALL_samples[[3]] <- UCS[, match_indices]

C <- 3
p <- length(match_indices)
sample_sizes <- c(dim(OV)[1], dim(UCEC)[1], dim(UCS)[1])

################################################################################
### Fitting JGL ################################################################
################################################################################

library(JGL)
library(mvtnorm)
library(psych)



sample <- ALL_samples
total_lambda_divisions <- 20
start_lambda <- 0.01
end_lambda <- 1
lambda_ones <- exp(seq(from = log(start_lambda), to = log(end_lambda), length.out= total_lambda_divisions))
lambda_twos <- exp(seq(from = log(start_lambda), to = log(end_lambda), length.out = total_lambda_divisions))
AIC_all_possible <- matrix(10 ^ 10, total_lambda_divisions, total_lambda_divisions)
for(i in 1:total_lambda_divisions) {
  for(j in 1:total_lambda_divisions) {
    JGL_result <- JGL(sample, penalty="group", lambda1 = lambda_ones[i], lambda2 = lambda_twos[j])
    Precision_estimated <- JGL_result$theta
    if(length(Precision_estimated[[1]]) == p^2) {
      num_non_zero_edges <- rep(0,C)
      AIC_sum <- 0
      for(c in 1:C){
        num_non_zero_edges[c] <- (length(which(abs(Precision_estimated[[c]])>0.0001))-20)/2
        S_mat <- t(sample[[c]]) %*% sample[[c]]
        AIC_sum <- AIC_sum + sample_sizes[c] * tr(S_mat %*% Precision_estimated[[c]]) - 
        sample_sizes[c] * log(det(Precision_estimated[[c]])) + 2 * num_non_zero_edges[c]
        }
      AIC_all_possible[i,j] <- AIC_sum
    }
  }
}

min_index <- which(AIC_all_possible == min(AIC_all_possible), arr.ind = TRUE)
lambda_opt <- c(lambda_ones[min_index[1]], lambda_twos[min_index[2]])

JGL_result_final <- JGL(ALL_samples, penalty="group", lambda1 = lambda_opt[1], lambda2=lambda_opt[2])
Precision_estimated_array <- JGL_result_final$theta

(length(which(abs(Precision_estimated_array[[1]]) > 10 ^ -3)) - p) / 2  # Number of non-zeros in Prec. Mat. 1

(length(which(abs(Precision_estimated_array[[2]]) > 10 ^ -3)) - p) / 2  # Number of non-zeros in Prec. Mat. 2

(length(which(abs(Precision_estimated_array[[3]]) > 10 ^ -3)) - p) / 2  # Number of non-zeros in Prec. Mat. 3

################################################################################
### PLOT: JGL group ############################################################
################################################################################

tol <- 1e-3
cancer_names <- c("Ovarian Cancer (AIC)", 
                  "Uterine Corpus Endometrial Carcinoma (AIC)",
                  "Uterine Carcinosarcoma (AIC)")

for (k in 1:3) {
  mat <- Precision_estimated_array[[k]]
  rownames(mat) <- colnames(mat) <- proteins_here
  
  mat_melt <- melt(mat)
  colnames(mat_melt) <- c("Row", "Col", "Value")
  
  # Reverse column order
  mat_melt$Col <- factor(mat_melt$Col, levels = rev(proteins_here))
  
  mat_melt$Color <- ifelse(abs(mat_melt$Value) < tol, "Zero",
                           ifelse(mat_melt$Value > 0, "Positive", "Negative"))
  
  mat_melt$Label <- ifelse(abs(mat_melt$Value) > tol,
                           sprintf("%.2f", mat_melt$Value),
                           "")
  
  plot_here <- ggplot(mat_melt, aes(x = Col, y = Row, fill = Color)) +
    geom_tile(color = "grey90") +
    geom_text(aes(label = Label), size = 3) + #, fontface = "bold") +
    scale_fill_manual(values = c("Zero" = "white", "Positive" = "salmon", "Negative" = "skyblue")) +
    scale_x_discrete(position = "top") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    ggtitle(cancer_names[k])
  
  print(plot_here)
  ggsave(filename = paste0("precision_heatmap_", k, ".jpg"), plot = plot_here,
         width = 8, height = 7, dpi = 300)
}



# Read images
img1 <- image_read("precision_heatmap_1.jpg")
img2 <- image_read("precision_heatmap_2.jpg")
img3 <- image_read("precision_heatmap_3.jpg")


combined <- image_append(c(img1, img2, img3))
print(combined)
image_write(combined, path = "precision_heatmaps_combined.jpg", format = "jpg")

################################################################################
### PLOT: pared_JGL group ######################################################
################################################################################

result <- pared_JGL(sample_list = ALL_samples, method = "group", Pareto_budget = 50)
result$summary_table
result$figure




lambda_opt_pared <- c(0.101, 0.27)
JGL_result_pared <- JGL(ALL_samples, penalty="group", lambda1 = lambda_opt_pared[1], lambda2=lambda_opt_pared[2])
Precision_estimated_array_pared <- JGL_result_pared$theta

(length(which(abs(Precision_estimated_array_pared[[1]]) > 10 ^ -3)) - p) / 2  # Number of non-zeros in Prec. Mat. 1

(length(which(abs(Precision_estimated_array_pared[[2]]) > 10 ^ -3)) - p) / 2  # Number of non-zeros in Prec. Mat. 2

(length(which(abs(Precision_estimated_array_pared[[3]]) > 10 ^ -3)) - p) / 2  # Number of non-zeros in Prec. Mat. 3


tol <- 1e-3
cancer_names <- c("Ovarian Cancer (pared)", 
                  "Uterine Corpus Endometrial Carcinoma (pared)",
                  "Uterine Carcinosarcoma (pared)")

for (k in 1:3) {
  mat <- Precision_estimated_array_pared[[k]]
  rownames(mat) <- colnames(mat) <- proteins_here
  
  mat_melt <- melt(mat)
  colnames(mat_melt) <- c("Row", "Col", "Value")
  
  # Reverse column order
  mat_melt$Col <- factor(mat_melt$Col, levels = rev(proteins_here))
  
  mat_melt$Color <- ifelse(abs(mat_melt$Value) < tol, "Zero",
                           ifelse(mat_melt$Value > 0, "Positive", "Negative"))
  
  mat_melt$Label <- ifelse(abs(mat_melt$Value) > tol,
                           sprintf("%.2f", mat_melt$Value),
                           "")
  
  plot_here <- ggplot(mat_melt, aes(x = Col, y = Row, fill = Color)) +
    geom_tile(color = "grey90") +
    geom_text(aes(label = Label), size = 3) + #, fontface = "bold") +
    scale_fill_manual(values = c("Zero" = "white", "Positive" = "salmon", "Negative" = "skyblue")) +
    scale_x_discrete(position = "top") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    ggtitle(cancer_names[k])
  
  print(plot_here)
  ggsave(filename = paste0("precision_heatmap_pared_", k, ".jpg"), plot = plot_here,
         width = 8, height = 7, dpi = 300)
}

# Read images
img1 <- image_read("precision_heatmap_pared_1.jpg")
img2 <- image_read("precision_heatmap_pared_2.jpg")
img3 <- image_read("precision_heatmap_pared_3.jpg")


combined <- image_append(c(img1, img2, img3))
print(combined)
image_write(combined, path = "precision_heatmaps_pared_combined.jpg", format = "jpg")