
## Single Cell Analysis of Natural Marine Phytoplankton Ionomes by LA-TOF-ICP-MS

# Isolate Multielement Data of Individual Phytoplankton Cells using Masks Segmented by a 
# U-Net Neural Network Model vs Manual Segmentation of Optical Light Microscope Imagery.

# SRM NIST61x series glasses (NIST-610 and NIST-612) used as External Standards and 
# Raw Counts were Normalized to P as an Internal Standard for Semi-Quantitative 
# LA-TOF-ICP-MS Multi-Elemental Data of Phytoplankton on a Single Cell level.

# Author: Joe Furby
# Date: 14/04/2025

# Set your working directory to Github repository folder
base_dir <- file.path(getwd(), "data")

library(abind)
library(agricolae) 
library(fields)
library(ggfortify) 
library(imager)
library(patchwork)
library(randomcoloR)
library(terra)
library(tidyverse)
library(vegan)

#### Read in functions ####

## Renumbers the ablation lines .csv file names in the chosen directory (e.g. line_1 --> line_001).
rename_laser_csv_files <- function(sample_dirs) {
  for (dir in sample_dirs) {
    if (dir.exists(dir)) {
      files <- list.files(path = dir, pattern = "\\.csv$", full.names = TRUE)  # Get all CSV files
      files_sorted <- files[order(as.numeric(sub(".*_(\\d+)\\.csv$", "\\1", files)))]
      count <- 1
      for (file in files_sorted) {
        file_name <- basename(file)
        num_match <- sub(".*_(\\d+)\\.csv$", "\\1", file_name)
        if (!is.na(as.numeric(num_match))) {
          new_name <- sprintf("line_%03d.csv", count)  # Format as "line_####.csv"
          count <- count + 1
        } else {
          new_name <- sprintf("line_%03d.csv", count)  # Default sequential numbering
          count <- count + 1
        }
        new_path <- file.path(dir, new_name)
        if (file != new_path) {
          file.rename(file, new_path)
        }
      }
    } else {
      message("Directory does not exist: ", dir)
    }
  }
}

## Reads all ablation line .csv files in a directory and processes them into a 3D array.
# Inputs: directory = directory of ablation line .csv files, isotopes = list of isotopes to process
# Outputs: return list (1:standards) of 3D arrays [1:shots, 1:isotopes, 1:line]
read_and_process_data <- function(directory, isotopes) {
  files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  data_list <- lapply(files, function(file) {
    data <- read.csv(file, skip = 11, header = TRUE)[,-c(1:3)]
    colnames(data) <- gsub("^X", "", colnames(data))
    matched_cols <- intersect(isotopes, colnames(data))
    data <- data[, matched_cols, drop = FALSE]
    if (nrow(data) == 0) return(NULL)
    as.matrix(data)
  })
  abind::abind(data_list, along = 3)
}

## Removes isotope-specific cps outliers using the interquartile range (IQR) method.
# Outliers are values outside the range [Q1 - 1.5*IQR, Q3 + 1.5*IQR]
remove_outliers_iqr <- function(data) {
  Q1 <- quantile(data, 0.25, na.rm = TRUE)
  Q3 <- quantile(data, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  data_clean <- data
  data_clean[data < lower_bound] <- NA
  data_clean[data > upper_bound] <- NA
  return(data_clean)
}

## Performs a model_I linear regression between two variables (X, Y) forced through (0, 0).
# Inputs: X = Standard Concentration Data [ug/g], Y = Standard Mean LA-ICP-ToFMS Data [cps].
# Outputs: Returns a list (intercept, slope, r, n, p, R2)
model_I_reg <- function(X, Y) {
  good <- !is.na(X) & !is.na(Y)
  X <- X[good]
  Y <- Y[good]
  fit <- lm(Y ~ X + 0)
  list(
    intercept = 0,  
    slope = coef(fit)[1],
    r = cor(X, Y),
    n = length(X),
    p = summary(fit)$coefficients[1, 4],  # Correct row for p-value
    R2 = summary(fit)$r.squared
  )
}

## Reads in Mask Images (.tiff) and Converts into a Matrix using the Terra Package (Hijmans, 2025).
# Input: mask_dirs = list of mask image file directories.
# Output: Returns a list of mask matrices (0 = background, 1:cell_n = labelled cell masks)
process_masks <- function(mask_dirs) {
  masks <- lapply(mask_dirs, function(dir) {
    mask <- rast(dir)  
    m <- matrix(values(mask), nrow = nrow(mask), ncol = ncol(mask), byrow = TRUE)
    m <- aperm(m, c(2, 1))
    m <- m[, c(ncol(m):1)]
    return(m)
  })
  return(masks)
}

## Function expands labeled masks to neighboring pixels specified by the kernel size in a raster or matrix. 
# Inputs: mask = SpatRast Mask Layer, kernel_size = Dilation Factor (odd integer, default = 3). 
# Output: A matrix with dilated regions using morphological dilation.
dilate_mask <- function(mask, kernel_size = 3) {
  mask_labels <- unique(values(mask))
  mask_labels <- mask_labels[mask_labels > 0]
  kernel_size <- max(3, ifelse(kernel_size %% 2 == 0, kernel_size + 1, kernel_size))
  struct_element <- matrix(1, nrow = kernel_size, ncol = kernel_size)
  dilated_mask <- mask
  for (n in mask_labels) {
    binary_mask <- ifel(mask == n, 1, 0)
    dilated_binary_mask <- focal(binary_mask, w = struct_element, fun = max, pad = TRUE, padValue = 0)
    dilated_mask <- ifel(dilated_binary_mask == 1, n, dilated_mask)
  }
  m <- matrix(values(dilated_mask), nrow = nrow(dilated_mask), ncol = ncol(dilated_mask), byrow = TRUE)
  m <- aperm(m, c(2, 1))
  m <- m[, ncol(m):1]
  return(m)
}

## Calculates model segmentation success metrics on a pixel basis against a check matrix
# Inputs: actual_mask = Ground truth binary mask, predicted_mask = Model segmented binary mask for the same region.
# Outputs: TP, TN, FP, FN, precision, recall, dsc (Dice Similarity Coefficient), 
#          m = manual segmentation, n = model segmentation, add = actual + predicted, subtract = actual - predicted
evaluation_metrics <- function(actual_mask, predicted_mask) {
  m <- ifelse(actual_mask > 0, 1, 0)
  n <- ifelse(predicted_mask > 0, 1, 0)
  add <- m + n
  subtract <- m - n
  TP <- sum(add == 2)
  TN <- sum(add == 0)
  FP <- sum(subtract == -1)
  FN <- sum(subtract == 1)
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  dsc <- (2 * precision * recall) / (precision + recall)
  return(list(TP = TP, TN = TN, FP = FP, FN = FN, 
              Precision = precision, Recall = recall, DSC = dsc, 
              m=m, n=n, add = add, subtract = subtract))
}


## Iterates through all samples segmentation matrices, calculates segmentation evaluation metrics, and saves the results.
# Inputs: check_masks = list of Ground truth segmentation matrices, model_masks = list of model segmentation matrices, 
#         dilated_model_masks = list of dilated model segmentation matrices.
# Outputs: model_results (and d_model_results) = model (and dilated model) segmentation evaluation metrics for the model predictions.
#          add_matrices (dadd_matrices) = list of actual + model matrices, subtract_matrices (dsubtract_matrices) = list of actual - model matrices.
evaluate_segmentation <- function(check_masks, model_masks, dilated_model_masks) {
  num_samples <- length(check_masks)
  col_names <- c("TP", "TN", "FP", "FN", "Precision", "Recall", "DSC")
  model_results <- matrix(nrow = num_samples, ncol = length(col_names))
  d_model_results <- matrix(nrow = num_samples, ncol = length(col_names))
  m_matrices <- list()
  n_matrices <- list()
  dn_matrices <- list()
  add_matrices <- list()
  dadd_matrices <- list()
  subtract_matrices <- list()
  dsubtract_matrices <- list()
  for (i in 1:num_samples) {
    model_metrics <- evaluation_metrics(check_masks[[i]], model_masks[[i]])
    d_model_metrics <- evaluation_metrics(check_masks[[i]], dilated_model_masks[[i]])
    model_results[i, ] <- c(model_metrics$TP, model_metrics$TN, model_metrics$FP, model_metrics$FN, 
                            model_metrics$Precision, model_metrics$Recall, model_metrics$DSC)
    d_model_results[i, ] <- c(d_model_metrics$TP, d_model_metrics$TN, d_model_metrics$FP, d_model_metrics$FN, 
                              d_model_metrics$Precision, d_model_metrics$Recall, d_model_metrics$DSC)
    m_matrices[[i]] <- model_metrics$m
    n_matrices[[i]] <- model_metrics$n
    dn_matrices[[i]] <- d_model_metrics$n
    add_matrices[[i]] <- model_metrics$add
    subtract_matrices[[i]] <- model_metrics$subtract
    dadd_matrices[[i]] <- d_model_metrics$add
    dsubtract_matrices[[i]] <- d_model_metrics$subtract
  }
  row_names <- paste0("S", 1:num_samples)
  rownames(model_results) <- row_names
  colnames(model_results) <- col_names
  rownames(d_model_results) <- row_names
  colnames(d_model_results) <- col_names
  return(list(model_results = model_results, d_model_results = d_model_results, 
              m_matrices = m_matrices, n_matrices = n_matrices, dn_matrices = dn_matrices,
              add_matrices = add_matrices, dadd_matrices = dadd_matrices, 
              subtract_matrices = subtract_matrices, dsubtract_matrices = dsubtract_matrices))
}

## Function performs one-way ANOVA and Tukey's HSD test for P-normalized Cellular Quotas across Sampling Stations.
# Inputs: variables_to_analyze = Vector of element:P to analyze, ionomes = Data frame containing Station and Ionome data.
# Output: A data frame with Tukey group classifications for each station per variable.
perform_anova_and_tukey_all <- function(variables_to_analyze, ionomes) {
  tukey_groups <- data.frame(Station = levels(ionomes$Station))
  for (var in variables_to_analyze) {
    formula <- as.formula(paste("`", var, "` ~ Station", sep = ""))
    amod <- aov(formula, data = ionomes)
    tukey_results <- HSD.test(amod, "Station")
    groups <- tukey_results$groups
    group_names <- rownames(groups)
    station_order <- levels(ionomes$Station)
    group_labels <- rep(NA, length(station_order))
    for (i in 1:length(station_order)) {
      station <- station_order[i]
      group_labels[i] <- groups[station == group_names, "groups"]
    }
    tukey_groups[[var]] <- group_labels
  }
  return(tukey_groups)
}

## Function creates a dynamic box plot with flexible Y-axis limits, breaks, and labels.
# Inputs: element = element:P to plot, data = log10 Cellular Ionomes with Sampling Station and Segmentation.
#         tukey_groups = Sampling Station Tukey group classifications,
#         y_limits = List of Y-axis limits for each element, y_breaks = List of Y-axis breaks for each element.
# Output: A ggplot2 box plot of log10 Cellular Ionomes Between Different Segmentation Approaches
plot_boxplot <- function(element, data, tukey_groups, y_limits, y_breaks) {
  label_positions <- tukey_groups %>%
    dplyr::select(Station, all_of(element)) %>%
    rename(Group = all_of(element))
  upper_limit <- y_limits[[element]][2]
  segmentation_levels <- unique(data$Segmentation)
  color_mapping <- setNames(c("grey", "red"), segmentation_levels)
  p <- ggplot(data, aes(x = Station, y = .data[[element]], fill = Segmentation)) +
    geom_boxplot(outlier.shape = NA) +  
    scale_fill_manual(values = color_mapping) + 
    geom_text(data = label_positions, aes(x = Station, y = upper_limit - 0.5, label = Group), 
              inherit.aes = FALSE, color = "black", size = 5, vjust = -0.5) +
    labs(title = paste(element), x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1), legend.position = "none") +
    scale_y_continuous(
      limits = y_limits[[element]], 
      breaks = y_breaks[[element]],
      labels = function(x) parse(text = paste0("10^", x))
    )
  return(p)
}

#### Calibration against SRM NIST61x glasses (NIST-610 and NIST-612) ####

## Enter Standard Data Directory and Select Standards to Process
std_dir <- file.path(base_dir, "LA TOF-ICP-MS", "Calibration")
standards <- c("NIST_610_scan", "NIST_612_scan")
std_dirs <- file.path(std_dir, standards)
std <- length(standards)
#rename_laser_csv_files(std_dirs)

## Choose which Isotopes to Extract from LA-TOF-ICP-MS Data 
isotopes <- c("23Na.cps","24Mg.cps","28Si.cps","31P.cps","55Mn.cps","56Fe.cps","63Cu.cps","64Zn.cps")
elements <- c("Na","Mg","Si","P","Mn","Fe","Cu","Zn")
ele <- length(elements)

# Read and process standard data
std_list <- lapply(std_dirs, function(dir) {
  std_data <- read_and_process_data(dir, isotopes)
  if (!is.null(std_data)) aperm(std_data, c(1, 3, 2)) else NULL
})
# Function to compute outlier limits and mean cps for each element across both standards
outlier_limits_and_mean <- lapply(std_list, function(std_data) {
  limit <- matrix(NA, nrow = 2, ncol = ele)
  mean_cps <- numeric(ele)
  for (elem in 1:ele) {
    cleaned_data <- remove_outliers_iqr(std_data[,,elem])
    limit[1, elem] <- min(cleaned_data, na.rm = TRUE)
    limit[2, elem] <- max(cleaned_data, na.rm = TRUE)
    mean_cps[elem] <- mean(cleaned_data, na.rm = TRUE)
  }
  colnames(limit) <- elements
  return(list(limit = limit, mean_cps = mean_cps))
})
outlier_limits <- lapply(outlier_limits_and_mean, function(x) x$limit)
mean_cps <- sapply(outlier_limits_and_mean, function(x) x$mean_cps)

# Plot the 2D element heat maps removing outliers
plot_limits <- matrix(NA, nrow = 2, ncol = ele)
for (elem in 1:ele) {
  lower_limits <- sapply(outlier_limits, function(x) x[1, elem])
  upper_limits <- sapply(outlier_limits, function(x) x[2, elem])
  plot_limits[1, elem] <- min(lower_limits, na.rm = TRUE)
  plot_limits[2, elem] <- max(upper_limits, na.rm = TRUE)
}
png(file.path(base_dir,"Output","NIST_heatmaps.png"), width = 900, height = 400, res = 150)
par(mfrow = c(2 * ele, 1), mar = c(0.1, 5, 0.1, 0.1))
for (i in seq_len(ele)) {
  for (s in 1:2) {
    image(std_list[[s]][,,i], zlim = plot_limits[, i], col = hcl.colors(20), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    if (s == 1) {
      mtext(elements[i], side = 3, line = -2, col = "black", adj = -0.12, cex = 1)
    }
    mtext(610 + 2 * (s - 1), side = 3, line = -1, col = "black", adj = -0.06, cex = 0.5)
  }
}
dev.off()

# Plot the NIST glass element cps distributions as histograms
png(file.path(base_dir,"Output","NIST_Histogram_Panel.png"), width = 2500, height = 1000)
par(mfrow = c(std, ele), mar = c(5, 5, 3, 2), cex = 1.4, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5, cex.sub = 1.4)
for (s in 1:std) {
  for (e in 1:ele) {
    clean_data <- remove_outliers_iqr(std_list[[s]][,,e])
    hist(
      clean_data, breaks = 20, col = "lightblue", border = "white",
      main = paste(elements[e], "-", standards[s]),
      xlab = paste(elements[e], "cps"), ylab = "Frequency",
      xlim = c(outlier_limits[[s]][1, e], outlier_limits[[s]][2, e])
    )
    abline(v = mean_cps[e, s], col = "red", lty = 2, lwd = 2)
    legend("topright", legend = paste("Mean =", round(mean_cps[e, s], 2)), bty = "n", cex = 1.2)
  }
}
dev.off()

# Load in the NIST glass element concentrations
concs <- read.csv(file.path(base_dir,"LA TOF-ICP-MS","Calibration","NIST Concentrations.csv"))[,-1]
colnames(concs) <- standards
rownames(concs) <- elements
colnames(mean_cps) <- standards
rownames(mean_cps) <- elements

# Perform a model-I regression (intensity[cps] = Factor * conc[ug/g] + 0) and plot the calibration curves
results <- data.frame(Element = elements, Slope = NA, R2 = NA)
plot_list <- list()
for (i in seq_along(elements)) {
  element <- elements[i]
  cps <- as.numeric(mean_cps[i, ])  
  conc <- as.numeric(concs[i, ])  
  reg <- model_I_reg(conc, cps)
  results[i, c("Slope", "R2")] <- c(reg$slope, reg$R2)
  data <- data.frame(
    Concentration = conc,
    CPS = cps
  )
  p <- ggplot(data, aes(x = Concentration, y = CPS)) +
    geom_point() +
    geom_abline(intercept = reg$intercept, slope = reg$slope, color = "red", linetype = "solid", linewidth = 1) +
    labs(
      title = paste(element, "Calibration"),
      x = "Concentration (µg/g)",
      y = "Mean CPS"
    ) +
    theme_bw() + 
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.spacing.y = unit(0.2, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 50, 10)
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(conc, na.rm = TRUE) * 1.1)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(cps, na.rm = TRUE) * 1.1)) +  
    annotate(
      "label", x = max(conc, na.rm = TRUE) * 0.1, y = max(cps, na.rm = TRUE) * 0.9,
      label = paste("NIST: R² =", round(reg$R2, 4), "\nSlope =", round(reg$slope, 4)),
      color = "black", fill = "white", hjust = 0, vjust = 1, size = 4, label.padding = unit(0.2, "lines")
    )
  plot_list[[i]] <- p
}
panel_plot <- wrap_plots(plotlist = plot_list, nrow = 2, ncol = 4)
ggsave(file.path(base_dir, "Output","Calibration_Plots.png"), plot = panel_plot, width = 16, height = 10)

# Save regression results (Factors and R^2)
write.csv(results, file.path(base_dir, "LA TOF-ICP-MS","Calibration","NIST Calibration Results.csv"), row.names = FALSE)


#### Loading in LA ICP-ToFMS element data, Cell Segmentation Masks, and Optical Images ####

# Enter Sample Data Directory and Select Samples to Process
sample_dir <- file.path(base_dir,"LA TOF-ICP-MS") 
samples <- c("STATION1","STATION2","STATION3","STATION4","STATION5","STATION6")
sample_dirs <- file.path(sample_dir, samples)
#rename_laser_csv_files(sample_dirs)
  
# Choose which Isotopes to Extract from LA ICP-ToFMS Data 
isotopes <- c("23Na.cps","24Mg.cps","28Si.cps","31P.cps","55Mn.cps","56Fe.cps","63Cu.cps","64Zn.cps")
elements <- c("Na","Mg","Si","P","Mn","Fe","Cu","Zn")
iso <- length(isotopes)
ele <- length(elements)
n   <- length(samples)

# Calibration Factors from a model-I regression of external SRM NIST61x glasses
factors <- read.csv(file.path(base_dir,"LA TOF-ICP-MS","Calibration","NIST Calibration Results.csv"))[,2]

# Define isotopic masses for elements of interest to convert concentration data to molar mass 
Mr <- c(23, 24, 28, 31, 55, 56, 63, 64)

## Load in and process the La ICP-ToFMS data
# Convert raw isotope intensity (cps) to concentration (NIST-ug/g) by dividing by calibration factor and atomic mass
# A list of 1:samples (3D arrays [1:shots per line, 1:ablation lines, 1:elements])
sample_list <- lapply(file.path(sample_dirs), function(dir) {
  element_data <- read_and_process_data(dir, isotopes)
  if (is.null(element_data)) return(NULL)
  element_data <- aperm(element_data, c(1, 3, 2)) 
  element_data <- element_data[, ncol(element_data):1, , drop = FALSE]
  for (i in 1:ele) {
    element_data[,,i] <- element_data[,,i] / (Mr[i]*factors[i])
  }
  cat("Directory:", dir, "Dimensions:", if (!is.null(element_data)) dim(element_data) else "NULL", "\n")
  return(element_data)
})

## Cell segmentation using the Aligned Optical Light Microscope Images
# Load in Optically Segmented Mask Layers Image Directories and Process into 2D Matrices (with and without dilation)
validation_dir <- file.path(base_dir, "U-Net Model Validation") 
mask_dir <- file.path(base_dir, "Segmentation Masks")
manual_mask_names <- c("S1_manual.tiff","S2_manual.tiff","S3_manual.tiff","S4_manual.tiff","S5_manual.tiff","S6_manual.tiff")
model_mask_names <- c("S1_model.tiff","S2_model.tiff","S3_model.tiff","S4_model.tiff","S5_model.tiff","S6_model.tiff")

manual_v_masks <- process_masks(file.path(validation_dir, manual_mask_names))
man_v_masks <- lapply(file.path(validation_dir, manual_mask_names), rast)
dilated_manual_v_masks <- lapply(man_v_masks, dilate_mask)
model_v_masks <- process_masks(file.path(validation_dir, model_mask_names))
mod_v_masks <- lapply(file.path(validation_dir, model_mask_names), rast)
dilated_model_v_masks <- lapply(mod_v_masks, dilate_mask)

manual_masks <- process_masks(file.path(mask_dir, manual_mask_names))
man_masks <- lapply(file.path(mask_dir, manual_mask_names), rast)
dilated_manual_masks <- lapply(man_masks, dilate_mask)
model_masks <- process_masks(file.path(mask_dir, model_mask_names))
mod_masks <- lapply(file.path(mask_dir, model_mask_names), rast)
dilated_model_masks <- lapply(mod_masks, dilate_mask)

# Identifying cellular material as pixels with  Mg, Si, and P, exceeding the 95th Percentile
get_threshold <- function(matrix_2d) {
  return(quantile(matrix_2d, probs = 0.95, na.rm = TRUE))
}
MgSiP_95 <- lapply(1:length(sample_list), function(station_index) {
  sample_data <- sample_list[[station_index]]
  Mg_map <- sample_data[,,2]
  Si_map <- sample_data[,,3]
  P_map  <- sample_data[,,4]
  Si_thresh <- get_threshold(Si_map)
  Mg_thresh <- get_threshold(Mg_map)
  P_thresh  <- get_threshold(P_map)
  e_map <- matrix(0, nrow = nrow(Si_map), ncol = ncol(Si_map))
  e_map[(Si_map > Si_thresh)] <- e_map[(Si_map > Si_thresh)] + 1
  e_map[(Mg_map > Mg_thresh)] <- e_map[(Mg_map > Mg_thresh)] + 1
  e_map[(P_map > P_thresh)] <- e_map[(P_map > P_thresh)] + 1
  e_map[e_map < 3] <- 0
  e_map[e_map >= 3] <- 1
  return(e_map)
})

## Load in Light Microscope Images to be used in figures
# Enter Image Directories and Image File Names
image_dir <- file.path(base_dir, "Optical Images")
image_names <- c("S1.jpg","S2.jpg","S3.jpg","S4.jpg","S5.jpg","S6.jpg")
images <- lapply(file.path(image_dir, image_names), function(dir) {
  i <- jpeg::readJPEG(dir)
  i
})
cell_names <- c("S1_cell.jpg","S2_cell.jpg","S3_cell.jpg","S4_cell.jpg","S5_cell.jpg","S6_cell.jpg")
cell_images <- lapply(file.path(image_dir, cell_names), function(dir) {
  i <- jpeg::readJPEG(dir)
  i
})

#### U-Net Semantic Segmentation Model Validation ####
# Segmentation Validation Indices range from 0 - worst to 1 - best
# Precision is a measure of how many predicted cell pixels segmented by the model were correct
# Recall is a measure of how many cell pixels were detected by the model
# Dice Similarity Coefficient (DSC) balances Precision and Recall as an overall performance metric

# Evaluate the U-net model segmentation against manually segmented images on a pixel basis (with and without dilation)
# Create a panel plot displaying the results of the validation for a region of one validation image
optical_validation <- evaluate_segmentation(manual_v_masks, model_v_masks, dilated_model_v_masks)
print(optical_validation$model_results)
print(optical_validation$d_model_results)

png(file.path(base_dir, "Output","S1 Optical Model Validation.png"), width = 4000, height = 2000, res = 500)
im <- jpeg::readJPEG(file.path(base_dir,"Optical Images","S1_validation.jpg"))
x_cor <- 190:462
y_cor <- 0:272
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 4, byrow = TRUE), widths = c(1, 1, 1, 1)) 
par(mar = c(1, 1, 1, 1))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))  
rasterImage(im, 0, 0, 1, 1)  
mtext("a", side = 3, line = 0, col = "black", adj = 0, cex = 1, font = 2)
image(x_cor, y_cor, optical_validation$add_matrices[[1]][x_cor, y_cor], zlim = c(0, 2), col = hcl.colors(3), xaxt = "n", yaxt = "n")
mtext("c", side = 3, line = 0, col = "black", adj = 0, cex = 1, font = 2)
image(x_cor, y_cor, optical_validation$dadd_matrices[[1]][x_cor, y_cor], zlim = c(0, 2), col = hcl.colors(3), xaxt = "n", yaxt = "n")
mtext("d", side = 3, line = 0, col = "black", adj = 0, cex = 1, font = 2)
image.plot(zlim = c(0, 2), col = hcl.colors(3), legend.only = TRUE, 
           breaks = c(0, 1, 2, 3),
           axis.args = list(at = c(0.5, 1.5, 2.5), labels = c("TN", "F", "TP")), 
           legend.width = 2, legend.mar = 0)
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
image(x_cor, y_cor, manual_v_masks[[1]][x_cor, y_cor], zlim = c(0, 1), 
      col = gray.colors(100, start = 0, end = 1), xaxt = "n", yaxt = "n")
mtext("b", side = 3, line = 0, col = "black", adj = 0, cex = 1, font = 2)
image(x_cor, y_cor, optical_validation$subtract_matrices[[1]][x_cor, y_cor], zlim = c(-1, 1), col = hcl.colors(3), xaxt = "n", yaxt = "n")
mtext("e", side = 3, line = 0, col = "black", adj = 0, cex = 1, font = 2)
image(x_cor, y_cor, optical_validation$dsubtract_matrices[[1]][x_cor, y_cor], zlim = c(-1, 1), col = hcl.colors(3), xaxt = "n", yaxt = "n")
mtext("f", side = 3, line = 0, col = "black", adj = 0, cex = 1, font = 2)
image.plot(zlim = c(-1, 1), col = hcl.colors(3), legend.only = TRUE, 
           breaks = c(-1, 0, 1, 2),
           axis.args = list(at = c(-0.5, 0.5, 1.5), labels = c("FP", "T", "FN")), 
           legend.width = 2, legend.mar = 0)
dev.off()

## Evaluate the U-net model segmentation against the MgSiP_95 element derived cell layer on a pixel basis (with and without dilation)
# Create a panel plot displaying the results of the validation for a region of one validation image
element_validation <- evaluate_segmentation(MgSiP_95, model_masks, dilated_model_masks)
print(element_validation$model_results)
print(element_validation$d_model_results)

png(file.path(base_dir,"Output","S2 Element Model Validation.png"), width = 4000, height = 2000, res = 500)
im <- jpeg::readJPEG(file.path(base_dir,"Optical Images","S2_cropped.jpg"))
x_cor <- 170:270
y_cor <- 0:100
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 4, byrow = TRUE), widths = c(1, 1, 1, 1)) 
par(mar = c(1, 1, 1, 1))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))  
rasterImage(im, 0, 0, 1, 1)  
mtext("a", side = 3, line = 0, col = "black", adj = 0, cex = 1, font = 2)
image(x_cor, y_cor, element_validation$add_matrices[[2]][x_cor, y_cor], zlim = c(0, 2), col = hcl.colors(3), xaxt = "n", yaxt = "n")
mtext("c", side = 3, line = 0, col = "black", adj = 0, cex = 1, font = 2)
image(x_cor, y_cor, element_validation$dadd_matrices[[2]][x_cor, y_cor], zlim = c(0, 2), col = hcl.colors(3), xaxt = "n", yaxt = "n")
mtext("d", side = 3, line = 0, col = "black", adj = 0, cex = 1, font = 2)
image.plot(zlim = c(0, 2), col = hcl.colors(3), legend.only = TRUE, 
           breaks = c(0, 1, 2, 3),
           axis.args = list(at = c(0.5, 1.5, 2.5), labels = c("TN", "F", "TP")), 
           legend.width = 2, legend.mar = 0)
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
image(x_cor, y_cor, MgSiP_95[[2]][x_cor, y_cor], zlim = c(0, 1), 
      col = gray.colors(100, start = 0, end = 1), xaxt = "n", yaxt = "n")
mtext("b", side = 3, line = 0, col = "black", adj = 0, cex = 1, font = 2)
image(x_cor, y_cor, element_validation$subtract_matrices[[2]][x_cor, y_cor], zlim = c(-1, 1), col = hcl.colors(3), xaxt = "n", yaxt = "n")
mtext("e", side = 3, line = 0, col = "black", adj = 0, cex = 1, font = 2)
image(x_cor, y_cor, element_validation$dsubtract_matrices[[2]][x_cor, y_cor], zlim = c(-1, 1), col = hcl.colors(3), xaxt = "n", yaxt = "n")
mtext("f", side = 3, line = 0, col = "black", adj = 0, cex = 1, font = 2)
image.plot(zlim = c(-1, 1), col = hcl.colors(3), legend.only = TRUE, 
           breaks = c(-1, 0, 1, 2),
           axis.args = list(at = c(-0.5, 0.5, 1.5), labels = c("FP", "T", "FN")), 
           legend.width = 2, legend.mar = 0)
dev.off()

#### Plotting 2D False Colour Heatmaps, Mask Maps, and Optical Images ####

# Plot Manually and Model Segmented Mask Locations (optional: With Mask Numbers)
mask_layers <- list(dilated_manual_masks, dilated_model_masks)
layer_names <- c("Dilated_Manual_Masks", "Dilated_Model_Masks")
for (layer in 1:length(mask_layers)) {
  masks <- mask_layers[[layer]]
  clean_name <- gsub("[^a-zA-Z0-9_-]", "_", layer_names[layer])  # replace unsafe chars with "_"
  png(filename = file.path(base_dir, "Output", paste0(clean_name, ".png")), 
      width = 1000, height = 1250, res = 150)  
  layout(matrix(1:6, ncol = 1))
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  for (i in 1:6) {
    num_mask <- max(masks[[i]], na.rm = TRUE)
    color_palette <- distinctColorPalette(num_mask)
    image(masks[[i]], zlim = c(1, num_mask), col = color_palette, 
          main = "", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    mtext(paste0("S", i), side = 3, line = -2.5, col = "black", adj = 0.02, cex = 1.5)
    ## Optional: Annotate each mask with its number
    # for (cell in 1:num_mask) {
    #   cell_pos <- which(masks[[i]] == cell, arr.ind = TRUE)[1, ]
    #   text(cell_pos[1] * (1/nrow(masks[[i]])), cell_pos[2] * (1/ncol(masks[[i]])), 
    #        labels = cell, col = "black", cex = 1, font = 2)
    # }
  }
  dev.off()
}

#Plot the Aligned Light Microscope Images (Arithmetic: A*A) for Each Station
png(file.path(base_dir, "Output","Aligned Sample Optical Images.png"), width = 1000, height = 1250, res = 150) 
layout(matrix(1:6, ncol = 1))
par(mar = c(0.1, 0.1, 0.1, 0.1))
for (i in 1:6) {
  img <- images[[i]]
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))  
  rasterImage(img, 0, 0, 1, 1) 
  mtext(paste0("S",i), side = 3, line = -2.5, col = "white", adj = 0.05, cex = 1.5)
}
dev.off()

# Create a Panel Plot focusing on a select few cells with Light Microscope Image, Segmented Masks, 
# and False-Colour 2D Element Concentration Heat Maps (5th to 95th percentile ranges)
limits <- matrix(NA, nrow = 2, ncol = 8)  
for (elem in 1:ele) {
  all_values <- c() 
  for (station in seq_along(sample_list)) {
    element_data <- sample_list[[station]][, , elem]  
    masked_values <- element_data[model_masks[[station]] != 0]
    all_values <- c(all_values, masked_values[!is.na(masked_values)])
  }  
  limits[, elem] <- quantile(all_values, probs = c(0.05, 0.95), na.rm = TRUE)
}
# Set X and Y coordinates and Mask Numbers of Chosen Cell(s) of Interest.
x_cor <- matrix(c(   # Define X coordinate ranges
  290, 199, 388, 415,  # Lower bounds
  322, 217, 400, 435   # Upper bounds
), nrow = 2, byrow = TRUE)
y_cor <- matrix(c(   # Define Y coordinate ranges
   62,  12,   3,   1,  # Lower bounds
   79,  42,  22,   9   # Upper bounds
), nrow = 2, byrow = TRUE)
manual_ranges <- matrix(c(   # Define manual mask numbers
    4,  28,  47,  54,  # Lower bounds
   11,  41,  60,  62   # Upper bounds
), nrow = 2, byrow = TRUE)
model_ranges <- matrix(c(   # Define model mask numbers
    4,  28,  46,  51,  # Lower bounds
    9,  42,  60,  62   # Upper bounds
), nrow = 2, byrow = TRUE)

selected_stations <- c(1, 2, 4, 6)
for (s in 1:4) {
  element_data <- sample_list[[selected_stations[s]]]
  manual_mask <- dilated_manual_masks[[selected_stations[s]]]
  model_mask <- dilated_model_masks[[selected_stations[s]]]
  cell_image <- cell_images[[selected_stations[s]]]
  x <- x_cor[1, s]:x_cor[2, s]
  y <- y_cor[1, s]:y_cor[2, s]
  mask_palette <- distinctColorPalette(10)
  palette <- hcl.colors(100)
  x_len <- length(x)
  y_len <- length(y)
  aspect_ratio <- x_len / y_len
  base_height <- 500
  base_width <- base_height * aspect_ratio * 3
  filename <- paste0("Cell_Element_Maps_", s, ".png")
  png(file.path(base_dir, "Output", filename), width = base_width, height = base_height, res = 150)
  layout(matrix(1:12, ncol = 6))
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
  rasterImage(cell_image, 0, 0, 1, 1)
  mtext("Optical", side = 3, line = -1.5, col = "white", adj = 0.1, cex = 1)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  image(x, y, manual_mask[x, y], zlim = c(manual_ranges[1, s], manual_ranges[2, s]),
        col = mask_palette, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  mtext("Manual", side = 3, line = -1.5, col = "black", adj = 0.1, cex = 1)
  
  image(x, y, model_mask[x, y], zlim = c(model_ranges[1, s], model_ranges[2, s]),
        col = mask_palette, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  mtext("Model", side = 3, line = -1.5, col = "black", adj = 0.1, cex = 1)
  element_indices <- c(1, 5, 2, 6, 3, 7, 4, 8)
  for (i in seq_along(element_indices)) {
    index <- element_indices[i]
    clipped <- pmin(pmax(element_data[x, y, index], limits[1, index]), limits[2, index])
    image(x, y, clipped, zlim = c(limits[1, index], limits[2, index]),
          col = palette, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    mtext(elements[i], side = 3, line = -2, col = "white", adj = 0.1, cex = 1)
  }
  dev.off()
}

#### Calculating Ionomes of Sample Marine Phytoplankton ####

# Choose which segmentation approach to isolate cellular material
mask_lists <- list(Manual = dilated_manual_masks, Model = dilated_model_masks)

# Remove masks that contains P less than the 5th percentile and do not likely contain cellular material 
# and masks that contain Mn|Fe|Zn greater than the mean + 15*sdev as contaminated material
sum_p_all_model <- c()  
tm_values_all <- vector("list", 4) 
for (s in seq_along(sample_list)) {
  element_data <- sample_list[[s]]  
  mask_layer <- mask_lists$Model[[s]]  
  if (is.null(element_data) || is.null(mask_layer)) next  
  sum_p_values <- sapply(1:max(mask_layer, na.rm = TRUE), function(cell_id) {   
    mask_indices <- which(mask_layer == cell_id, arr.ind = TRUE)
    sum(element_data[,,4][mask_indices], na.rm = TRUE)
  })
  sum_p_all_model <- c(sum_p_all_model, sum_p_values)
  for (i in 1:4) {
    tm_values_all[[i]] <- c(tm_values_all[[i]], element_data[,,(i+4)][mask_layer > 0])
  }
}   # Set a minimum P threshold as the 5th percentile of P quotas across all cell masks
p_0.05 <- quantile(sum_p_all_model, 0.05, na.rm = TRUE)
# Set an upper TM threshold to remove cells with particulate TM contamination
tm_mean <- sapply(tm_values_all, function(x) mean(x, na.rm = TRUE))
tm_sd <- sapply(tm_values_all, function(x) sd(x, na.rm = TRUE))
tm_limits_global <- tm_mean + 15 * tm_sd

# Calculate phytoplankton cell ionomes as summing the element concentrations from all pixels 
# as isolated by the segmentation masks normalized to the P content. Cells that contain insufficient
# P or contain trace metal contaminant 'hot-spots' are removed.
processed_cells <- list()  
cell_labels <- list()  
for (segmentation_type in names(mask_lists)) {  
  mask_list <- mask_lists[[segmentation_type]]  
  
  for (s in seq_along(sample_list)) {
    element_data <- sample_list[[s]]  
    mask_layer <- mask_list[[s]]  
    if (is.null(element_data) || is.null(mask_layer)) next  
    cell_n <- max(mask_layer, na.rm = TRUE)  
    station_name <- paste0("S", s)
    valid_data <- list()
    valid_row_names <- list()
    
    for (cell_id in 1:cell_n) {
      mask_indices <- which(mask_layer == cell_id, arr.ind = TRUE)  
      size <- length(mask_indices)
      # Checks cell mask for minimum P quota
      p_cell <- element_data[,,4][mask_indices]
      sum_p <- sum(p_cell, na.rm = TRUE)
      if (sum_p < p_0.05){      
        message(paste(segmentation_type, station_name, "Cell", cell_id, "insufficient P")) 
      next  # Skip this cell
      }
      # Checks cell mask for Trace Metal contamination 
      tm_cell <- sapply(5:8, function(idx) element_data[,,idx][mask_indices], simplify = FALSE)
      tm_cell <- do.call(cbind, tm_cell)
      TM_hotspot <- c()
      if (any(tm_cell[, 1] > tm_limits_global[1])) TM_hotspot <- c(TM_hotspot, "Mn")
      if (any(tm_cell[, 2] > tm_limits_global[2])) TM_hotspot <- c(TM_hotspot, "Fe")
      if (any(tm_cell[, 4] > tm_limits_global[4])) TM_hotspot <- c(TM_hotspot, "Zn")
      if (length(TM_hotspot) > 0) {
        message(paste(segmentation_type, station_name, "Cell", cell_id, 
                      "Contaminated by", paste(TM_hotspot, collapse = ", ")))
      next  # Skip this cell
      }
      sum_elements <- sapply(c(1:3, 5:8), function(idx) sum(element_data[,,idx][mask_indices], na.rm = TRUE))
      cellular_quotas <- sum_elements / sum_p
      valid_data[[length(valid_data) + 1]] <- c(cellular_quotas, size)
      valid_row_names[[length(valid_row_names) + 1]] <- paste0(station_name, "_Cell", cell_id)
    }
    if (length(valid_data) > 0) {
      processed_cells[[segmentation_type]][[s]] <- do.call(rbind, valid_data)
      cell_labels[[segmentation_type]][[s]] <- unlist(valid_row_names)
    }
  }
}
# Process cell ionome data from the different segmentation methods, and format into a structured dataframe, 
# assigning labels using the sampling station and mask number, and exports a CSV.
combined <- list() 
for (segmentation_type in names(processed_cells)) {
  cell_ionomes <- do.call(rbind, processed_cells[[segmentation_type]])
  cell_IDs <- unlist(cell_labels[[segmentation_type]])
  rownames(cell_ionomes) <- cell_IDs
  colnames(cell_ionomes) <- c("Na:P", "Mg:P", "Si:P", "Mn:P", "Fe:P", "Cu:P", "Zn:P", "Mask Size")
  df <- as.data.frame(cell_ionomes)
  df_ionomes <- df[complete.cases(df), ]
  df_ionomes$Station <- sub("_.*", "", rownames(df_ionomes))  
  df_ionomes$Cell <- sub(".*_", "", rownames(df_ionomes))
  df_ionomes$Station <- factor(df_ionomes$Station, levels = unique(df_ionomes$Station))
  levels(df_ionomes$Station) <- paste("S", 1:length(levels(df_ionomes$Station)), sep = "")
  df_ionomes$Segmentation <- segmentation_type
  combined[[segmentation_type]] <- df_ionomes
}

manual_cell_info <- read.csv(file.path(base_dir,"Species and Volume","Manual Cell Data.csv"), header = TRUE)
combined$Manual$Species <- manual_cell_info[, 3]
combined$Manual$Volume <- manual_cell_info[, 4]

model_cell_info <- read.csv(file.path(base_dir,"Species and Volume","Model Cell Data.csv"), header = TRUE)
combined$Model$Species <- model_cell_info[, 3]
combined$Model$Volume <- model_cell_info[, 4]

# Filter out Manually and U-Net Segmented Masks that could not be Manually Identified
combined$Manual <- combined$Manual[combined$Manual$Species != 0, ]
combined$Model <- combined$Model[combined$Model$Species != 0, ]

ionomes_combined <- rbind(combined$Manual, combined$Model)
ionomes_data <- ionomes_combined[, c("Station", "Cell", "Segmentation", "Species", "Volume",
                                     "Na:P", "Mg:P", "Si:P", "Mn:P", "Fe:P", "Cu:P", "Zn:P")]
write.csv(ionomes_data, file.path(file.path(base_dir, "Output","Total Cell Ionomes.csv")), row.names = FALSE)

#### Phytoplankton Ionome Stats ####

# log10 transform the cellular ionomes and make sampling station and segmentation factors
ionomes <- ionomes_combined %>%
  dplyr::select(Station, Cell, `Mask Size`, Segmentation, `Na:P`, `Mg:P`, `Si:P`, `Mn:P`, `Fe:P`, `Cu:P`, `Zn:P`) %>%
  mutate(across(`Na:P`:`Zn:P`, ~ ifelse(. > 0, ., NA))) %>% 
  mutate(across(`Na:P`:`Zn:P`, log10)) %>%
  filter(if_all(`Na:P`:`Zn:P`, is.finite))
ionomes$Station <- as.factor(ionomes$Station)
ionomes$Segmentation <- as.factor(ionomes$Segmentation) 

# Run  two-way ANOVA's for each log10(element:P) and mask size across sampling station and segmentation approach
anova_NaP <- aov(`Na:P` ~ Segmentation * Station, data = ionomes)
anova_MgP <- aov(`Mg:P` ~ Segmentation * Station, data = ionomes)
anova_SiP <- aov(`Si:P` ~ Segmentation * Station, data = ionomes)
anova_MnP <- aov(`Mn:P` ~ Segmentation * Station, data = ionomes)
anova_FeP <- aov(`Fe:P` ~ Segmentation * Station, data = ionomes)
anova_CuP <- aov(`Cu:P` ~ Segmentation * Station, data = ionomes)
anova_ZnP <- aov(`Zn:P` ~ Segmentation * Station, data = ionomes)
anova_masksize <- aov(`Mask Size` ~ Segmentation * Station, data = ionomes)
summary(anova_NaP)
summary(anova_MgP)
summary(anova_SiP)
summary(anova_MnP)
summary(anova_FeP)
summary(anova_CuP)
summary(anova_ZnP)
summary(anova_masksize)

# Load in Cell Species and Volume data to compile with U-net model segmented cell ionomes
# Log10 transform the element:P data and make sampling station and cell species factors
model <- ionomes_combined %>%
  filter(Segmentation == "Model") %>%
  dplyr::select(Station, Cell, Volume, Segmentation, Species, `Na:P`, `Mg:P`, `Si:P`, `Mn:P`, `Fe:P`, `Cu:P`, `Zn:P`) %>%
  mutate(across(Species:`Zn:P`, ~ ifelse(. > 0, ., NA))) %>%
  mutate(across(`Na:P`:`Zn:P`, log10)) %>%
  filter(if_all(Species:`Zn:P`, is.finite))
model$Station <- as.factor(model$Station)  
model$Species <- as.factor(model$Species)

# Run two-way ANOVA's for each log10(element:P) and Cell Volume across Sampling Station and Species
anova_NaP <- aov(`Na:P` ~ Species * Station, data = model)
anova_MgP <- aov(`Mg:P` ~ Species * Station, data = model)
anova_SiP <- aov(`Si:P` ~ Species * Station, data = model)
anova_MnP <- aov(`Mn:P` ~ Species * Station, data = model)
anova_FeP <- aov(`Fe:P` ~ Species * Station, data = model)
anova_CuP <- aov(`Cu:P` ~ Species * Station, data = model)
anova_ZnP <- aov(`Zn:P` ~ Species * Station, data = model)
anova_Vol <- aov(Volume ~ Species * Station, data = model)
summary(anova_NaP)
summary(anova_MgP)
summary(anova_SiP)
summary(anova_MnP)
summary(anova_FeP)
summary(anova_CuP)
summary(anova_ZnP)
summary(anova_Vol)

# Run a one-way ANOVA on the U-net model segmented cell ionomes across sampling station
# and perform a post-HOC Tukey HSD test to identify similar groups of sampling station for each element:P
variables_to_analyze <- c("Na:P", "Mg:P", "Si:P", "Mn:P", "Fe:P", "Cu:P", "Zn:P")
tukey_groups <- perform_anova_and_tukey_all(variables_to_analyze, model)

#### Cell Volume Distributions and Community Composition ####

# Load in Cell Volume and Species for the U-net segmented cell masks
cell_info <- read.csv(file.path(base_dir, "Species and Volume", "Model Cell Data.csv"), header = TRUE)
data <- data.frame(cell_info)
# Assign identified Species labels to the data set key and make species a factor
species_labels <- c(
  "0" = "Unknown",
  "1" = expression(italic("C. atlanticus")),
  "2" = expression(italic("C. brevis")),
  "3" = expression(italic("C. convolutus")),
  "4" = expression(italic("C. dichaeta")),
  "5" = "Pseudonitzschia")
data$Species <- as.factor(data$Species)
data$Station <- as.factor(data$Station)
data$Station_numeric <- as.numeric(as.factor(data$Station))

# Calculate species counts per sampling station, species community percentages, 
# totals cells per station, and orders species into predefined factor levels.
species_counts <- data %>%
  group_by(Station, Species) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(Station, desc(count))  
species_counts <- species_counts %>%
  group_by(Station) %>%
  mutate(percent = count / sum(count),
         total_cells = sum(count)) 
species_counts$Species <- factor(species_counts$Species, 
                                 levels = c("1", "2", "3", "4", "5", "0"))
# Set a common theme to use for Cell Volume box plots and Community Composition bars
common_theme <- theme_minimal() + 
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),  # Add black border
    axis.line = element_line(color = "black"),  # Ensure visible axis lines
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.background = element_rect(fill = "white", color = NA)  # Ensure white background
  )
# Create a stacked bar chart plot to display changes in phytoplankton community across the sampling transect
species_plot <- ggplot(species_counts, aes(x = factor(Station), y = percent, fill = Species)) +
  geom_bar(stat = "identity") +
  labs(x = "Station", y = expression("Community composition (>5µm)"), fill = "Species") +
  scale_fill_brewer(palette = "Set1", labels = species_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  geom_text(aes(label = paste("n=", total_cells), y = 1), vjust = -0.5, size = 4, color = "black") +
  theme_minimal() + common_theme

# Fit a linear model to the log10 transformed cell volume data aross the sampling transect
lm_fit <- lm(log10(Volume) ~ Station_numeric, data = subset(data, Species != 0))
m <- coef(lm_fit)[2]  # Slope
b <- coef(lm_fit)[1]  # Intercept
equation_text <- paste0("y = ", round(10^b), "x^", format(round(m, 3), nsmall = 3))
r_squared_text <- paste0("R² = ", round(summary(lm_fit)$r.squared, 3))
summary(lm_fit)
vol_plot <- ggplot(subset(data, Species != 0), aes(x = Station, y = log10(Volume))) + 
  geom_boxplot(outlier.shape = NA, fill = "#EEE9E9") +
  geom_point(position = position_jitter(width = 0.1), color = "black", alpha = 0.5) + 
  geom_smooth(aes(x = Station_numeric, y = log10(Volume)), method = "lm", color = "red", 
              linewidth = 1, linetype = "dashed", se = FALSE) +  
  labs(y = expression('Cell volume (µm'^'3'~')'), x = "Station") +
  scale_y_continuous(breaks = seq(2, 5, by = 1), 
                     labels = c(expression('10'^'2'), expression('10'^'3'), 
                                expression('10'^'4'), expression('10'^'5'))) +
  coord_cartesian(ylim = c(2, 5)) +  
  scale_x_discrete(breaks = levels(data$Station)) +
  theme_bw() + 
  common_theme +  
  geom_text(aes(x = 1, y = 4.8, label = equation_text), color = "red", size = 5, hjust = 0) +  
  geom_text(aes(x = 1, y = 4.6, label = r_squared_text), color = "red", size = 5, hjust = 0)
final_plot <- vol_plot + species_plot + plot_layout(ncol = 2)
ggsave(file.path(base_dir,"Output","Cell Volumes and Community Composition.png"), plot = final_plot, width = 10, height = 6, dpi = 300)

# Create a distance matrix, performs hierarchical clustering, beta dispersion analysis, heatmap, and PERMANOVA for community comparison across the trasnect
species_counts_unique <- species_counts %>%
  group_by(Station, Species) %>%
  summarise(count = sum(count), .groups = "drop")
species_matrix <- species_counts_unique %>%
  pivot_wider(names_from = Species, values_from = count, values_fill = 0) %>%
  column_to_rownames("Station")
colnames(species_matrix) <- c("C.atlanticus","C.brevis","C.convolutus","C.dichaeta","Pseudonitzschia","Unknown")
distance_matrix <- vegdist(species_matrix, method = "bray")
stations <- factor(c("S1", "S2", "S3", "S4", "S5", "S6"))  # Adjust if more stations
hclust_result <- hclust(distance_matrix, method = "average")
plot(hclust_result)
betadisper_result <- betadisper(distance_matrix, stations)
anova(betadisper_result)  # Perform permutation test on dispersion
heatmap(as.matrix(distance_matrix), main = "Bray-Curtis Distance Matrix", scale = "none")
permanova_result <- adonis2(distance_matrix ~ stations, permutations = 999)
print(permanova_result)


#### Cell Ionome Boxplots ####

# Create box plots to show log-10 transformed element:P distributions across sampling station 
# and segmentation approach, labelled with similarity group lettering from a post-HOC Tukey-HSD test
y_limits <- list("Na:P" = c(-2.0, 2.2), "Mg:P" = c(-2.0, 2.2), "Si:P" = c( 1.0, 4.0),
                 "Mn:P" = c(-5.0,-1.3), "Fe:P" = c(-4.0, 0.5), "Cu:P" = c(-5,-1.7), "Zn:P" = c(-4.0, 0.5))
y_breaks <- list("Na:P" = seq(-2, 2, 1), "Mg:P" = seq(-2, 2, 1), "Si:P" = seq( 1, 4, 1),
                 "Mn:P" = seq(-5,-2, 1), "Fe:P" = seq(-4, 0, 1), "Cu:P" = seq(-5,-2, 1), "Zn:P" = seq(-4, 0, 1))
plot_list <- lapply(variables_to_analyze, function(element) {
  plot_boxplot(element, ionomes, tukey_groups, y_limits, y_breaks)
})
final_plot <- wrap_plots(plot_list, ncol = 4) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(base_dir, "Output","Cell Ionome Boxplots.png"), plot = final_plot, width = 10, height = 6, dpi = 300)


#### Station, Species, and Size PCA plots ####

# Perform Principle Component Analysis (PCA) on cellular ionome ratios and visualized 
# against the first two principal components, labelled by station, species, and cell size.

ionomes <- model %>% select(`Na:P`:`Zn:P`)
pca_ionomes <- prcomp(ionomes, scale. = TRUE)
ionomes_pca_df <- data.frame(pca_ionomes$x, Station = model$Station, Species = model$Species, Volume = model$Volume)

# PCA Plot for Station
ionomes_pca_df$Station <- as.factor(ionomes_pca_df$Station)
pca_station_plot <- autoplot(pca_ionomes, data = ionomes_pca_df, colour = 'Station',shape = 'Station', 
                             size = 5,loadings = TRUE, loadings.colour = 'black',loadings.label = TRUE,loadings.label.size = 5,loadings.label.colour = 'black',loadings.lwd = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top") +
  scale_shape_manual(values = c(15, 16, 17, 18, 19, 20)) +
  labs(title = "(a) Station")
p_station <- ggplot_build(pca_station_plot)$plot
p_station$layers[[2]]$aes_params$size <- 1.2

# PCA Plot for Species
ionomes_pca_df$Species <- factor(ionomes_pca_df$Species,levels = c(1, 2, 3, 4, 5), 
                                 labels = c("C. atlanticus","C. brevis","C. convolutus","C. dichaeta","Pseudonitzschia"))
pca_species_plot <- autoplot(pca_ionomes, data = ionomes_pca_df,colour = 'Species',shape = 'Species',   
                             size = 5,loadings = TRUE,loadings.colour = 'black',loadings.label = TRUE,loadings.label.size = 5,loadings.label.colour = 'black',loadings.lwd = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "top",legend.text = element_text(face = "italic"), legend.direction = "vertical") +
  scale_shape_manual(values = c(15, 16, 17, 18, 19)) +
  labs(title = "(b) Species")
p_species <- ggplot_build(pca_species_plot)$plot
p_species$layers[[2]]$aes_params$size <- 1.2  # Increase vector thickness

# PCA Plot for Volume
ionomes_pca_df <- ionomes_pca_df %>%
  mutate(Volume = cut(Volume, 
                      breaks = c(0, 500, 1000, 2000, 3000, 5000, Inf),
                      labels = c("<500 µm³","500-1000 µm³","1000-2000 µm³","2000-3000 µm³","3000-5000 µm³",">5000 µm³"),include.lowest = TRUE))
pca_volume_plot <- autoplot(pca_ionomes, data = ionomes_pca_df,colour = 'Volume',shape = 'Volume', 
                            size = 5,loadings = TRUE,loadings.colour = 'black',loadings.label = TRUE,loadings.label.size = 5,loadings.label.colour = 'black',loadings.lwd = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "top") +
  scale_shape_manual(values = c(15, 16, 17, 18, 19, 20)) +
  labs(title = "(c) Volume")
p_volume <- ggplot_build(pca_volume_plot)$plot
p_volume$layers[[2]]$aes_params$size <- 1.2

# Generate a PCA panel plot: a = Station, b = Species, c = Volume
final_plot <- (p_station | p_species | p_volume) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.direction = "vertical")
ggsave(file.path(base_dir,"Output","PCA Ionome plots Station Species Volume.png"), plot = final_plot, width = 15, height = 6, dpi = 300)
