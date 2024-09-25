
run_program <- function(pdb_id,entropy_file_path, energy_file_path){
  
  
  
  # Load required libraries
  library(bio3d)
  library(cluster)
  library(ggplot2)
  library(umap)
  library(readxl)
  library(dplyr)
  library(geometry)
  library(MASS)  # For Mahalanobis distance calculation
  library(FactoMineR)  # For PCA
  library(dbscan)  # For DBSCAN clustering
  library(sf)  # For minimum bounding rectangle (MBR)
  library(ggrepel)  # For better label placement
  library(neuralnet)
  library(gridExtra)
  library(reshape2)
  set.seed(123)  # Set a random seed
 
  # Unified function to process a PDB ID
  process_pdb <- function(pdb_id, chain_id = NULL, num_clusters = 5, output_dir = getwd() ) {
    
    # Fetch the PDB file from the internet
    pdb <- read.pdb(pdb_id)
    
    # Check if PDB is read correctly
    if (is.null(pdb)) {
      stop(paste("Failed to read PDB ID:", pdb_id))
    }
    
    # Automatically select the chain(s) in the PDB file
    chain_ids <- unique(pdb$atom$chain)
    if (length(chain_ids) == 0) {
      stop(paste("No chains found in PDB ID:", pdb_id))
    }
    
    # Prompt user to select a chain if more than one is present
    if (length(chain_ids) > 1) {
      cat("Multiple chains found in PDB ID:", pdb_id, "\n")
      cat("Available chains:", paste(chain_ids, collapse = ", "), "\n")
      chain_id <- readline(prompt = "Please enter the chain ID you wish to process: ")
      
      # Check if the input chain ID is valid
      if (!(chain_id %in% chain_ids)) {
        stop(paste("Invalid chain ID selected:", chain_id))
      }
    } else {
      # If only one chain, select it automatically
      chain_id <- chain_ids[1]
      cat("Processing the only available chain:", chain_id, "for PDB ID:", pdb_id, "\n")
    }
    
    # Extract the specified chain data
    chain_data <- trim.pdb(pdb, chain=chain_id)
    
    # Check if chain data is extracted correctly
    if (is.null(chain_data)) {
      stop(paste("Failed to extract chain data for chain:", chain_id, "in PDB ID:", pdb_id))
    }
    
    
    # Calculate torsion angles
    torsion_angles <- torsion.pdb(chain_data)
    
    # Check if torsion angles are calculated correctly
    if (is.null(torsion_angles)) {
      stop(paste("Failed to calculate torsion angles for chain:", chain_id, "in PDB ID:", pdb_id))
    }
    
    # Create a data frame with torsion angles
    torsion_df <- data.frame(
      alpha = torsion_angles$alpha,
      omega = torsion_angles$omega
    )
    torsion_df <- cbind(torsion_angles$tbl, torsion_df)
    
    # Handle missing values by converting NAs to 0
    torsion_df[is.na(torsion_df)] <- 0
    
    # Remove rows where the sum is zero
    torsion_df <- torsion_df[rowSums(torsion_df) != 0, ]
    
    # Remove columns where the sum is zero
    torsion_df <- torsion_df[, colSums(torsion_df) != 0]
    
    # Check if the torsion_df is empty after removing rows/columns with sum zero
    if (nrow(torsion_df) == 0 || ncol(torsion_df) == 0) {
      stop("No valid torsion angles available after removing rows and columns with sum zero.")
    }
    
    # Write the torsion angles to a CSV file
    torsion_file <- file.path(output_dir, paste0("torsion_", pdb_id, ".csv"))
    write.csv(torsion_df, torsion_file, row.names = FALSE)
    message(paste("Torsion angles saved to:", torsion_file))
    
    # Transpose the data for correlation
    torsion_data_t <- t(torsion_df)
    
    # Calculate Pearson correlation matrix
    pearson_matrix <- cor(torsion_data_t, method = "pearson", use = "complete.obs")
    
    # Perform UMAP embedding
    umap_result <- tryCatch({
      umap(pearson_matrix)
    }, error = function(e) {
      stop("UMAP embedding failed. This could be due to insufficient variation in the data or invalid correlation matrix.")
    })
    
    umap_coordinates <- umap_result$layout
    
    # Create a data frame with UMAP coordinates
    umap_df <- data.frame(x = umap_coordinates[, 1], y = umap_coordinates[, 2])
    
    
    kmeans_result <- kmeans(umap_coordinates, centers = num_clusters)
    
    # Add cluster assignments to the UMAP coordinates
    umap_df$cluster <- as.factor(kmeans_result$cluster)
    
    # Calculate silhouette scores
    silhouette_scores <- silhouette(kmeans_result$cluster, dist(umap_coordinates))
    
    # Calculate average silhouette score
    avg_silhouette_score <- mean(silhouette_scores[, 'sil_width'])
    
    # Merge UMAP and k-means plot into one and add silhouette score to the title
    umap_kmeans_plot <- ggplot(umap_df, aes(x = x, y = y, color = cluster)) +
      geom_point() +
      ggtitle(paste("UMAP Clustering with k-means for", pdb_id, "- Avg Silhouette:", round(avg_silhouette_score, 2))) +
      xlab("UMAP 1") +
      ylab("UMAP 2") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_point(aes(x = 0, y = 0), color = "blue", size = 3) +
      theme_minimal() +
      coord_fixed(ratio = 1) +
      theme(axis.title = element_text(face = "bold", size = 12))  # Make axis labels bold
    
    umap_kmeans_plot_file <- file.path(output_dir, paste0(pdb_id, "_umap_kmeans_plot.png"))
    ggsave(umap_kmeans_plot_file, umap_kmeans_plot, width = 8, height = 6)
    message(paste("UMAP k-means plot saved to:", umap_kmeans_plot_file))
    silhouette_df <- data.frame("PDB" = pdb_id, "avg_silhouette_score" = avg_silhouette_score)
    # Save results in variables dynamically
    assign(paste0("tor_",pdb_id), torsion_df, envir = .GlobalEnv)
    assign(paste0("pearson_matrix_", pdb_id), pearson_matrix, envir = .GlobalEnv)
    assign(paste0("umap_result_", pdb_id), umap_result, envir = .GlobalEnv)
    assign(paste0("kmeans_result_", pdb_id), kmeans_result, envir = .GlobalEnv)
    assign(paste0("clustered_data_", pdb_id, "_umap"), umap_df, envir = .GlobalEnv)
    assign("silhouette_df", silhouette_df, envir = .GlobalEnv)
  }
  
  
  
  
  load_and_merge_data <- function(pdb_id, entropy_file_path, energy_file_path) {
    cat("\nProcessing PDB ID:", pdb_id, "\n")
    
    # Load data from the provided Excel file paths
    entropy_data <- tryCatch({
      read_excel(entropy_file_path, sheet = "Sheet1")
    }, error = function(e) {
      warning(paste("Error loading entropy data for PDB ID:", pdb_id))
      return(NULL)
    })
    
    energy_data <- tryCatch({
      read_excel(energy_file_path, sheet = "Sheet1")
    }, error = function(e) {
      warning(paste("Error loading energy data for PDB ID:", pdb_id))
      return(NULL)
    })
    
    if (is.null(entropy_data) || is.null(energy_data)) {
      warning(paste("Missing entropy or energy data for PDB ID:", pdb_id))
      return(NULL)
    }
    
    # Check if tor_ and clustered_data_ objects exist
    tor_name <- paste0("tor_", pdb_id)
    clustered_data_name <- paste0("clustered_data_", pdb_id, "_umap")
    
    if (exists(tor_name, envir = .GlobalEnv) && exists(clustered_data_name, envir = .GlobalEnv)) {
      # Load tor and clustered data
      tor_data <- data.frame(row_names = rownames(get(tor_name)), get(tor_name), stringsAsFactors = FALSE)
      clustered_data <- get(clustered_data_name, envir = .GlobalEnv)
      
      if (nrow(tor_data) == 0) {
        warning(paste("tor_data is empty for PDB ID:", pdb_id))
        return(NULL)
      }
      if (nrow(clustered_data) == 0) {
        warning(paste("clustered_data is empty for PDB ID:", pdb_id))
        return(NULL)
      }
      
      # Trim whitespace from relevant columns
      tor_data$row_names <- trimws(tor_data$row_names)
      rownames(clustered_data) <- trimws(rownames(clustered_data))
      entropy_data$residue <- trimws(entropy_data$residue)
      energy_data$Variable <- trimws(energy_data$Variable)
      
      # Rename columns to ensure consistency
      colnames(entropy_data)[colnames(entropy_data) == "residue"] <- "Variable"
      
      # Print column names for diagnostic purposes
      cat("Column names for", tor_name, ":\n", colnames(tor_data), "\n")
      cat("Column names for", clustered_data_name, ":\n", colnames(clustered_data), "\n")
      cat("Column names for entropy_data:\n", colnames(entropy_data), "\n")
      cat("Column names for energy_data:\n", colnames(energy_data), "\n")
      
      # Merge tor data with clustered data based on 'Variable' column
      merged_data <- merge(tor_data, clustered_data, by.x = "row_names", by.y = "row.names")
      
      if (nrow(merged_data) == 0) {
        warning(paste("Merged data is empty after merging tor_data and clustered_data for PDB ID:", pdb_id))
        return(NULL)
      }
      
      # Rename columns in merged_data if needed
      colnames(merged_data)[colnames(merged_data) == "row_names"] <- "Variable"
      
      # Merge with entropy_data
      merged_data <- merge(merged_data, entropy_data, by = "Variable")
      
      if (nrow(merged_data) == 0) {
        warning(paste("Merged data is empty after merging with entropy_data for PDB ID:", pdb_id))
        return(NULL)
      }
      
      # Merge with energy_data
      merged_data <- merge(merged_data, energy_data, by = "Variable")
      
      if (nrow(merged_data) == 0) {
        warning(paste("Merged data is empty after merging with energy_data for PDB ID:", pdb_id))
        return(NULL)
      }
      
      # Store the merged data frame in a variable within the global environment
      assign("merged_df", merged_data, envir = .GlobalEnv)
      
      # Return the merged data frame (optional, if you want to use it immediately)
      return(merged_data)
    } else {
      warning(paste("Missing tor or clustered data for PDB ID:", pdb_id))
      return(NULL)  # Return NULL if required data is missing
    }
  }
  
  
  
  
  extract_cluster_sums_and_averages <- function(pdb_data) {
    # Group by cluster and summarize entropy and energy
    df_sums <- pdb_data %>%
      group_by(cluster) %>%
      summarize(sum_entropy = sum(PackingEntropy, na.rm = TRUE),
                sum_energy = sum(Energy, na.rm = TRUE),
                .groups = 'drop')  # Avoid warning about grouped data
    
    # Group by cluster and calculate average entropy and energy
    df_averages <- pdb_data %>%
      group_by(cluster) %>%
      summarize(avg_entropy = mean(PackingEntropy, na.rm = TRUE),
                avg_energy = mean(Energy, na.rm = TRUE),
                .groups = 'drop')  # Avoid warning about grouped data
    
    # Create a combined data frame for sums and averages
    cluster_sum_average <- merge(df_sums, df_averages, by = "cluster")
    
    # Store the result in the global environment
    assign("cluster_sum_average", cluster_sum_average, envir = .GlobalEnv)
    
    # Optionally return the result (if you want to use it immediately)
    return(cluster_sum_average)
  }
  
  transform_to_row <- function(pdb_id, data) {
    
    # Create a named vector with the required format
    result_1 <- c(
      pdb_id = pdb_id,
      entropy_sum_1 = data[1, "sum_entropy"],
      entropy_sum_2 = data[2, "sum_entropy"],
      entropy_sum_3 = data[3, "sum_entropy"],
      entropy_sum_4 = data[4, "sum_entropy"],
      entropy_sum_5 = data[5, "sum_entropy"],
      energy_sum_1 = data[1, "sum_energy"],
      energy_sum_2 = data[2, "sum_energy"],
      energy_sum_3 = data[3, "sum_energy"],
      energy_sum_4 = data[4, "sum_energy"],
      energy_sum_5 = data[5, "sum_energy"],
      entropy_avg_1 = data[1, "avg_entropy"],
      entropy_avg_2 = data[2, "avg_entropy"],
      entropy_avg_3 = data[3, "avg_entropy"],
      entropy_avg_4 = data[4, "avg_entropy"],
      entropy_avg_5 = data[5, "avg_entropy"],
      energy_avg_1 = data[1, "avg_energy"],
      energy_avg_2 = data[2, "avg_energy"],
      energy_avg_3 = data[3, "avg_energy"],
      energy_avg_4 = data[4, "avg_energy"],
      energy_avg_5 = data[5, "avg_energy"]
    )
    
    # Convert the result to a data frame (single row)
    result_row <- as.data.frame(t(result_1), stringsAsFactors = FALSE)
    
    # Assign to the global environment
    assign("result_row", result_row, envir = .GlobalEnv)
    
    return(result_row)  # Optional: return the result as well
  }
  
  
  Convex_hull <- function(entropy_energy, pdb_id) {
    library(ggplot2)
    library(ggrepel)
    library(pracma)
    
    # Extract the corresponding dataframe
    df <- entropy_energy
    
    # Check if the dataframe is not empty
    if (is.null(df) || nrow(df) == 0) {
      stop(paste("No data found for PDB ID:", pdb_id))
    }
    
    # Calculate the convex hull
    hull_indices <- chull(df$sum_entropy, df$sum_energy)
    hull_points <- df[hull_indices, c("sum_entropy", "sum_energy")]
    
    # Calculate the area of the convex hull
    hull_area <- abs(pracma::polyarea(hull_points$sum_entropy, hull_points$sum_energy))
    
    # Calculate the centroid
    centroid_entropy <- mean(df$sum_entropy)
    centroid_energy <- mean(df$sum_energy)
    
    # Calculate the MBR area
    mbr_width <- max(df$sum_entropy) - min(df$sum_entropy)
    mbr_height <- max(df$sum_energy) - min(df$sum_energy)
    mbr_area <- mbr_width * mbr_height
    
    # Calculate the ratio of MBR area to convex hull area
    mbr_to_hull_ratio <- mbr_area / hull_area
    
    # Create the plot
    p <- ggplot(df, aes(x = sum_entropy, y = sum_energy)) +
      geom_point(color = "black", alpha = 0.6) +
      geom_polygon(data = hull_points, aes(x = sum_entropy, y = sum_energy), fill = "red", alpha = 0.3, color = "darkred") +
      geom_rect(aes(xmin = min(df$sum_entropy), xmax = max(df$sum_entropy), 
                    ymin = min(df$sum_energy), ymax = max(df$sum_energy)), 
                fill = NA, color = "blue", linetype = "dashed") +
      geom_point(aes(x = centroid_entropy, y = centroid_energy), color = "purple", size = 4, shape = 4) +
      geom_text_repel(aes(x = centroid_entropy, y = centroid_energy, label = "Centroid"), 
                      color = "purple", fontface = "bold") +
      labs(title = paste("Geometric Measures for", pdb_id),
           x = "Sum of Entropy", y = "Sum of Energy") +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 14, color = "black"),  # X-axis label bold and black
        axis.title.y = element_text(face = "bold", size = 14, color = "black"),  # Y-axis label bold and black
        axis.text = element_text(size = 12, face = "bold", color = "black"),     # Axis tick labels bold and black
        axis.line = element_line(size = 1.2, color = "black"),                   # Make axis lines bold and black
        legend.position = "none",
        axis.title = element_blank()                                             # Remove axis titles for paper
      ) +
      annotate("text", x = min(df$sum_entropy), y = max(df$sum_energy), 
               label = paste("Convex Hull Area:", round(hull_area, 2)),
               hjust = 0, vjust = 1, color = "darkred", fontface = "bold") +
      annotate("text", x = min(df$sum_entropy), y = max(df$sum_energy) - (max(df$sum_energy) - min(df$sum_energy))*0.05, 
               label = paste("MBR Area:", round(mbr_area, 2)),
               hjust = 0, vjust = 1, color = "blue", fontface = "bold")
    
    # Save the plot
    ggsave(paste0(pdb_id, "_enhanced_geometric_measures_plot.png"), plot = p, width = 10, height = 8, dpi = 300)
    
    # Create and return a summary of calculations
    result <- data.frame(
      pdb_id = pdb_id,
      hull_area = hull_area,
      mbr_area = mbr_area,
      mbr_to_hull_ratio = mbr_to_hull_ratio,
      centroid_entropy = centroid_entropy,
      centroid_energy = centroid_energy
    )
    
    # Print the summary results
    print(result)
    
    # Return the result
    return(result)
  }
  
  
  
  
  
  compute_eigenvalues_from_merged_data <- function(merged_df) {
    
    # Step 1: Split dataframe by 'cluster' column
    sub_dataframe_list <- split(merged_df, merged_df$cluster)
    
    # Step 2: Define columns to keep
    columns_to_keep <- c("phi", "psi", "chi1", "chi2", "chi3", "chi4", "alpha", "omega")
    
    # Step 3: Process each sub-dataframe
    process_dataframe <- function(df) {
      if (!is.null(rownames(df))) {
        df$row_names <- rownames(df)
      } else {
        df$row_names <- seq_len(nrow(df))
      }
      
      df_reduced <- df[, c("row_names", intersect(columns_to_keep, colnames(df)))]
      rownames(df_reduced) <- df_reduced$row_names
      df_reduced$row_names <- NULL
      
      return(df_reduced)
    }
    
    sub_dataframe_list_reduced <- lapply(sub_dataframe_list, process_dataframe)
    
    # Step 4: Compute correlation matrices
    compute_cor_matrix_abs_no_na <- function(df) {
      df_numeric <- df[sapply(df, is.numeric)]
      cor_matrix <- cor(df_numeric, use = "pairwise.complete.obs")
      cor_matrix_abs <- abs(cor_matrix)
      cor_matrix_abs[is.na(cor_matrix_abs)] <- 0
      cor_df <- as.data.frame(cor_matrix_abs)
      cor_df$variable <- rownames(cor_df)
      
      return(cor_df)
    }
    
    cor_matrix_list <- lapply(sub_dataframe_list_reduced, compute_cor_matrix_abs_no_na)
    
    # Step 5: Compute eigenvalues
    compute_eigenvalues <- function(cor_matrix, matrix_name) {
      cor_matrix_numeric <- cor_matrix[, !colnames(cor_matrix) %in% "variable", drop = FALSE]
      eigen_values <- eigen(as.matrix(cor_matrix_numeric))$values
      
      result <- data.frame(
        matrix_name = rep(matrix_name, length(eigen_values)),
        eigenvalue = eigen_values,
        rank = seq_along(eigen_values)
      )
      
      return(result)
    }
    
    eigen_list <- lapply(names(cor_matrix_list), function(name) {
      compute_eigenvalues(cor_matrix_list[[name]], name)
    })
    
    # Step 6: Combine into a single dataframe
    eigen_df <- do.call(rbind, eigen_list)
    
    # Step 7: Reshape to wide format
    eigen_df_wide <- reshape(eigen_df, 
                             idvar = "matrix_name",
                             timevar = "rank",
                             direction = "wide")
    
    colnames(eigen_df_wide) <- c("matrix_name", paste0("eigenvalue_", seq_len(ncol(eigen_df_wide) - 1)))
    
    # Step 8: Sort columns
    eigen_df_wide <- eigen_df_wide[, c("matrix_name", sort(colnames(eigen_df_wide)[-1]))]
    assign("eigen_df", eigen_df_wide, envir = .GlobalEnv)
    return(eigen_df_wide)
  }
  
  
  
  merge_dataframes <- function(cluster_sum_average, eigen_df) {
    # Merge the dataframes based on the 'cluster' column in cluster_sum_average and 'matrix_name' column in eigendf
    merged_df <- merge(cluster_sum_average, eigen_df, by.x = "cluster", by.y = "matrix_name", all = FALSE)
    assign("merged_df", merged_df, envir = .GlobalEnv)
    return(merged_df)
  }
  
  
  # Function to standardize the data
  standardize_data <- function(df, dependent_var, independent_vars) {
    df %>%
      mutate(across(all_of(c(dependent_var, independent_vars)), scale))
  }
  
  create_ann_model <- function(df, dependent_var, independent_vars) {
    # Standardize the data
    df <- standardize_data(df, dependent_var, independent_vars)
    
    # Create the formula for the ANN model
    formula <- as.formula(paste(dependent_var, "~", paste(independent_vars, collapse = " + ")))
    
    # Split the data into training and testing sets
    set.seed(123)  # for reproducibility
    train_indices <- sample(1:nrow(df), 0.7 * nrow(df))
    train_data <- df[train_indices, ]
    test_data <- df[-train_indices, ]
    
    # Fit the ANN model
    ann_model <- neuralnet(formula, data = train_data, hidden = c(5, 3), linear.output = TRUE)
    
    # Predict on test data
    predictions <- compute(ann_model, test_data[, independent_vars])$net.result
    
    # Calculate metrics
    residuals <- test_data[[dependent_var]] - predictions
    rss <- sum(residuals^2)
    rmse <- sqrt(mean(residuals^2))
    n <- nrow(test_data)
    k <- length(independent_vars)
    ser <- sqrt(rss / abs(n - k - 1))
    
    # Scale RSS and SER as in the original code
    scaled_rss <- rss * 10^5
    scaled_ser <- ser * 100
    
    # Create a data frame for the metric values
    metric_values <- data.frame(
      rmse = rmse,
      rss = rss,
      scaled_rss = scaled_rss,
      ser = ser,
      scaled_ser = scaled_ser
    )
    
    # Print the metric values
    print(metric_values)
    
    # Assign a name to store the results in the global environment (optional)
    assign("ann_metrics", metric_values, envir = .GlobalEnv)
    
    return(metric_values)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  process_pdb(pdb_id)
  load_and_merge_data(pdb_id,entropy_file_path, energy_file_path)
  extract_cluster_sums_and_averages(merged_df)
  Convex_hull(cluster_sum_average, pdb_id)
  transform_to_row(pdb_id, cluster_sum_average)
  compute_eigenvalues_from_merged_data(merged_df)
  merge_dataframes(cluster_sum_average, eigen_df)
  # Example usage
  dependent_var <- "sum_entropy"
  independent_vars <- c("eigenvalue_1", "eigenvalue_2", "eigenvalue_3", "eigenvalue_4",
                        "eigenvalue_5", "eigenvalue_6", "eigenvalue_7", "eigenvalue_8")
  # Function to standardize the data
  standardize_data <- function(df, dependent_var, independent_vars) {
    df %>%
      mutate(across(all_of(c(dependent_var, independent_vars)), scale))
  }
  # Assuming you have a dataframe named merged_df
  create_ann_model(merged_df, dependent_var, independent_vars)
  final_result <- cbind(silhouette_df ,result_row, convex_hull_result ,ann_metrics)
  assign("Final_result", final_result, envir = .GlobalEnv)
  analysis_result <- final_result[,c("PDB", "avg_silhouette_score","hull_area","pca_variance", "mbr_area", "rmse", "ser", "scaled_ser" )]
  assign("analysis_result", analysis_result, envir = .GlobalEnv)
  print(analysis_result)
}
  
run_program("5bon", "D:\\2024 account\\all entropy data\\entropy_5bon.xlsx", "D:\\2024 account\\all pdb files\\energy file for all pdb\\energy_5bon.xlsx")
run_program("3s94", "D:\\2024 account\\all entropy data\\entropy_3s94.xlsx", "D:\\2024 account\\all pdb files\\energy file for all pdb\\energy_3s94.xlsx")

run_program("2ZZP", "D:\\2024 account\\new validation data\\all entropy\\entropy_2zzp.xlsx", "D:\\2024 account\\new validation data\\all energy\\energy_2zzp.xlsx")




entropy_file_path <- paste0("D:\\2024 account\\new validation data\\all entropy\\entropy_", pdb_id, ".xlsx")
energy_file_path <- paste0("D:\\2024 account\\new validation data\\all energy\\energy_", pdb_id, ".xlsx")


