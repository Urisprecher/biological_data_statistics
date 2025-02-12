###Statistical Analysis Pipeline
## takes as input a folder with a csv file that must contain at least one descriptive column with 2 or more groups to compare and other numeric columns. 
## all output files will be saved in the main folder under a results folder of your choice. 


##run these lines  
perform_statistical_analysis <- function() {
  #directory to default
  setwd("~")
  cat("\033[1;32m==== Starting statistical analysis!====\033[0m\n")
  message("Please ensure your csv file is correctly formatted and located in a specified folder:)")
  message("Loading required libraries, please wait...")
  #required packages
  required_packages <- c("ggplot2", "dplyr", "pwr", "tidyr", "tibble", "broom",
                         "purrr", "MASS", "caret", "Metrics", "combinat",
                         "devtools", "FSA", "dabestr", "AICcmodavg", "cli",
                         "BayesFactor", "car", "nnet", "pROC", "parallel", "future", "furrr", "corrplot", "mixOmics")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      if (!require(pkg, character.only = TRUE)) {
        stop(paste("Package", pkg, "could not be installed. Check library paths and permissions."))
      }
    }
  }
  #main directory choice and change
  main_dir <- readline("Enter the main folder path where your csv file is located (or press Enter to use the current directory): ")
  if (main_dir != "") {
    if (dir.exists(main_dir)) {
      setwd(main_dir)
      message("Directory set to: ", main_dir)
    } else {
      stop("Error: The specified main directory does not exist.")
    }
  } else {
    message("Proceeding with the current directory: ", getwd())
  }
  #output directory name
  output_dir <- readline(prompt = "Enter name of results folder of your choice: ")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
    message(paste("Directory", output_dir, "created."))
  } else {
    message(paste("Directory", output_dir, "already exists."))
  }
  #subdirectory for multi group statistics
  all_stat_dir <- file.path(output_dir, "all_stats")
  if (!dir.exists(all_stat_dir)) {
    dir.create(all_stat_dir)
  }
  #Load csv file
  load_data <- function(main_dir) {
    #all csv files in the main directory
    csv_files <- list.files(path = main_dir, pattern = "\\.csv$", full.names = TRUE)
    #are there are any csv files in the folder?
    if (length(csv_files) == 0) {
      stop("No CSV files found in the specified directory. Please check the folder and try again.")
    }
    #loading the first csv file found in the directory
    tryCatch({
      data <- read.csv(csv_files[1])
      message("Data loaded successfully from: ", csv_files[1])
      return(data)
    }, error = function(e) {
      stop("Failed to load data. Ensure the file is properly formatted.")
    })
  }
  #Load data
  data <- load_data(getwd())
  cli_alert_success("main folder, outpur folder and csv files were handled successfully!")
  cat("\033[1;31m==== Step 1 - processing data... ====\033[0m\n")
  #optionally remove columns
  cat("\033[1;32m==== do you want to remove any columns from the data ? ====\033[0m\n")
  print("Current columns in the data frame:")
  print(names(data))
  columns_to_remove_pre <- readline(prompt = "Enter names of columns to remove (comma-separated) before statistical analysis, or press enter to skip: ")
  if (columns_to_remove_pre != "") {
    columns_to_remove_pre <- strsplit(columns_to_remove_pre, ",")[[1]]
    if (any(!columns_to_remove_pre %in% names(data))) {
      stop("Some specified columns do not exist in the data.verify the column names.")
    }
    data <- data[, !(names(data) %in% columns_to_remove_pre)]
    write.csv(data, file.path(output_dir, "data_columns_removal.csv"), row.names = FALSE)
  }

  perform_correlation_analysis <- function(data, output_folder) {
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    numeric_data <- data[, sapply(data, is.numeric), drop = TRUE]
    numeric_data <- scale(numeric_data)
    #handle NAs
    correlation_matrix <- cor(numeric_data, use = "pairwise.complete.obs", method = "pearson")
    #save correlation matrix 
    write.csv(correlation_matrix, file.path(output_folder, "correlation_matrix.csv"), row.names = TRUE)
    #plot heatmap using corrplot
    png(file.path(output_folder, "correlation_heatmap.png"), width = 800, height = 800)
    corrplot(correlation_matrix, method = "pie", type = "upper", order="hclust", addCoef.col = "black",
             col = colorRampPalette(c("blue", "white", "red"))(200),
             tl.col = "black", tl.srt = 45)
    dev.off()
    list(correlation_matrix = correlation_matrix)
  }
  plot_folder_cor <- file.path(dirname(all_stat_dir), "correlation_data")
  if (!dir.exists(plot_folder_cor)) {
    dir.create(plot_folder_cor)
  }
  message("generating correlation data...")
  cor_frame <- as.data.frame(data)
  perform_correlation_analysis(data = cor_frame, output_folder = plot_folder_cor)
  
  ###linear regression analysis
  linear_reg_res_dir <- file.path(all_stat_dir, "stat_linear_regression")
  dir.create(linear_reg_res_dir, showWarnings = FALSE)
  
  run_simple_regression <- as.logical(readline(prompt = "Do you want to run per feature linear regression analysis? (TRUE/FALSE): "))
  
  if (run_simple_regression) {
    #index column with validation for linear reg 
    choose_lin_index_column <- function(data) {
      cat("Available columns:\n")
      print(names(data))
      lin_index_col <- readline(prompt = "Enter the name of the index column (should be continous): ")
      if (!(lin_index_col %in% names(data))) {
        stop("The specified index column does not exist.check and try again.")
      }
      return(lin_index_col)
    }
    lin_index_col <- choose_lin_index_column(data)
    
    # linear regression analysis per feature 
    perform_simple_regression_analysis <- function(data, lin_index_col) {
      tryCatch({
        sub_data_filtered <- data
        #dependent variable numeric validation
        if (!is.numeric(sub_data_filtered[[lin_index_col]])) {
          stop("Error: The dependent variable must be numeric.")
        }
        # numeric columns and handling non-numeric columns
        feature_cols <- setdiff(names(sub_data_filtered), lin_index_col)
        #converting non-numeric columns to numeric 
        convert_non_numeric_columns <- function(data) {
          non_numeric_cols <- names(data)[!sapply(data, is.numeric)]
          
          for (col in non_numeric_cols) {
            unique_vals <- unique(data[[col]])
            
            if (length(unique_vals) == 2) {
              #binary columns to 0/1
              message(paste("Converting binary column:", col, "to numeric (0/1)"))
              data[[col]] <- ifelse(data[[col]] == unique_vals[1], 0, 1)
            } else {
              #categorical columns with multiple levels factorize and convert to numeric
              message(paste("Converting categorical column:", col, "to numeric factors"))
              data[[col]] <- as.numeric(as.factor(data[[col]]))
            }
          }
          
          return(data)
        }
        
        #convert data before regression analysis
        sub_data_filtered <- convert_non_numeric_columns(sub_data_filtered)
        
        model_list <- list()
        plot_path_lm <- file.path(all_stat_dir, "linear_reg_plots")
        
        #directory for plots
        if (!dir.exists(plot_path_lm)) {
          dir.create(plot_path_lm, recursive = TRUE, showWarnings = FALSE)
        }
        
        #results data frame
        regression_results <- data.frame()
        print("Starting linear regression analysis")
        
        #loop through features to fit and plot linear models
        for (feature in feature_cols) {
          tryCatch({
            sub_data_filtered <- sub_data_filtered[!is.na(sub_data_filtered[[feature]]), ]
            
            #linear model
            formula <- as.formula(paste(lin_index_col, "~", feature))
            model <- lm(formula, data = sub_data_filtered)
            model_list[[feature]] <- model
            model_name <- paste0("model_", feature)
            
            #plots
            plot_file <- file.path(plot_path_lm, paste0(model_name, ".png"))
            p <- ggplot(sub_data_filtered, aes_string(x = feature, y = lin_index_col)) +
              geom_point() +
              geom_smooth(method = "lm", se = TRUE, color = "blue") +
              labs(
                title = paste("Regression Model for", feature),
                x = feature,
                y = lin_index_col,
                caption = "Regression line with confidence intervals"
              ) +
              theme_minimal()
            ggsave(plot_file, plot = p, width = 8, height = 6)
            
            #residuals
            residuals_plot <- augment(model) %>%
              ggplot(aes(.fitted, .resid)) +
              geom_point() +
              geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
              labs(
                title = paste("Residual Plot for", feature),
                x = "Fitted Values",
                y = "Residuals",
                caption = "Residuals should be randomly scattered around 0"
              ) +
              theme_minimal()
            residuals_file <- file.path(plot_path_lm, paste0(model_name, "_residuals.png"))
            ggsave(residuals_file, plot = residuals_plot, width = 8, height = 6)
            
            #Cooks distance
            cooks_plot <- ggplot(data.frame(CookD = cooks.distance(model)), aes(seq_along(CookD), CookD)) +
              geom_bar(stat = "identity", fill = "red", alpha = 0.7) +
              labs(
                title = paste("Cook's Distance for", feature),
                x = "Observation Index",
                y = "Cook's Distance",
                caption = "High values may indicate influential points"
              ) +
              theme_minimal()
            cooks_file <- file.path(plot_path_lm, paste0(model_name, "_cooks.png"))
            ggsave(cooks_file, plot = cooks_plot, width = 8, height = 6)
            
            #summary metrics
            summary_model <- summary(model)
            summary_row <- data.frame(
              Model = model_name,
              R_squared = summary_model$r.squared,
              Adjusted_R_squared = summary_model$adj.r.squared,
              p_value = summary_model$coefficients[2, 4],
              AIC_value = AIC(model)
            )
            regression_results <- rbind(regression_results, summary_row)
            message(paste("Completed model for feature:", feature))
            
          }, error = function(e) {
            message(paste("Error with feature:", feature, "-", e$message))
          })
        }
        
        return(list(regression_results = regression_results, model_list = model_list))
        
      }, error = function(e) {
        message("Error in perform_simple_regression_analysis:", e$message)
      })
    }
    #RUN
    simp_regression_analysis <- perform_simple_regression_analysis(data, lin_index_col)
    
    if (!is.null(simp_regression_analysis)) {
      simp_regression_results <- simp_regression_analysis$regression_results
      simp_model_list <- simp_regression_analysis$model_list
      
      regression_results_path <- file.path(linear_reg_res_dir, "simple_regression_results.csv")
      tryCatch({
        write.csv(simp_regression_results, regression_results_path, row.names = FALSE)
        print("Simple regression results saved successfully")
      }, error = function(e) {
        message("Failed to save simple regression results:", e$message)
      })
      
      #top 10 models by AIC
      tryCatch({
        top_aic_models_simp <- simp_regression_results %>% arrange(AIC_value) %>% head(10)
        
        plot_path_aic <- file.path(linear_reg_res_dir, "simp-top_aic_models.png")
        ggplot(top_aic_models_simp, aes(x = reorder(Model, AIC_value), y = AIC_value)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          coord_flip() +
          labs(title = "Top AIC Models", x = "Models", y = "AIC") +
          theme_minimal() +
          scale_y_reverse() +
          theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
        ggsave(plot_path_aic)
        
      }, error = function(e) {
        message("Error creating top AIC plot:", e$message)
      })
      
      #closing open graphic devices
      while (!is.null(dev.list())) dev.off()
    }}
  cat("\033[1;32m==== per feature linear regression analysis completed, moving on...  ====\033[0m\n")
  
  ## lin comb try 
  run_combination_linear <- as.logical(tolower(readline(prompt="Do you want to run the combination linear regression analysis? (TRUE/FALSE): ")))
  
  if (run_combination_linear) {
    #index column with validation for linear reg 
    choose_lin_index_column <- function(data) {
      cat("Available columns:\n")
      print(names(data))
      lin_index_col <- readline(prompt = "Enter the name of the index column (including groups to be compared): ")
      if (!(lin_index_col %in% names(data))) {
        stop("The specified index column does not exist.check and try again.")
      }
      return(lin_index_col)
    }
    lin_index_col <- choose_lin_index_column(data)
    
    perform_combination_linear_analysis <- function(data, lin_index_col) {
      sub_data_filtered <- data
      #dependent variable is numeric
      if (!is.numeric(sub_data_filtered[[lin_index_col]])) {
        stop("Error: The dependent variable must be numeric.")
      }
      #numeric columns and handle non-numeric columns
      feature_cols <- setdiff(names(sub_data_filtered), lin_index_col)
      #converting non-numeric columns to numeric format
      convert_non_numeric_columns <- function(data) {
        non_numeric_cols <- names(data)[!sapply(data, is.numeric)]
        
        for (col in non_numeric_cols) {
          unique_vals <- unique(data[[col]])
          
          if (length(unique_vals) == 2) {
            #binary column to 0/1
            message(paste("Converting binary column:", col, "to numeric (0/1)"))
            data[[col]] <- ifelse(data[[col]] == unique_vals[1], 0, 1)
          } else {
            #categorical column with multiple levels factorize and convert to numeric
            message(paste("Converting categorical column:", col, "to numeric factors"))
            data[[col]] <- as.numeric(as.factor(data[[col]]))
          }
        }
        
        return(data)
      }
      #convert data before regression analysis
      sub_data_filtered <- convert_non_numeric_columns(sub_data_filtered)
      #plots path
      plot_path_lm <- file.path(all_stat_dir, "linear_reg_plots_comb")
      if (!dir.exists(plot_path_lm)) {
        dir.create(plot_path_lm, recursive = TRUE)
      }
      
      #parallel processing
      plan(multisession, workers = availableCores() - 2)
      
      regression_results <- data.frame(
        Model = character(),
        AIC_value = numeric(),
        p_value = numeric(),
        R_squared = numeric(),
        Adjusted_R_squared = numeric(),
        Residuals = numeric(),
        Fitted = numeric(),
        Cooks_Distance = numeric(),
        stringsAsFactors = FALSE
      )
      
      print("Running combination linear regression")
      
      chunk_size <- 100 
      
      feature_combinations <- unlist(
        lapply(2:length(feature_cols), function(i) combn(feature_cols, i, simplify = FALSE)),
        recursive = FALSE
      )
      
      #combinations into chunks
      combination_chunks <- split(feature_combinations, ceiling(seq_along(feature_combinations) / chunk_size))
      
      process_chunk <- function(chunk) {
        results <- lapply(chunk, function(combo) {
          tryCatch({
            combo_string <- paste(combo, collapse = "+")
            formula <- as.formula(paste(lin_index_col, "~", combo_string))
            #rows with NA in selected features drop
            select_cols <- c(lin_index_col, combo)
            complete_cases <- sub_data_filtered[complete.cases(sub_data_filtered[, select_cols]), select_cols]
            complete_cases[combo] <- scale(complete_cases[combo])
            if (nrow(complete_cases) == 0) {
              warning(paste("No cases for combination:", combo_string))
              return(NULL)
            }
            
            #model
            model <- lm(formula, data = complete_cases)
            summary_model <- summary(model)
            
            
            list(
              Model = paste(combo, collapse = " + "),
              AIC_value = AIC(model),
              p_value = summary_model$coefficients[2, 4],
              R_squared = summary_model$r.squared,
              Adjusted_R_squared = summary_model$adj.r.squared,
              Residuals = resid(model),
              Fitted = fitted(model),
              Cooks_Distance = cooks.distance(model),
              Complete_Cases = complete_cases
            )
          }, error = function(e) {
            message("Error in linear regression for combination", paste(combo, collapse = ", "), ":", e$message)
            NULL
          })
        })
        return(results[!sapply(results, is.null)])
      }
      
      all_results <- future_map(combination_chunks, process_chunk) %>% unlist(recursive = FALSE)
      
      #sort results
      tryCatch({
        
        
        regression_results <- do.call(rbind, lapply(all_results, function(res) {
          data.frame(
            Model = res$Model,
            AIC_value = res$AIC_value,
            p_value = res$p_value,
            R_squared = res$R_squared,
            Adjusted_R_squared = res$Adjusted_R_squared
          )
        })) %>% distinct()
        
        top_50 <- regression_results %>% arrange(AIC_value, p_value) %>% head(50)
        #top 50 models
        linear_results_path <- file.path(linear_reg_res_dir, "combination_linear_results.csv")
        write.csv(top_50, linear_results_path, row.names = FALSE)
        print("Top 50 Combination linear regression results saved.")
        
        #top models
        top_20 <- top_50 %>% head(20)
        
        #plot top 20 models by AIC
        plot_path_aic <- file.path(linear_reg_res_dir, "top_comb_linear_aic_models.png")
        top_20_aic_plot <- ggplot(top_20, aes(x = reorder(Model, -AIC_value), y = AIC_value)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          coord_flip() +
          labs(
            title = "Top Combination Models by AIC",
            x = "Models",
            y = "AIC Value"
          ) +
          theme_minimal()+
          scale_y_reverse() 
        
        ggsave(plot_path_aic, plot = top_20_aic_plot, width = 48, height = 26)
        while (!is.null(dev.list())) dev.off()
        
        ##metadata 
        metadata <- data.frame(
          Model = character(),
          Hash = character(),
          stringsAsFactors = FALSE
        )
        
        
        #plotting for Cooks Distance, residuals, and regression plots
        lapply(1:nrow(top_20), function(i) {
          model_row <- top_20[i, ]
          model_result <- all_results[[which(sapply(all_results, function(x) x$Model == model_row$Model))]]
          complete_cases <- model_result$Complete_Cases
          cooks_distances <- model_result$Cooks_Distance
          fitted_values <- model_result$Fitted
          residuals <- model_result$Residuals
          #hash for the model name
          model_hash <- digest::digest(model_row$Model, algo = "crc32")
          
          #regression
          combo <- strsplit(model_row$Model, " \\+ ")[[1]]
          dependent_var <- lin_index_col
          independent_vars <- combo[1]
          
          p_reg <- ggplot(complete_cases, aes_string(x = independent_vars, y = dependent_var)) +
            geom_point(alpha = 0.6) +
            geom_smooth(method = "lm", se = TRUE, color = "blue") +
            labs(
              title = paste("Regression Model:", model_row$Model),
              x = paste(independent_vars),
              y = dependent_var,
              caption = "Regression line with confidence intervals"
            ) +
            theme_minimal()
          
          model_hash <- digest::digest(model_row$Model, algo = "crc32")
          ggsave(file.path(plot_path_lm, paste0("regression_plot_", model_hash, ".png")), plot = p_reg, width = 8, height = 6)
          
          #cooks distance 
          cooks_df <- data.frame(
            Observation = seq_along(cooks_distances),
            Cooks_Distance = cooks_distances
          )
          
          
          p_cooks <- ggplot(cooks_df, aes(x = Observation, y = Cooks_Distance)) +
            geom_bar(stat = "identity", fill = "red", alpha = 0.7) +
            labs(
              title = paste(model_row$Model),
              x = "Observation",
              y = "Cook's Distance",
              caption = "High values may indicate influential points"
            ) +
            theme_minimal()
          
          ggsave(file.path(plot_path_lm, paste0("cooks_distance_plot_", model_hash, ".png")), plot = p_cooks, width = 8, height = 6)
          
          #residuals
          residuals_df <- data.frame(
            Fitted = fitted_values,
            Residuals = residuals
          )
          
          p_residual <- ggplot(residuals_df, aes(x = Fitted, y = Residuals)) +
            geom_point(alpha = 0.6) +
            geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
            labs(
              title = paste(model_row$Model),
              x = "Fitted Values",
              y = "Residuals",
              caption = "Residuals should be randomly scattered around 0"
            ) +
            theme_minimal()
          
          ggsave(file.path(plot_path_lm, paste0("residual_plot_", model_hash, ".png")), plot = p_residual, width = 8, height = 6)
          
          #metadata add
          metadata <<- rbind(metadata, data.frame(Model = model_row$Model, Hash = model_hash))
        })
      }, error = function(e) {
        message("Error in logistic plotting for combination: ", e$message)
        return(NULL)
      })
      # metadata to csv
      metadata_path <- file.path(plot_path_lm, "model_metadata.csv")
      write.csv(metadata, metadata_path, row.names = FALSE)
      print(paste("Model metadata saved to:", metadata_path))
      
      return(regression_results)
    }
    
    #RUN
    combination_linear_analysis <- perform_combination_linear_analysis(data, lin_index_col)
    
  }
  cat("\033[1;32m==== combination feature linear regression analysis completed, moving on...  ====\033[0m\n")
  
  cat("\033[1;32m==== Choose the column that represents the groups you would like to compare in the rest of the analysis ====\033[0m\n")
  #select an index column with validation
  choose_index_column <- function(data) {
    cat("Available columns:\n")
    print(names(data))
    index_col <- readline(prompt = "Enter the name of the index column (including groups to be compared): ")
    if (!(index_col %in% names(data))) {
      stop("The specified index column does not exist.check and try again.")
    }
    return(index_col)
  }
  index_col <- choose_index_column(data)
  # unique values in index column to verify
  if (index_col %in% colnames(data)) {
    print(unique(data[[index_col]]))
  } else {
    stop("Index column could not be found in data.")
  }
  #separate numeric columns, retaining the index column
  numeric_data <- data[, sapply(data, is.numeric)]
  data <- cbind(data[[index_col]], numeric_data)
  colnames(data)[1] <- index_col
  norm_function <- function(data, index_col) {
    tryCatch({
      original_data <- data
      apply_norm <- readline(prompt = "Do you want to normalize your data? (yes/no): ")
      if (tolower(apply_norm) == "yes") {
        cat("Normalizing data using CLR.\n")
        feature_data <- data[, !(names(data) %in% index_col)]
        data <- logratio.transfo(feature_data, logratio = 'CLR')
      }
      else {
        cat("No normalization performed.\n")
        return(original_data)
      }
      return(data)
    }, error = function(e) {
      cat("\033[1;31mAn error occurred during normalization:\033[0m", e$message, "\n")
      cat("Returning original data and skipping this step.\n")
      return(original_data)
    })
  }
  data <- norm_function(data, index_col)
  write.csv(data, file.path(output_dir, "post_norm_data.csv"), row.names = FALSE)
  cat("\033[1;32m==== normalization step completed, moving on... ====\033[0m\n")
  cli_alert_success("Data processing complete. Proceeding to statistical analysis.")
  cat("\033[1;31m==== Step 2 - multi group testing began... ====\033[0m\n")
  ## outlier segment 
  outlier_detection <- function(data, index_col, threshold = 3) {
    tryCatch({
      dummy_frame <- data
      original_data <- data
      apply_outlier <- readline(prompt = "Would you like to apply outlier detection? (yes/no): ")
      
      if (tolower(apply_outlier) == "yes") {
        method <- readline(prompt = "Choose the outlier detection method (mad/iqr): ")
        threshold <- as.numeric(readline(prompt = "Enter the threshold for outlier removal (default is 3): "))
        if (is.na(threshold)) threshold <- 3
        
        outlier_dir <- file.path(output_dir, "outliers_results")
        dir.create(outlier_dir, showWarnings = FALSE)
        
        feature_data <- data[, !(names(data) %in% index_col)]
        outliers <- data.frame()
        log_summary <- data.frame(Feature = character(), Num_Outliers = numeric(), stringsAsFactors = FALSE)
        
        mad_outliers <- function(x) {
          mad_val <- mad(x, constant = 1.4826)
          lower_bound <- median(x) - threshold * mad_val
          upper_bound <- median(x) + threshold * mad_val
          which(x < lower_bound | x > upper_bound)
        }
        
        iqr_outliers <- function(x) {
          iqr_val <- IQR(x, na.rm = TRUE)
          lower_bound <- quantile(x, 0.25, na.rm = TRUE) - threshold * iqr_val
          upper_bound <- quantile(x, 0.75, na.rm = TRUE) + threshold * iqr_val
          which(x < lower_bound | x > upper_bound)
        }
        
        detect_outliers <- switch(tolower(method),
                                  "mad" = mad_outliers,
                                  "iqr" = iqr_outliers,
                                  stop("Invalid method chosen. Please choose either 'mad' or 'iqr'."))
        
        for (feature in colnames(feature_data)) {
          outlier_rows <- detect_outliers(feature_data[[feature]])
          
          if (length(outlier_rows) > 0) {
            feature_outliers <- dummy_frame[outlier_rows, ]
            feature_outliers$Feature <- feature
            feature_outliers$Original_Row_Index <- row.names(dummy_frame)[outlier_rows]
            outliers <- rbind(outliers, feature_outliers)
            
            data[outlier_rows, feature] <- NA
            log_summary <- rbind(log_summary, data.frame(Feature = feature, Num_Outliers = length(outlier_rows)))
            cat(paste0("Detected ", length(outlier_rows), " outliers in feature '", feature, "'\n"))
          }
        }
        
        if (nrow(outliers) > 0) {
          write.csv(outliers, file.path(outlier_dir, "removed_outliers.csv"), row.names = TRUE)
          write.csv(log_summary, file.path(outlier_dir, "outlier_summary_log.csv"), row.names = FALSE)
          cat("Outliers and summary log have been saved in the 'outliers_results' directory.\n")
          
          review_outliers <- readline(prompt = "Do you want to proceed with these adjustments (yes to confirm, no to revert)? (yes/no): ")
          if (tolower(review_outliers) == "no") {
            cat("Reverting changes and returning original data.\n")
            return(original_data)  # Original data without outlier adjustments
          } else {
            remove_rows <- readline(prompt = "Do you want to remove rows with NA values? (yes/no): ")
            if (tolower(remove_rows) == "yes") {
              data <- data[complete.cases(data), ]
            } else {
              remove_cols <- readline(prompt = "Do you want to remove columns with NA values? (yes/no): ")
              if (tolower(remove_cols) == "yes") {
                data <- data[, colSums(is.na(data)) == 0]
              }
            }
          }
        } else {
          cat("No outliers detected.\n")
        }
      } else {
        cat("Outlier detection not applied.\n")
      }
      
      return(data)
    }, error = function(e) {
      cat("\033[1;31mAn error occurred during outlier detection:\033[0m", e$message, "\n")
      cat("Returning original data and skipping this step.\n")
      return(original_data)
    })
  }
  
  data <- outlier_detection(data, index_col)
  write.csv(data, file.path(output_dir, "full_data_with_out.csv"), row.names = FALSE)
  cat("\033[1;32m==== Outlier detection completed, moving on... ====\033[0m\n")
  #######
  perform_multi_reg_all <- readline(prompt = "Would you like to run multi_reg on all groups present in the index column? (yes/no): ")
  if (tolower(perform_multi_reg_all) == "yes") {
    multi_reg_dir <- file.path(all_stat_dir, "multi_reg_dir")
    dir.create(multi_reg_dir, showWarnings = FALSE)
    
    perform_multi_regression <- function(data, index_col, dependent_vars, results_dir) {
      cat("Starting multinomial regression with user-defined reference group.\n")
      
      #choose the reference group
      cat("Available groups:\n")
      unique_groups <- unique(data[[index_col]])
      for (i in seq_along(unique_groups)) {
        cat(i, ": ", unique_groups[i], "\n")
      }
      reference_index <- as.numeric(readline(prompt = "Enter the index of the reference group (by number): "))
      if (is.na(reference_index) || reference_index > length(unique_groups)) {
        stop("Invalid reference group index provided.")
      }
      reference_group <- unique_groups[reference_index]
      cat("Selected reference group: ", reference_group, "\n")
      
      #reference group setup
      data[[index_col]] <- factor(data[[index_col]], levels = c(reference_group, setdiff(unique_groups, reference_group)))
      
      # data prep
      tryCatch({
        sub_data <- data[, c(index_col, dependent_vars)]
        cat("Data preparation complete.\n")
      }, error = function(e) {
        stop("Error preparing data for regression: ", conditionMessage(e))
      })
      
      #results dir
      if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
      
      # multinomial model
      tryCatch({
        formula <- as.formula(paste(index_col, "~ ."))
        multinom_model <- multinom(formula, data = sub_data)
        model_summary <- summary(multinom_model)
        cat("Model fitting complete.\n")
      }, error = function(e) {
        stop("Error fitting multinomial model: ", conditionMessage(e))
      })
      
      #model summary
      tryCatch({
        summary_path <- file.path(results_dir, "multinomial_model_summary.txt")
        sink(summary_path)
        print(model_summary)
        sink()
        cat("Model summary saved.\n")
      }, error = function(e) {
        stop("Error saving model summary: ", conditionMessage(e))
      })
      
      #LRT for significance
      cat("Performing Likelihood Ratio Tests (LRT).\n")
      lrt_results <- list()
      for (feature in dependent_vars) {
        reduced_formula <- as.formula(paste(index_col, "~ . -", feature))
        reduced_model <- multinom(reduced_formula, data = sub_data)
        lrt <- anova(reduced_model, multinom_model, test = "Chisq")
        lrt_results[[feature]] <- lrt
      }
      write.csv(do.call(rbind, lrt_results), file.path(results_dir, "lrt_results.csv"), row.names = TRUE)
      cat("LRT results saved.\n")
      
      #CV
      cat("Performing cross-validation.\n")
      folds <- createFolds(sub_data[[index_col]], k = 3, list = TRUE, returnTrain = FALSE)
      cv_accuracy <- sapply(folds, function(test_indices) {
        train_data <- sub_data[-test_indices, ]
        test_data <- sub_data[test_indices, ]
        cv_model <- multinom(formula, data = train_data)
        predicted <- predict(cv_model, newdata = test_data)
        mean(predicted == test_data[[index_col]])
      })
      
      cat("Fold sizes:\n")
      for (i in seq_along(folds)) {
        train_data_p <- sub_data[-folds[[i]], ]
        test_data_p <- sub_data[folds[[i]], ]
        cat(sprintf("Fold %d - Training size: %d, Test size: %d\n", i, nrow(train_data_p), nrow(test_data_p)))
      }
      group_length <- length(unique(sub_data[[index_col]]))
      minimum_samples_threshold <- 10 * group_length
      if (nrow(test_data_p) < minimum_samples_threshold) {
        warning(sprintf("Test set size in Fold %d is too small for reliable evaluation.", i))
      }
      
      avg_cv_accuracy <- mean(cv_accuracy)
      write.csv(data.frame(Fold = seq_along(cv_accuracy), Accuracy = cv_accuracy), file.path(results_dir, "cv_results.csv"), row.names = FALSE)
      cat("Cross-validation complete. Average Accuracy: ", round(avg_cv_accuracy * 100, 2), "%\n")
      
      #predictions and accuracy calculation
      tryCatch({
        predicted_groups <- predict(multinom_model, newdata = sub_data)
        accuracy <- mean(predicted_groups == sub_data[[index_col]])
        cat("Model Accuracy:", round(accuracy * 100, 2), "%\n")
      }, error = function(e) {
        warning("Error in calculating predictions or accuracy: ", conditionMessage(e))
        accuracy <- NA  # Set accuracy to NA if error occurs
      })
      
      #confusion matrix and plot
      tryCatch({
        conf_matrix <- confusionMatrix(predicted_groups, sub_data[[index_col]])
        conf_matrix_path <- file.path(results_dir, "confusion_matrix.txt")
        sink(conf_matrix_path)
        print(conf_matrix)
        sink()
        
        conf_matrix_plot <- as.data.frame(conf_matrix$table)
        colnames(conf_matrix_plot) <- c("True Class", "Predicted Class", "Count")
        ggplot(conf_matrix_plot, aes(x = `True Class`, y = `Predicted Class`, fill = Count)) +
          geom_tile() +
          geom_text(aes(label = Count)) +
          scale_fill_gradient(low = "white", high = "blue") +
          ggtitle("Confusion Matrix") +
          theme_minimal()
        ggsave(file.path(results_dir, "confusion_matrix_plot.png"))
        cat("Confusion matrix and plot saved.\n")
      }, error = function(e) {
        warning("Error creating or saving confusion matrix and plot: ", conditionMessage(e))
      })
      
      #coefficients, odds ratios & feature importance plotting
      tryCatch({
        #coefficients and calculate odds ratios
        raw_coefs <- coef(multinom_model)
        coef_multinom <- exp(raw_coefs)  # Odds ratios
        #standard errors for Wald tests
        se <- sqrt(diag(vcov(multinom_model)))
        
        #wald test z-scores and p-values
        z_scores <- raw_coefs / se
        p_values <- 2 * (1 - pnorm(abs(z_scores)))  # Two-tailed p-values
        
        #results data frame
        results_df <- as.data.frame(as.table(coef_multinom))
        colnames(results_df) <- c("Group", "Feature", "Odds Ratio")
        results_df$Coefficient <- as.vector(raw_coefs)
        results_df$Z_Score <- as.vector(z_scores)
        results_df$P_Value <- as.vector(p_values)
        results_df$Significance <- ifelse(results_df$P_Value < 0.05, "Significant", "Not Significant")
        
        write.csv(results_df, file.path(results_dir, "model_coef_results.csv"), row.names = FALSE)
        cat("Detailed model coefficients and significance testing results saved.\n")
      }, error = function(e) {
        warning("Error generating or saving model coefficients or feature importance plot: ", conditionMessage(e))
      })
      
      
      tryCatch({
        #feature importance as the mean absolute value of coefficients across groups
        coef_abs <- abs(coef(multinom_model)) 
        feature_names <- colnames(coef_abs)  
        importance_values <- colMeans(coef_abs)
        
        importance_df <- data.frame(Feature = feature_names, Importance = importance_values)
        importance_df <- importance_df[order(-importance_df$Importance), ]  # Sort by importance
        
        #importance plot
        ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          coord_flip() +
          ggtitle("Feature Importance") +
          xlab("Feature") + ylab("Average Absolute Coefficient") +
          theme_minimal()
        
        ggsave(file.path(results_dir, "feature_importance_plot.png"))
        cat("Feature importance plot saved.\n")
      }, error = function(e) {
        warning("Error generating or saving model coefficients or feature importance plot: ", conditionMessage(e))
      })
      
      
      return(list(summary = model_summary, accuracy = accuracy, coefficients = coef_multinom))
    }
    
    dependent_vars <- colnames(data)[-1]
    results <- perform_multi_regression(data, index_col, dependent_vars, multi_reg_dir)
  }
  
  cat("\033[1;32m==== Multi regression completed, moving on... ====\033[0m\n")
  ###############
  #box plots for each feature against the index column
  perform_bp_all <- readline(prompt = "Would you like to plot box plots on all groups present in the index column? (yes/no): ")
  if (tolower(perform_bp_all) == "yes") {
    
    plot_box_plots <- function(data, index_col, output_dir) {
      boxplot_dir <- file.path(output_dir, "boxplots")
      dir.create(boxplot_dir, showWarnings = FALSE)
      
      #user to specify group order
      cat("Available groups:\n")
      unique_groups <- unique(data[[index_col]])
      for (i in seq_along(unique_groups)) {
        cat(i, ": ", unique_groups[i], "\n")
      }
      group_order_indices <- as.numeric(strsplit(readline(prompt="Enter the indices of the groups in the desired order, separated by a comma: "), ",")[[1]])
      if (any(is.na(group_order_indices)) || any(group_order_indices > length(unique_groups))) {
        cat("Invalid group order input. Proceeding with the default order.\n")
        ordered_groups <- unique_groups
      } else {
        ordered_groups <- unique_groups[group_order_indices]
      }
      
      data[[index_col]] <- factor(data[[index_col]], levels = ordered_groups)
      
      for (feature in colnames(data)[-1]) {
        p <- ggplot(data, aes_string(x = index_col, y = feature)) +
          geom_boxplot(aes(fill = !!sym(index_col))) +
          geom_jitter(width = 0.2, size = 0.5) +
          labs(title = paste("Box Plot for", feature), x = index_col, y = feature) +
          theme(axis.text.x = element_text(size = 8)) +
          scale_fill_brewer(palette = "Set3")
        
        boxplot_path <- file.path(boxplot_dir, paste("boxplot_", gsub(" ", "_", feature), ".png", sep = ""))
        ggsave(boxplot_path, plot = p, width = 8, height = 6)
        boxplot_path_2 <- file.path(boxplot_dir, paste("boxplot_", gsub(" ", "_", feature), ".pdf", sep = ""))
        ggsave(boxplot_path_2, plot = p, width = 8, height = 6, device = "pdf")
      }
    }
    
    plot_box_plots(data, index_col, output_dir)
  }
  
  cat("\033[1;32m==== Box plots of all features generated, moving on... ====\033[0m\n")
  #multi group statistical testing
  perform_stat_tests <- readline(prompt = "Would you like to run statistical tests on all groups present in the index column? (yes/no): ")
  if (tolower(perform_stat_tests) == "yes") {
    multi_group_test_dir <- file.path(all_stat_dir, "multi_group_test")
    dir.create(multi_group_test_dir, showWarnings = FALSE)
    #choosing statistical test
    cat("Choose a statistical test to run:\n1: ANOVA (assuming normality) & Tukey\n2: Kruskal-Wallis (non-parametric) & Dunn's test\n")
    test_choice <- suppressWarnings(as.numeric(readline(prompt = "Enter the number of the test you want to perform: ")))
    
    if (is.na(test_choice) || !(test_choice %in% c(1, 2))) {
      cat("Invalid choice. Please select 1 or 2 for the test type.\n")
    } else {
      #result lists
      results_list <- list()
      posthoc_results_list <- list()
      
      #run over features, running the selected tests
      for (feature in names(data)[-which(names(data) == index_col)]) {
        cat("\nProcessing feature:", feature, "\n")
        
        tryCatch({
          if (test_choice == 1) {
            # ANOVA test
            cat("Running ANOVA...\n")
            aov_results <- aov(data[[feature]] ~ data[[index_col]])
            aov_summary <- summary(aov_results)
            results_list[[feature]] <- aov_summary
            
            # Post-hoc test (Tukey HSD)
            tukey_results <- TukeyHSD(aov_results)
            posthoc_results_list[[feature]] <- tukey_results
            cat("ANOVA and Tukey HSD completed for", feature, "\n")
            
          } else if (test_choice == 2) {
            # Kruskal-Wallis test
            cat("Running Kruskal-Wallis...\n")
            kruskal_results <- kruskal.test(data[[feature]] ~ data[[index_col]])
            results_list[[feature]] <- kruskal_results
            
            # Post-hoc test (Dunn's test)
            dunn_results <- dunnTest(data[[feature]] ~ data[[index_col]], method = "bonferroni")
            posthoc_results_list[[feature]] <- dunn_results$res
            cat("Kruskal-Wallis and Dunn's test completed for", feature, "\n")
          }
          
        }, error = function(e) {
          cat("Error in processing feature:", feature, "\n", e$message, "\n")
        })
      }
      
      #save results if tests were successful
      if (length(results_list) > 0) {
        if (test_choice == 1) {
          #ANOVA results
          aov_results_df <- do.call(rbind, lapply(names(results_list), function(feature) {
            cbind(Feature = feature, as.data.frame(results_list[[feature]][[1]]))
          }))
          write.csv(aov_results_df, file.path(multi_group_test_dir, "anova_results.csv"), row.names = FALSE)
          
          #Tukey results
          tukey_results_df <- do.call(rbind, lapply(names(posthoc_results_list), function(feature) {
            cbind(Feature = feature, as.data.frame(posthoc_results_list[[feature]][[1]]))
          }))
          write.csv(tukey_results_df, file.path(multi_group_test_dir, "tukey_hsd_results.csv"), row.names = TRUE)
          
        } else if (test_choice == 2) {
          #Kruskal] results
          kruskal_results_df <- do.call(rbind, lapply(names(results_list), function(feature) {
            data.frame(Feature = feature, Statistic = results_list[[feature]]$statistic, P_value = results_list[[feature]]$p.value)
          }))
          write.csv(kruskal_results_df, file.path(multi_group_test_dir, "kruskal_wallis_results.csv"), row.names = FALSE)
          
          # Dunn results
          dunn_results_df <- do.call(rbind, lapply(names(posthoc_results_list), function(feature) {
            cbind(Feature = feature, posthoc_results_list[[feature]])
          }))
          write.csv(dunn_results_df, file.path(multi_group_test_dir, "dunn_test_results.csv"), row.names = FALSE)
        }
        
        cat("Statistical test results saved in:", multi_group_test_dir, "\n")
      } else {
        cat("No valid results to save; check data or feature selection.\n")
      }
    }
  }
  cat("\033[1;32m==== Multi group statistical test completed, moving on... ====\033[0m\n")
  perform_manova <- readline(prompt = "Would you like to perform MANOVA on all groups present in the index column? (yes/no): ")
  if (tolower(perform_manova) == "yes") {
    manova_test_dir <- file.path(all_stat_dir, "manova_test")
    dir.create(manova_test_dir, showWarnings = FALSE)
    
    perform_manova <- function(data, index_col, dependent_vars, test_statistic = "Pillai") {
      cat("Preparing data for MANOVA...\n")
      
      #verify index column 
      if (!index_col %in% colnames(data)) {
        cat("Error: Specified index column does not exist in the data.\n")
        return(NULL)
      }
      
      #at least two unique groups in the index column
      unique_groups <- length(unique(data[[index_col]]))
      if (unique_groups < 2) {
        cat("Error: MANOVA requires at least two unique groups in the index column.\n")
        return(NULL)
      }
      
      # sufficient dependent variables
      if (length(dependent_vars) < 2) {
        cat("Error: MANOVA requires at least two dependent variables. Please provide additional numeric features.\n")
        return(NULL)
      }
      
      tryCatch({
        sub_data <- data[, c(index_col, dependent_vars)]
        sub_data[[index_col]] <- as.factor(sub_data[[index_col]])
      }, error = function(e) {
        cat("Error preparing data for MANOVA: ", e$message, "\n")
        return(NULL)
      })
      
      #dependent variables to a numeric matrix
      feat_matrix <- tryCatch({
        as.matrix(sub_data[, dependent_vars])
      }, error = function(e) {
        cat("Error converting dependent variables to matrix: ", e$message, "\n")
        return(NULL)
      })
      
      cat("Running MANOVA...\n")
      #MANOVA run
      manova_results <- tryCatch({
        manova_result <- manova(feat_matrix ~ sub_data[[index_col]])
        
        #MANOVA results
        manova_summary <- tryCatch({
          summary(manova_result, test = test_statistic)
        }, error = function(e) {
          cat("Error summarizing MANOVA results: ", e$message, "\n")
          return(NULL)
        })
        
        if (!is.null(manova_summary)) {
          overall_p_value <- manova_summary$stats[1, "Pr(>F)"]
          
          #ANOVA for individual feature testing
          anova_results <- tryCatch({
            summary.aov(manova_result)
          }, error = function(e) {
            cat("Error performing ANOVA for individual features: ", e$message, "\n")
            return(NULL)
          })
          
          #significant features if ANOVA succeeds
          if (!is.null(anova_results)) {
            significant_genes_df <- data.frame(Gene = character(), P.Value = numeric(), Feature = character(), stringsAsFactors = FALSE)
            
            for (i in names(anova_results)) {
              p_values <- anova_results[[i]]$`Pr(>F)`
              feature_names <- rownames(anova_results[[i]])
              significant <- which(p_values < 0.1)  # Threshold for significance
              
              if (length(significant) > 0) {
                temp_df <- data.frame(
                  Gene = feature_names[significant],
                  P.Value = p_values[significant],
                  Feature = i,
                  stringsAsFactors = FALSE
                )
                significant_genes_df <- rbind(significant_genes_df, temp_df)
              }
            }
            
            write.csv(significant_genes_df, file.path(manova_test_dir, "anova_significant_results.csv"), row.names = FALSE)
          }
          
          manova_summary_df <- data.frame(Statistic = test_statistic, P.Value = overall_p_value)
          write.csv(manova_summary_df, file.path(manova_test_dir, "manova_summary.csv"), row.names = FALSE)
          
          return(list(ManovaSummary = manova_summary, AnovaResults = anova_results))
        }
      }, error = function(e) {
        cat("Error performing MANOVA: ", e$message, "\n")
        return(NULL)
      })
      
      return(manova_results)
    }
    
    #Run
    dependent_vars <- colnames(data)[-1]  
    ## thresh for features 
    var_threshold <- 0.4  
    high_variance_vars <- names(which(apply(data[, dependent_vars], 2, var) > var_threshold))
    dependent_vars <- intersect(dependent_vars, high_variance_vars)
    
    manova_results <- perform_manova(data, index_col, dependent_vars)
  }
  cat("\033[1;32m==== MANOVA analysis completed, moving on... ====\033[0m\n")
  ## boot all groups
  boot_dat_all <- data
  all_groups <- unique(boot_dat_all[[index_col]])
  
  perform_boot_all <- readline(prompt = "Would you like to run bootstrapping analysis tests on all groups present in the index column? (yes/no): ")
  if (tolower(perform_boot_all) == "yes") {
    
    #groups to test from the index column for all boot
    choose_groups_wide <- function(boot_dat_all, index_col) {
      unique_groups <- unique(boot_dat_all[[index_col]])
      cat("Available groups:\n")
      for (i in seq_along(unique_groups)) {
        cat(i, ": ", unique_groups[i], "\n")
      }
      group_indices <- as.numeric(strsplit(readline(prompt="Enter the indices of the groups to test for bootstrap confidence interval plots (order by plot needs), separated by a comma: "), ",")[[1]])
      
      if (any(is.na(group_indices)) || any(group_indices > length(unique_groups))) {
        cat("Invalid group selection. Please ensure the indices correspond to the available groups.\n")
        return(NULL)
      }
      chosen_groups <- unique_groups[group_indices]
      return(chosen_groups)
    }
    
    chosen_groups_wide <- choose_groups_wide(boot_dat_all, index_col)
    
    if (is.null(chosen_groups_wide)) {
      cat("No valid groups selected. Skipping bootstrap analysis.\n")
    } else {
      #plot bootstrap confidence intervals
      plot_bootstrap_ci <- function(boot_dat_all, index_col, output_dir) {
        ci_plot_dir <- file.path(output_dir, "boot_ci_plots_all")
        dir.create(ci_plot_dir, showWarnings = FALSE)
        
        for (feature in colnames(boot_dat_all)[-1]) {
          #safeguard against NA values in the current feature
          feature_data <- boot_dat_all[!is.na(boot_dat_all[[feature]]), ]
          if (nrow(feature_data) == 0) {
            cat(sprintf("Skipping feature '%s' as it contains only NA values.\n", feature))
            next
          }
          
          tryCatch({
            #compute bootstrap confidence intervals
            dabest_data <- load(data = feature_data, x = !!sym(index_col), y = !!sym(feature), idx = chosen_groups_wide)
            dabest_data_diff <- mean_diff(dabest_data)
            
            #create and save plot
            ci_plot <- dabest_plot(dabest_data_diff, rawplot.type = "swarmplot")
            ci_plot_path <- file.path(ci_plot_dir, paste("boot_ci_plot_all_", gsub(" ", "_", feature), ".png", sep = ""))
            ggsave(ci_plot_path, plot = ci_plot, width = 8, height = 6)
            ci_plot_path_2 <- file.path(ci_plot_dir, paste("boot_ci_plot_all_", gsub(" ", "_", feature), ".pdf", sep = ""))
            ggsave(ci_plot_path_2, plot = ci_plot, width = 8, height = 6, device = "pdf")
            cat(sprintf("Saved bootstrap CI plot for feature '%s' at '%s'\n", feature, ci_plot_path))
          }, error = function(e) {
            cat(sprintf("Error generating bootstrap plot for feature '%s': %s\n", feature, e$message))
            error_log_path <- file.path(ci_plot_dir, "bootstrap_error_log.txt")
            write(sprintf("Feature: %s | Error: %s\n", feature, e$message), file = error_log_path, append = TRUE)
          })
        }
      }
      
      #RUN
      plot_bootstrap_ci(boot_dat_all, index_col, output_dir)
    }
  }
  
  cat("\033[1;32m==== Bootstrap analysis completed, moving on... ====\033[0m\n")
  cli_alert_success("multi group statistical testing handled successfully!")
  cat("\033[1;31m==== Step 3 - 2 group statistical testing!  ====\033[0m\n")
  #choose groups for 2 groups comparision
  choose_groups <- function(data, index_col) {
    unique_groups <- unique(data[[index_col]])
    cat("Available groups:\n")
    for (i in seq_along(unique_groups)) {
      cat(i, ": ", unique_groups[i], "\n")
    }
    group_input <- readline(prompt="Enter the indices of 2 groups to test (comma-separated) or type 'done' to exit: ")
    
    if (tolower(group_input) == "done") {
      return(NULL)  # Signal to exit
    }
    
    group_indices <- as.numeric(strsplit(group_input, ",")[[1]])
    
    if (length(group_indices) == 2 && all(group_indices %in% seq_along(unique_groups))) {
      return(unique_groups[group_indices])
    } else {
      cat("\033[1;31mInvalid selection. Please enter two valid indices.\033[0m\n")
    }
  }
  
  # Loop for multiple 2-group comparisons
  repeat {
    chosen_groups <- choose_groups(data, index_col)
    if (is.null(chosen_groups)) {
      cat("\033[1;32mExiting two-group comparisons. Moving to the next analysis step...\033[0m\n")
      break
    }
    
    #sub-data based on chosen groups
    sub_data <- data %>% filter(!!sym(index_col) %in% chosen_groups)
    sub_data_dir <- file.path(output_dir, paste("sub_data_", paste(chosen_groups, collapse = "_"), sep = ""))
    dir.create(sub_data_dir, showWarnings = FALSE)
    write.csv(sub_data, file.path(sub_data_dir, "sub_data.csv"), row.names = FALSE)
    cat("\033[1;32m==== groups chosen succesfully, moving on.... ====\033[0m\n")
    ## DABSETR BOOT PLOTS per 2 groups
    boot_dat_two <- sub_data
    all_groups <- unique(boot_dat_two[[index_col]])
    
    perform_boot_two <- readline(prompt = "Would you like to run bootstrapping analysis tests on two groups present in the index column? (yes/no): ")
    if (tolower(perform_boot_two) == "yes") {
      
      #plot bootstrap confidence intervals
      plot_bootstrap_ci <- function(boot_dat_two, index_col, output_dir) {
        ci_plot_dir <- file.path(sub_data_dir, "boot_ci_plots")
        dir.create(ci_plot_dir, showWarnings = FALSE)
        
        for (feature in colnames(boot_dat_two)[-1]) {
          #safeguard against NA values in the current feature
          feature_data <- boot_dat_two[!is.na(boot_dat_two[[feature]]), ]
          if (nrow(feature_data) == 0) {
            cat(sprintf("Skipping feature '%s' as it contains only NA values.\n", feature))
            next
          }
          
          tryCatch({
            #compute bootstrap confidence intervals
            dabest_data <- load(data = feature_data, x = !!sym(index_col), y = !!sym(feature), idx = chosen_groups)
            dabest_data_diff <- mean_diff(dabest_data)
            
            #create and save the plot
            ci_plot <- dabest_plot(dabest_data_diff, rawplot.type = "swarmplot")
            ci_plot_path <- file.path(ci_plot_dir, paste("boot_ci_plot_", gsub(" ", "_", feature), ".png", sep = ""))
            ggsave(ci_plot_path, plot = ci_plot, width = 8, height = 6)
            ci_plot_path_2 <- file.path(ci_plot_dir, paste("boot_ci_plot_", gsub(" ", "_", feature), ".pdf", sep = ""))
            ggsave(ci_plot_path_2, plot = ci_plot, width = 8, height = 6, device = "pdf")
            cat(sprintf("Saved bootstrap CI plot for feature '%s' at '%s'\n", feature, ci_plot_path))
          }, error = function(e) {
            cat(sprintf("Error generating bootstrap plot for feature '%s': %s\n", feature, e$message))
            error_log_path <- file.path(ci_plot_dir, "bootstrap_error_log.txt")
            write(sprintf("Feature: %s | Error: %s\n", feature, e$message), file = error_log_path, append = TRUE)
          })
        }
      }
      
      #RUN
      plot_bootstrap_ci(boot_dat_all, index_col, output_dir)
    }
    
    stats_info_dir <- file.path(sub_data_dir, "stats_info")
    dir.create(stats_info_dir, showWarnings = FALSE)
    cat("\033[1;32m==== Performed bootstrap analysis on 2 chosen groups succesfully, moving on.... ====\033[0m\n")
    #levene for equal variance 
    print("Testing for equal variance among chosen groups")
    perform_variance_testing <- function(data, index_col, chosen_groups) {
      variance_results <- data.frame(feature = character(), p_value = numeric(), stringsAsFactors = FALSE)
      #validate chosen_groups
      if (length(chosen_groups) < 2) {
        stop("Error: At least two groups are required for variance testing.")
      }
      sub_data_filtered <- data %>% filter(!!sym(index_col) %in% chosen_groups)
      
      for (feature in colnames(data)[-1]) {
        formula <- as.formula(paste(feature, "~", index_col))
        
        tryCatch({
          #Levene test to assess equal variance
          levene_test_result <- leveneTest(formula, data = sub_data_filtered)
          variance_results <- rbind(variance_results,
                                    data.frame(feature = feature,
                                               p_value = levene_test_result$`Pr(>F)`[1]))
          cat("Levene's test for", feature, "completed successfully.\n")
          
        }, error = function(e) {
          cat("Error in Levene's test for feature", feature, ":", e$message, "\n")
        })
      }
      
      return(variance_results)
    }
    variance_results <- perform_variance_testing(data, index_col, chosen_groups)
    
    #variance test results if they exist
    if (nrow(variance_results) > 0) {
      variance_results_path <- file.path(stats_info_dir, "variance_test_results.csv")
      write.csv(variance_results, variance_results_path, row.names = FALSE)
      cat("Variance test results saved to:", variance_results_path, "\n")
    } else {
      cat("No variance test results to save; check input data or group selection.\n")
    }
    
    #shapiro-wilk normality test 
    print("Testing for normality among chosen groups")
    
    perform_normality_testing <- function(data, index_col, chosen_groups) {
      normality_results <- data.frame(feature = character(), p_value_group1 = numeric(), p_value_group2 = numeric(), stringsAsFactors = FALSE)
      
      #ensure exactly two groups for pairwise comparison
      if (length(chosen_groups) != 2) {
        stop("Error: Normality testing requires exactly two groups.")
      }
      
      sub_data_filtered <- data %>% filter(!!sym(index_col) %in% chosen_groups)
      
      for (feature in colnames(data)[-1]) {
        tryCatch({
          #shapiro-wilk test 
          group1_data <- sub_data_filtered[sub_data_filtered[[index_col]] == chosen_groups[1], feature]
          group2_data <- sub_data_filtered[sub_data_filtered[[index_col]] == chosen_groups[2], feature]
          
          if (length(group1_data) < 3 || length(group2_data) < 3) {
            cat("Insufficient data for normality test on feature:", feature, "\n")
            next
          }
          
          shapiro_group1 <- shapiro.test(group1_data)
          shapiro_group2 <- shapiro.test(group2_data)
          
          normality_results <- rbind(normality_results,
                                     data.frame(feature = feature,
                                                p_value_group1 = shapiro_group1$p.value,
                                                p_value_group2 = shapiro_group2$p.value))
          cat("Shapiro-Wilk test for", feature, "completed successfully.\n")
          
        }, error = function(e) {
          cat("Error in Shapiro-Wilk test for feature", feature, ":", e$message, "\n")
        })
      }
      
      return(normality_results)
    }
    
    normality_results <- perform_normality_testing(data, index_col, chosen_groups)
    
    #normality test results if they exist
    if (nrow(normality_results) > 0) {
      normality_results_path <- file.path(stats_info_dir, "normality_test_results.csv")
      write.csv(normality_results, normality_results_path, row.names = FALSE)
      cat("Normality test results saved to:", normality_results_path, "\n")
    } else {
      cat("No normality test results to save; check input data or group selection.\n")
    }
    
    cat("\033[1;32m==== Performed variance and normality analysis on 2 chosen groups succesfully, moving on.... ====\033[0m\n")
    stat_test_dir <- file.path(sub_data_dir, "stat_test")
    dir.create(stat_test_dir, showWarnings = FALSE)
    #2 group stat test t-test, wilcoxon test or Mann-Whitney U test
    perform_statistical_testing <- function(data, index_col, chosen_groups) {
      #user to choose of the test options
      cat("Choose the test to perform based on normality results:\n")
      cat("1: t-test (assuming normality)\n")
      cat("2: Wilcoxon test (non-parametric, not assuming normality)\n")
      cat("3: Mann-Whitney U test (non-parametric, compares medians)\n")
      test_choice <- as.numeric(readline(prompt="Enter the number of the chosen test: "))
      
      #validate test selection
      if (!test_choice %in% c(1, 2, 3)) {
        stop("Invalid choice for the test. Please enter 1, 2, or 3.")
      }
      #index column is treated as a factor and subset to the two chosen groups
      data[[index_col]] <- factor(data[[index_col]], levels = chosen_groups)
      sub_data_filtered <- data %>% filter(!!sym(index_col) %in% chosen_groups)
      #filtered data has exactly two groups
      if (nlevels(sub_data_filtered[[index_col]]) != 2) {
        stop("The selected index column must have exactly two groups for this test.")
      }
      
      #results data frame
      test_results <- data.frame(feature = character(), p_value = numeric(), stringsAsFactors = FALSE)
      
      for (feature in colnames(data)[-1]) {
        cat("Processing feature:", feature, "\n")
        
        #group data
        group1_data <- sub_data_filtered[sub_data_filtered[[index_col]] == chosen_groups[1], feature]
        group2_data <- sub_data_filtered[sub_data_filtered[[index_col]] == chosen_groups[2], feature]
        
        #is data sufficient for testing
        if (length(group1_data) < 2 || length(group2_data) < 2) {
          cat("Skipping feature due to insufficient data in one or both groups.\n")
          next
        }
        
        tryCatch({
          if (test_choice == 1) {
            # in case of T - Levene test on the filtered subset for equal variances
            levene_test <- car::leveneTest(as.formula(paste(feature, "~", index_col)), data = sub_data_filtered)
            
            if (levene_test$`Pr(>F)`[1] < 0.05) {
              cat("Variances are unequal; performing Welch's t-test for", feature, "\n")
              test_result <- t.test(group1_data, group2_data, var.equal = FALSE)  # Welch's t-test
            } else {
              cat("Variances are equal; performing regular t-test for", feature, "\n")
              test_result <- t.test(group1_data, group2_data, var.equal = TRUE)  # Regular t-test
            }
            
          } else if (test_choice == 2) {
            # Wilcoxon test
            cat("Performing Wilcoxon test for", feature, "\n")
            test_result <- wilcox.test(group1_data, group2_data)
            
          } else if (test_choice == 3) {
            # Mann-Whitney U test
            cat("Performing Mann-Whitney U test for", feature, "\n")
            test_result <- wilcox.test(group1_data, group2_data, exact = FALSE)  # Mann-Whitney
            
          }
          
          #result to data frame
          test_results <- rbind(test_results, data.frame(feature = feature, p_value = test_result$p.value))
          
        }, warning = function(w) {
          cat("Warning for feature", feature, ":", conditionMessage(w), "\n")
        }, error = function(e) {
          cat("Error in statistical testing for feature", feature, ":", conditionMessage(e), "\n")
        })
      }
      #check if results were obtained
      if (nrow(test_results) == 0) {
        warning("No test results were generated. Verify data and chosen groups.")
      } else {
        cat("Statistical testing completed successfully.\n")
      }
      
      return(test_results)
    }
    
    #RUN
    test_results <- perform_statistical_testing(data, index_col, chosen_groups)
    cat("\033[1;32m==== Performed 2 group statistical test succesfully, moving on.... ====\033[0m\n")
    #apply FDR correction on p vals?
    tryCatch({
      apply_fdr <- as.logical(readline(prompt="Do you want to apply FDR correction to the p-values? (TRUE/FALSE): "))
      
      if (apply_fdr) {
        test_results$p_value_adj <- p.adjust(test_results$p_value, method = "fdr")
        p_value_histogram <- ggplot(test_results, aes(x = p_value_adj)) +
          geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
          labs(title = "Histogram of P-values-adj", x = "P-value-adj", y = "Frequency")
        
        histogram_path <- file.path(stat_test_dir, "p_value_adj_histogram.png")
        ggsave(histogram_path, plot = p_value_histogram, width = 8, height = 6)
        print("fdr applied- see adj p value histogram in results directory")
        
      }})
    
    #save statistical test results
    if (nrow(test_results) > 0) {
      test_results_path <- file.path(stat_test_dir, "statistical_test_results.csv")
      write.csv(test_results, test_results_path, row.names = FALSE)
      cat("2 group test results saved to:", test_results_path, "\n")
    } else {
      cat("No 2 group test results to save; check input data or group selection.\n")
    }
    
    #histogram of p-values
    tryCatch({
      p_value_histogram <- ggplot(test_results, aes(x = p_value)) +
        geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
        labs(title = "Histogram of P-values", x = "P-value", y = "Frequency")
      histogram_path <- file.path(stat_test_dir, "p_value_histogram.png")
      ggsave(histogram_path, plot = p_value_histogram, width = 8, height = 6)
      print("2 group statistical test applied- see  p value histogram in results directory")
    })
    cat("\033[1;32m==== 2 group statistical test, FDR, and p value hostograms completed, moving on.... ====\033[0m\n")
    summary_dir <- file.path(sub_data_dir, "stat_summary")
    dir.create(summary_dir, showWarnings = FALSE)
    #summarize DF 
    summarize_dataframe <- function(df, index_col) {
      # Validate input data
      if (!index_col %in% names(df)) {
        stop("Error: index_col does not exist in the DataFrame.")
      }
      summary_df <- data.frame(
        column_name = character(),
        mean = numeric(),
        variance = numeric(),
        num_samples = integer(),
        group_means = numeric(),
        group_variances = numeric(),
        group_num_samples = integer(),
        total_mean = numeric(),
        total_variance = numeric(),
        total_num_samples = integer(),
        calculated_effect_size = numeric(),
        S = numeric(),
        stringsAsFactors = FALSE
      )
      
      #summary statistics
      for (col_name in names(df)) {
        if (col_name != index_col) {
          tryCatch({
            #mean, variance, and number of samples for column
            col_mean <- mean(df[[col_name]], na.rm = TRUE)
            col_var <- var(df[[col_name]], na.rm = TRUE)
            col_num_samples <- sum(!is.na(df[[col_name]]))
            
            #mean, variance, and number of samples for groups
            group_means <- tapply(df[[col_name]], df[[index_col]], mean, na.rm = TRUE)
            group_variances <- tapply(df[[col_name]], df[[index_col]], var, na.rm = TRUE)
            group_num_samples <- tapply(df[[col_name]], df[[index_col]], function(x) sum(!is.na(x)))
            
            #two groups for comparison
            if (length(group_means) < 2) {
              warning(paste("Column", col_name, "does not have at least two groups. Skipping..."))
              next
            }
            
            #effect size and S value
            mean_dif <- group_means[1] - group_means[2]
            s_value <- sqrt(
              ((group_num_samples[1] - 1) * group_variances[1] +
                 (group_num_samples[2] - 1) * group_variances[2]) /
                (col_num_samples - 2)
            )
            
            #CIs
            v_value <- ((group_num_samples[1] - 1) * group_variances[1] +
                          (group_num_samples[2] - 1) * group_variances[2]) /
              (col_num_samples - 2)
            marg <- qt(0.975, df = col_num_samples - 1) * sqrt(v_value / group_num_samples[1] + v_value / group_num_samples[2])
            low_int <- mean_dif - marg
            up_int <- mean_dif + marg
            
            #append to summary DF
            summary_col <- data.frame(
              column_name = col_name,
              mean = col_mean,
              variance = col_var,
              num_samples = col_num_samples,
              group_means = group_means,
              group_variances = group_variances,
              group_num_samples = group_num_samples,
              total_mean = col_mean,
              total_variance = col_var,
              total_num_samples = col_num_samples,
              calculated_effect_size = abs(mean_dif) / s_value,
              S = s_value,
              V = v_value,
              mean_dif = mean_dif,
              marg = marg,
              low_int = low_int,
              up_int = up_int,
              group_sd = sqrt(group_variances)
            )
            
            summary_df <- rbind(summary_df, summary_col)
            
          }, error = function(e) {
            message(paste("Error in column", col_name, ":", e$message))
          })
        }
      }
      return(summary_df)
    }
    
    #RUN
    summary_df <- summarize_dataframe(sub_data, index_col)
    
    #save summary DF
    summary_df_path <- file.path(summary_dir, "summary_data.csv")
    tryCatch({
      write.csv(summary_df, summary_df_path, row.names = FALSE)
      message("Summary data saved successfully at ", summary_df_path)
    }, error = function(e) {
      message("Failed to save summary data: ", e$message)
    })
    
    #plot CIs
    two_row_ci_plot <- function(df, path) {
      plot_list <- list()
      n <- nrow(df)
      
      for (i in seq(1, n, by = 2)) {
        tryCatch({
          df_small <- df[i, ]
          up_int <- df_small$up_int
          low_int <- df_small$low_int
          effect_size <- df_small$mean_dif
          column_name <- df_small$column_name[1]
          
          plot <- ggplot(df_small, aes(x = 1, y = effect_size)) +
            geom_point(aes(y = mean_dif), shape = 21, colour = "black", fill = "white", size = 5, stroke = 5) +
            geom_errorbar(aes(ymin = low_int, ymax = up_int), width = 0.2) +
            geom_hline(yintercept = effect_size, linetype = "dotted") +
            labs(title = column_name) +
            theme_minimal()
          
          plot_list[[i]] <- plot
          
        }, error = function(e) {
          message("Error generating plot for row ", i, ": ", e$message)
        })
      }
      
      #save each plot and close open devices
      for (i in seq_along(plot_list)) {
        tryCatch({
          ggsave(paste0(path, "/plot_", i, ".png"), plot_list[[i]], width = 8, height = 6)
          if (length(dev.list()) > 10) {
            while (!is.null(dev.list())) dev.off()
          }
        }, error = function(e) {
          message("Error saving plot ", i, ": ", e$message)
        })
      }
    }
    
    # RUN
    ci_plot_dir <- file.path(summary_dir, "ci_plots")
    dir.create(ci_plot_dir, showWarnings = FALSE)
    two_row_ci_plot(summary_df, ci_plot_dir)
    cat("\033[1;32m==== 2 group statistical summary & Cohen's CI plots completed, moving on.... ====\033[0m\n")
    #function to add power calculations to the summary data frame
    add_power_to_df <- function(df) {
      n1 <- as.numeric(readline(prompt="Enter the value for n1: "))
      n2 <- as.numeric(readline(prompt="Enter the value for n2: "))
      sig.level <- as.numeric(readline(prompt="Enter the significance level: "))
      
      power_vals <- rep(NA, nrow(df))
      
      #run on each row and compute the power using pwr.t2n.test
      for (i in 1:nrow(df)) {
        d <- df$calculated_effect_size[i]
        power <- pwr.t2n.test(n1 = n1, n2 = n2, d = d, sig.level = sig.level)$power
        power_vals[i] <- power
      }
      
      df$power <- power_vals
      
      return(df)
    }
    
    summary_df_with_power <- add_power_to_df(summary_df)
    
    #save summary data frame with power calculations
    summary_df_with_power_path <- file.path(summary_dir, "summary_data_with_power.csv")
    write.csv(as.data.frame(summary_df_with_power), summary_df_with_power_path)
    cat("\033[1;32m==== power analysis completed and saved, moving on.... ====\033[0m\n")
    bayes_dir <- file.path(sub_data_dir, "stat_bayes")
    dir.create(bayes_dir, showWarnings = FALSE)
    ## Bayes factor addition
    perform_bayes_factor_analysis <- function(data, index_col, chosen_groups) {
      split_dataframe <- function(df) {
        index_col <- names(df)[1]
        df_col_names <- names(df)[-1]
        split_df_list <- list()
        for (col_name in df_col_names) {
          feat_df <- df[!is.na(df[[col_name]]), ]
          new_df <- feat_df[c(index_col, col_name)]
          split_df_list[[col_name]] <- new_df
        }
        return(split_df_list)
      }
      
      split_dfs <- split_dataframe(data)
      
      ttest_bf <- function(df) {
        group1 <- df[[2]][df[[1]] == chosen_groups[1]]
        group2 <- df[[2]][df[[1]] == chosen_groups[2]]
        result <- ttestBF(x = group1, y = group2, rscale = "medium")
        return(result)
      }
      
      bayes_results <- data.frame(feature = character(), bayes_factor = numeric(), stringsAsFactors = FALSE)
      error_log <- list()
      
      for (df_name in names(split_dfs)) {
        tryCatch({
          res <- ttest_bf(split_dfs[[df_name]])
          res_df <- as.data.frame(res)
          bayes_results <- rbind(bayes_results, data.frame(feature = df_name, bayes_factor = res_df$bf[1]))
        }, error = function(e) {
          warning(paste("Error processing feature:", df_name, "-", conditionMessage(e)))
          error_log[[df_name]] <- conditionMessage(e)
        })
      }
      
      # Log errors
      if (length(error_log) > 0) {
        error_log_path <- file.path(bayes_dir, "bayes_factor_errors.log")
        writeLines(paste(names(error_log), unlist(error_log), sep = ": "), con = error_log_path)
        cat("Errors logged for Bayes factor analysis.\n")
      }
      
      return(bayes_results)
    }
    
    #prep Bayes factor analysis
    bayes_dat <- data
    tryCatch({
      bayes_results <- perform_bayes_factor_analysis(bayes_dat, index_col, chosen_groups)
      #Bayes factor results
      bayes_results_path <- file.path(bayes_dir, "bayes_factor_results.csv")
      write.csv(bayes_results, bayes_results_path, row.names = FALSE)
      cat("\033[1;32m==== Bayes factor analysis completed and saved, moving on.... ====\033[0m\n")
    }, error = function(e) {
      stop("Critical error in Bayes factor analysis: ", conditionMessage(e))
    })
    #permutation testing 
    perform_permutation_testing <- function(data, index_col, chosen_groups, num_permutations) {
      perm_test_results <- data.frame(feature = character(), p_value = numeric(), stringsAsFactors = FALSE)
      
      #validate chosen_groups
      if (length(chosen_groups) != 2) {
        stop("Error: Permutation testing requires exactly two groups.")
      }
      
      #filter data 
      sub_data_filtered <- data %>% filter(!!sym(index_col) %in% chosen_groups)
      
      for (feature in colnames(data)[-1]) {
        tryCatch({
          group1_data <- sub_data_filtered[sub_data_filtered[[index_col]] == chosen_groups[1], feature]
          group2_data <- sub_data_filtered[sub_data_filtered[[index_col]] == chosen_groups[2], feature]
          
          #data is sufficient 
          if (length(group1_data) < 3 || length(group2_data) < 3) {
            cat("Insufficient data for permutation test on feature:", feature, "\n")
            next
          }
          
          original_diff <- abs(mean(group1_data) - mean(group2_data))
          
          #permutation testing
          perm_diffs <- numeric(num_permutations)
          combined_data <- c(group1_data, group2_data)
          
          for (i in 1:num_permutations) {
            permuted_data <- sample(combined_data)
            perm_group1 <- permuted_data[1:length(group1_data)]
            perm_group2 <- permuted_data[(length(group1_data) + 1):length(permuted_data)]
            perm_diffs[i] <- abs(mean(perm_group1) - mean(perm_group2))
          }
          
          perm_p_value <- mean(perm_diffs >= original_diff)
          perm_test_results <- rbind(perm_test_results, data.frame(feature = feature, p_value = perm_p_value))
          
          cat("Permutation test for", feature, "completed successfully.\n")
          
        }, error = function(e) {
          cat("Error in permutation test for feature", feature, ":", e$message, "\n")
        })
      }
      
      return(perm_test_results)
    }
    
    #input for permutation testing
    perform_permutations <- as.logical(tolower(readline(prompt="Do you want to perform permutation tests? (TRUE/FALSE): ")))
    
    if (perform_permutations) {
      num_permutations <- as.numeric(readline(prompt="Enter the number of permutations (e.g., 500, 1000, or 10000): "))
      
      # num_permutations
      if (is.na(num_permutations) || num_permutations <= 0) {
        stop("Error: Number of permutations must be a positive integer.")
      }
      
      #RUN
      perm_test_results <- perform_permutation_testing(data, index_col, chosen_groups, num_permutations)
      
      if (nrow(perm_test_results) > 0) {
        perm_test_results_path <- file.path(stat_test_dir, "permutation_test_results.csv")
        write.csv(perm_test_results, perm_test_results_path, row.names = FALSE)
        cat("Permutation test results saved to:", perm_test_results_path, "\n")
      } else {
        cat("No permutation test results to save; check input data or group selection.\n")
      }
    }
    cat("\033[1;32m==== Permuations analysis completed and saved, moving on.... ====\033[0m\n")
    #box plots for sub-data
    plot_box_plots <- function(data, index_col, output_dir) {
      boxplot_dir <- file.path(output_dir, "boxplots")
      dir.create(boxplot_dir, showWarnings = FALSE)
      
      for (feature in colnames(data)[-1]) {
        p <- ggplot(data, aes_string(x=index_col, y=feature)) +
          geom_boxplot(aes(fill = !!sym(index_col))) +
          geom_jitter(width = 0.2, size = 0.5) +
          labs(title = paste("Box Plot for", feature), x = index_col, y = feature) +
          theme(axis.text.x = element_text(size = 8)) +
          scale_fill_brewer(palette = "Set3")
        
        boxplot_path <- file.path(boxplot_dir, paste("boxplot_", gsub(" ", "_", feature), ".png", sep=""))
        ggsave(boxplot_path, plot = p, width = 8, height = 6)
        boxplot_path_2 <- file.path(boxplot_dir, paste("boxplot_", gsub(" ", "_", feature), ".pdf", sep=""))
        ggsave(boxplot_path_2, plot = p, width = 8, height = 6, device = "pdf")
      }}
    plot_box_plots(sub_data, index_col, sub_data_dir)
    cat("\033[1;32m==== Box plots for two groups completed and saved, moving on.... ====\033[0m\n")
    cat("\033[1;32m==== Moving to Regression analysis... ====\033[0m\n")
    cat("\033[1;31m==== Highly recommended to inspect correlations between features prior...  ====\033[0m\n")
    log_reg_res_dir <- file.path(sub_data_dir, "stat_log_regression")
    dir.create(log_reg_res_dir, showWarnings = FALSE)
    #prompt user to run simple logistic regression
    run_simple_logistic <- as.logical(tolower(readline(prompt="Do you want to run per feature logistic regression analysis? (TRUE/FALSE): ")))
    
    if (run_simple_logistic) {
      
      perform_simple_logistic_analysis <- function(data, index_col, chosen_groups) {
        # validate chosen groups
        if (length(chosen_groups) != 2) {
          stop("Error: Logistic regression requires exactly two groups.")
        }
        sub_data_filtered <- data %>% filter(!!sym(index_col) %in% chosen_groups)
        
        #index col in factor form
        sub_data_filtered[, 1] <- as.factor(sub_data_filtered[, 1])
        feature_cols <- names(sub_data_filtered)[-1]
        group_levels <- levels(factor(sub_data_filtered[[index_col]]))
        group_0 <- group_levels[1]
        group_1 <- group_levels[2]
        print(group_levels)
        model_list <- list()
        logistic_results <- data.frame(Model = character(), AIC_value = numeric(), p_value = numeric(), AUC = numeric(),
                                       Odds_Ratio = numeric(),
                                       Lower_CI = numeric(),
                                       Upper_CI = numeric(), stringsAsFactors = FALSE)
        plot_path_glm <- file.path(sub_data_dir, "logistic_reg_plots")
        if (!dir.exists(plot_path_glm)) {
          dir.create(plot_path_glm, recursive = TRUE)
        }
        
        print("Running logistic regression")
        for (feature in feature_cols) {
          tryCatch({
            sub_data_filtered <- sub_data_filtered[!is.na(sub_data_filtered[[feature]]), ]
            sub_data_filtered[[feature]] <- scale(sub_data_filtered[[feature]])
            #logistic regression model
            model <- glm(sub_data_filtered[[1]] ~ sub_data_filtered[[feature]], family = binomial)
            
            #model convergence
            if (!model$converged) {
              warning(paste("Model for feature", feature, "did not converge. Skipping."))
              next
            }
            #model to list
            model_name <- paste0("logistic_model_", feature)
            model_list[[model_name]] <- model
            plot_data <- augment(model)
            plot_data <- as.data.frame(plot_data)
            y_obs_axis_title <- paste0(" Observed Values (0 = ", group_0, ", 1 = ", group_1, ")")
            y_pred_axis_title <- paste0("Predicted Values (0 = ", group_0, ", 1 = ", group_1, ")")
            sub_data_filtered[[1]] <- type.convert(factor(sub_data_filtered[[1]], labels = 0:1), as.is = TRUE)
            sub_data_filtered$.fitted <- predict(model, type = "response")
            
            #AUC-ROC computation & plotting
            roc_curve <- roc(sub_data_filtered[[1]], fitted(model))
            auc_value <- auc(roc_curve)
            roc_plot <- ggroc(roc_curve) +
              labs(title = paste("ROC Curve for", feature),
                   x = "False Positive Rate",
                   y = "True Positive Rate") +
              theme_minimal()
            roc_file <- file.path(plot_path_glm, paste0(model_name, "_ROC.png"))
            ggsave(roc_file, plot = roc_plot, width = 8, height = 6)
            
            #observed values with sigmoidal curve plot
            plot_observed <- ggplot(sub_data_filtered, aes(x = sub_data_filtered[[feature]], y = sub_data_filtered[[1]])) +
              geom_jitter(color = "red", alpha = 0.6, width = 0.02, height = 0.002) +  
              stat_smooth(
                method = "glm",
                method.args = list(family = "binomial"),
                color = "green",
                se = TRUE, fullrange = TRUE, n = 1000, size = 1.5, linetype = "dashed"
              ) +  
              xlim(min(sub_data_filtered[[feature]]) - 1, max(sub_data_filtered[[feature]]) + 1) +
              labs(
                title = paste("Observed Values vs Feature for", feature),
                x = feature,
                y = y_obs_axis_title
              ) +
              scale_y_continuous(breaks = c(0, 1), labels = c(group_0, group_1)) +
              theme_minimal()
            
            #predicted probabilities with sigmoidal curve plot
            plot_predicted <- ggplot(sub_data_filtered, aes_string(x = sub_data_filtered[[feature]], y =  sub_data_filtered[[".fitted"]])) +
              geom_point(color = "black", alpha = 0.6) +  
              stat_smooth(
                method = "glm",
                method.args = list(family = "binomial"),
                color = "blue",
                se = TRUE, fullrange = TRUE, n = 1000, size = 0.5, linetype = "dashed") +
              xlim(min(sub_data_filtered[[feature]]) - 1, max(sub_data_filtered[[feature]]) + 1) +
              labs(title = paste("Predicted Probabilities vs Feature for", feature),
                   x = feature,
                   y = y_pred_axis_title) +
              theme_minimal()
            
            #saving plots
            ggsave(file.path(plot_path_glm, paste0("observed_values_", feature, ".png")), plot = plot_observed, width = 8, height = 6)
            ggsave(file.path(plot_path_glm, paste0("predicted_probabilities_", feature, ".png")), plot = plot_predicted, width = 8, height = 6)
            
            #CIs and odds ratios
            conf_intervals <- confint(model)
            odds_ratios <- exp(coef(model))
            odds_ratios_conf <- exp(conf_intervals)
            
            #model summary 
            summary_model <- summary(model)
            summary_row <- data.frame(Model = model_name,
                                      AIC_value = AIC(model),
                                      p_value = coef(summary_model)[2, "Pr(>|z|)"],
                                      AUC = auc_value,
                                      Odds_Ratio = odds_ratios[2],
                                      Lower_CI = odds_ratios_conf[2, 1],
                                      Upper_CI = odds_ratios_conf[2, 2])
            logistic_results <- rbind(logistic_results, summary_row)
            
          }, error = function(e) {
            cat("Error in logistic regression for feature", feature, ":", e$message, "\n")
          })
        }
        
        return(list(logistic_results = logistic_results, model_list = model_list))
      }
      
      #logistic regression and results
      simple_logistic_analysis <- perform_simple_logistic_analysis(data, index_col, chosen_groups)
      simple_logistic_results <- simple_logistic_analysis$logistic_results
      
      #save simple logistic regression results
      logistic_results_path <- file.path(log_reg_res_dir, "simple_logistic_results.csv")
      write.csv(simple_logistic_results, logistic_results_path, row.names = FALSE)
      print("Simple logistic regression results saved.")
      
      #sort by AIC and select top 10 models
      top_10_aic_models <- simple_logistic_results %>% arrange(AIC_value) %>% head(10)
      
      #plot top 10 models by AIC, with lower AIC values at the top
      plot_path_aic <- file.path(log_reg_res_dir, "top_logistic_aic_models.png")
      ggplot(top_10_aic_models, aes(x = reorder(Model, AIC_value), y = AIC_value)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        labs(title = "Top Logistic Regression Models by AIC",
             x = "Model",
             y = "AIC") +
        theme_minimal() +
        scale_y_reverse() + 
        theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
      
      ggsave(plot_path_aic, width = 8, height = 6)
      while (!is.null(dev.list())) dev.off()
    }
    cat("\033[1;32m==== per feature logistics regression analysis completed, moving on...  ====\033[0m\n")
    #user to run combination logistic regression
    run_combination_logistic <- as.logical(tolower(readline(prompt="Do you want to run the combination logistic regression analysis? (TRUE/FALSE): ")))
    
    if (run_combination_logistic) {
      
      perform_combination_logistic_analysis <- function(data, index_col, chosen_groups) {
        #validate chosen_groups
        if (length(chosen_groups) != 2) {
          stop("Error: Combination logistic regression requires exactly two groups.")
        }
        
        sub_data_filtered <- data %>% filter(!!sym(index_col) %in% chosen_groups)
        sub_data_filtered[[index_col]] <- as.factor(sub_data_filtered[[index_col]])
        feature_cols <- names(sub_data_filtered)[-1]
        
        #parallel processing
        plan(multisession, workers = availableCores() - 2)
        
        #logistic results storage
        logistic_results <- data.frame(
          Model = character(),
          AIC_value = numeric(),
          p_value = numeric(),
          AUC = numeric(),
          Odds_Ratio = numeric(),
          Lower_CI = numeric(),
          Upper_CI = numeric(),
          stringsAsFactors = FALSE
        )
        
        # plots dir
        plot_path_glm <- file.path(sub_data_dir, "log_reg_comb_plots")
        if (!dir.exists(plot_path_glm)) {
          dir.create(plot_path_glm, recursive = TRUE)
        }
        print("Running combination logistic regression")
        
        chunk_size <- 100  #adjust if needed
        
        feature_combinations <- unlist(
          lapply(2:length(feature_cols), function(i) combn(feature_cols, i, simplify = FALSE)),
          recursive = FALSE
        )
        
        #combinations into chunks
        combination_chunks <- split(feature_combinations, ceiling(seq_along(feature_combinations) / chunk_size))
        
        process_chunk <- function(chunk) {
          results <- lapply(chunk, function(combo) {
            tryCatch({
              combo_string <- paste(combo, collapse = "+")
              formula <- as.formula(paste(index_col, "~", combo_string))
              
              #rows with NA in selected features drop
              select_cols <- c(index_col, combo)
              complete_cases <- sub_data_filtered[complete.cases(sub_data_filtered[, select_cols]), select_cols]
              complete_cases[combo] <- scale(complete_cases[combo])
              if (nrow(complete_cases) == 0) {
                warning(paste("No cases for combination:", combo_string))
                return(NULL)
              }
              
              #model
              model <- glm(formula, data = complete_cases, family = binomial)
              if (!model$converged) {
                warning(paste("Model for combination", combo_string, "did not converge. Skipping."))
                return(NULL)
              }
              
              #AIC, p-value, and CIs
              conf_intervals <- confint(model)
              odds_ratios <- exp(coef(model))
              odds_ratios_conf <- exp(conf_intervals)
              summary_model <- summary(model)
              
              auc_value <- auc(roc(complete_cases[[index_col]], fitted(model)))
              
              list(
                Model = paste(combo, collapse = " + "),
                AIC_value = AIC(model),
                p_value = coef(summary_model)[2, "Pr(>|z|)"],
                AUC = auc_value,
                Odds_Ratio = odds_ratios[2],
                Lower_CI = odds_ratios_conf[2, 1],
                Upper_CI = odds_ratios_conf[2, 2],
                Model_Object = model,
                Complete_Cases = complete_cases
              )
            }, error = function(e) {
              message("Error in logistic regression for combination", paste(combo, collapse = ", "), ":", e$message)
              NULL
            })
          })
          return(results[!sapply(results, is.null)])
        }
        
        all_results <- future_map(combination_chunks, process_chunk) %>% unlist(recursive = FALSE)
        
        #sort results
        tryCatch({
          logistic_results <- do.call(rbind, lapply(all_results, function(res) {
            data.frame(
              Model = res$Model,
              AIC_value = res$AIC_value,
              p_value = res$p_value,
              AUC = res$AUC,
              Odds_Ratio = res$Odds_Ratio,
              Lower_CI = res$Lower_CI,
              Upper_CI = res$Upper_CI
            )
          }))
          
          top_50 <- logistic_results %>% arrange(AIC_value, p_value) %>% head(50)
          #top 50 models to csv
          logistic_results_path <- file.path(log_reg_res_dir, "combination_logistic_results.csv")
          write.csv(top_50, logistic_results_path, row.names = FALSE)
          print("Top 50 Combination logistic regression results saved.")
          
          #top 20 AIC models
          top_20 <- top_50 %>% head(20)
          
          #plot top 10 models by AIC, with lower AIC values at the top
          plot_path_aic <- file.path(log_reg_res_dir, "top_20_comb_logistic_aic_models.png")
          top_20_aic_plot <- ggplot(top_20, aes(x = reorder(Model, -AIC_value), y = AIC_value)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            coord_flip() +
            labs(
              title = "Top Logistic Combination Models by AIC",
              x = "Model",
              y = "AIC Value"
            ) +
            theme_minimal()+
            scale_y_reverse() 
          
          ggsave(plot_path_aic, plot = top_20_aic_plot, width = 48, height = 26)
          while (!is.null(dev.list())) dev.off()
          
          ##metadata 
          metadata <- data.frame(
            Model = character(),
            Hash = character(),
            stringsAsFactors = FALSE
          )
          
          #metadata dir
          metadata_path <- file.path(plot_path_glm, "model_metadata.csv")
          
          #observed, predicted and roc plots for top 20 models
          lapply(1:nrow(top_20), function(i) {
            model_row <- top_20[i, ]
            model_obj <- all_results[[which(sapply(all_results, function(x) x$Model == model_row$Model))]]$Model_Object
            complete_cases <- all_results[[which(sapply(all_results, function(x) x$Model == model_row$Model))]]$Complete_Cases
            group_levels <- levels(factor(complete_cases[[index_col]]))
            group_0 <- group_levels[1]
            group_1 <- group_levels[2]
            print(group_levels)
            
            y_obs_axis_title <- paste0("Observed Values (0 = ", group_0, ", 1 = ", group_1, ")")
            y_pred_axis_title <- paste0("Predicted Values (0 = ", group_0, ", 1 = ", group_1, ")")
            complete_cases[[index_col]] <- type.convert(factor(complete_cases[[index_col]], labels = 0:1), as.is = TRUE)
            complete_cases$.fitted <- predict(model_obj, type = "response")
            
            
            roc_obj <- roc(complete_cases[[index_col]], predict(all_results[[which(sapply(all_results, function(x) x$Model == model_row$Model))]]$Model_Object, type = "response"))
            roc_plot <- ggroc(roc_obj) +
              ggtitle(paste("ROC for Model:", model_row$Model)) +
              theme_minimal()
            
            
            plot_observed <- ggplot(complete_cases, aes(x = complete_cases[[2]], y = complete_cases[[index_col]])) +
              geom_jitter(color = "red", alpha = 0.6, width = 0.02, height = 0.002) +
              stat_smooth(
                method = "glm",
                method.args = list(family = "binomial"),
                color = "green",
                se = TRUE, fullrange = TRUE, n = 1000, size = 1.5, linetype = "dashed"
              ) +
              labs(
                title = paste("Observed vs Feature for", model_row$Model),
                x = paste(model_row$Model),
                y = y_obs_axis_title
              ) +
              scale_y_continuous(breaks = c(0, 1), labels = c(group_0, group_1)) +  
              theme_minimal()
            
            
            plot_predicted <- ggplot(complete_cases, aes(x = complete_cases[[2]], y = complete_cases$.fitted)) +
              geom_point(color = "black", alpha = 0.6) +
              stat_smooth(
                method = "glm",
                method.args = list(family = "binomial"),
                color = "blue",
                se = TRUE, fullrange = TRUE, n = 1000, size = 0.5, linetype = "dashed"
              ) +
              labs(
                title = paste("Predicted vs Feature for", model_row$Model),
                x = paste(model_row$Model),
                y = y_pred_axis_title
              ) +
              theme_minimal()
            
            
            #hash for the models to fix long names 
            model_hash <- digest::digest(model_row$Model, algo = "crc32")
            #metadata add
            metadata <<- rbind(metadata, data.frame(Model = model_row$Model, Hash = model_hash))
            roc_plot_path <- file.path(plot_path_glm, paste0(model_hash, "_ROC.png"))
            ggsave(roc_plot_path, plot = roc_plot, width = 16, height = 12)
            ggsave(file.path(plot_path_glm, paste0("observed_values_", model_hash, ".png")), plot = plot_observed, width = 8, height = 6)
            ggsave(file.path(plot_path_glm, paste0("predicted_values_", model_hash, ".png")), plot = plot_predicted, width = 8, height = 6)
          })
        }, error = function(e) {
          message("Error in logistic plotting for combination: ", e$message)
          return(NULL)
        })
        # metadata to csv
        write.csv(metadata, metadata_path, row.names = FALSE)
        print(paste("Model metadata saved to:", metadata_path))
        
        return(logistic_results)
      }
      
      #RUN combination logistic regression and results
      combination_logistic_analysis <- perform_combination_logistic_analysis(data, index_col, chosen_groups)
      
    }
    cat("\033[1;32m==== combination feature logistics regression analysis completed, moving on...  ====\033[0m\n")
    
    cli_alert_success(" two group statistical testing handled successfully!")
  }
  cat("\033[1;31m==== Analysis completed, review files and run again as needed:)  ====\033[0m\n")
  cat("\033[1;31m==== IMPORTANT! PLEASE REMOVE YOUR FILES FROM THE DIRECTORY ====\033[0m\n")
}
perform_statistical_analysis()
