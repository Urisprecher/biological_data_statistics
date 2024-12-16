### Batch effect analysis
## takes as input a folder that includes an additional folder inside it with meta data and feature data.
## all output files will be saved in the main folder including all steps in the analysis.

##run these lines  
detect_batch_effects <- function(
                                 output_folder_prompt = TRUE,
                                 rowname_prompt = TRUE,
                                 batch_param_prompt = TRUE,
                                 num_features = 20
) {
  

  
  cat("\033[1;32m==== Starting Batch effect analysis!====\033[0m\n")
  message("Please ensure files are correctly formatted and in a folder inside the main directory.")
  message("loading function and libraries, please wait...")
  
  #required libraries
  required_packages <- c("ggplot2", "pheatmap", "vegan", "clr", "RColorBrewer", "cluster", "pheatmap", "mixOmics", "ggplot2",
                         "MASS", "Matrix", "knitr", "xtable", "sva", "gridExtra", "limma", "variancePartition",
                         "pvca", "ruv", "lmerTest", "bapred", "readxl", "cli")


  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      if (!require(pkg, character.only = TRUE)) {
        stop(paste("Package", pkg, "could not be installed. Check library paths and permissions."))
      }
    }
  }
  missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  
  if(length(missing_packages)) {
    install.packages(missing_packages, dependencies = TRUE)
  }
  
  #functions
  Scatter_Density <- function(data = data, batch = batch, trt = NULL, expl.var = expl.var,
                              xlim = xlim, ylim = ylim, batch.legend.title = 'Batch', 
                              trt.legend.title = 'Treatment', density.lwd = 0.2,
                              title = NULL, title.cex = 1.5, legend.cex = 0.7, legend.title.cex =0.75){
    data = as.data.frame(data)
    batch = as.factor(batch)
    trt = as.factor(trt)
    if(nlevels(trt) >= 2){
      pMain <- ggplot(data = data, aes(x = data[ ,1], y = data[ ,2], colour = batch, shape = trt)) + 
        geom_point() + xlab(paste0('PC1: ', round(as.numeric(expl.var[1])*100), '% expl.var')) + 
        ylab(paste0('PC2: ', round(as.numeric(expl.var[2])*100), '% expl.var')) + 
        scale_color_manual(values = color.mixo(1:20)) + theme_bw() + xlim(xlim[1], xlim[2]) + 
        ylim(ylim[1], ylim[2]) + labs(colour = batch.legend.title, shape = trt.legend.title)
      
      pTop <- ggplot(data,aes(x = data[ ,1], fill = batch, linetype = trt)) + 
        geom_density(size = density.lwd, alpha = 0.5) + ylab('Density') + 
        theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.8)), 
              plot.title = element_text(hjust = 0.5, size = rel(title.cex)), 
              axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
              panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:20)) +
        xlim(xlim[1], xlim[2]) + labs(title = title)
      
      pRight <- ggplot(data, aes(x=data[ ,2], fill = batch, linetype = trt)) + 
        geom_density(size = density.lwd,alpha = 0.5) +  coord_flip() + ylab('Density') +
        theme(axis.title.x = element_text(size = rel(0.8)), 
              axis.title.y = element_blank(), axis.line = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(),
              panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:20)) +
        xlim(ylim[1], ylim[2])
      
    }else{
      pMain <- ggplot(data = data, aes(x = data[ ,1], y=data[ ,2], colour = batch)) + 
        geom_point(shape = 16) + xlab(paste0('PC1: ', round(as.numeric(expl.var[1])*100), '% expl.var')) + 
        ylab(paste0('PC2: ', round(as.numeric(expl.var[2])*100), '% expl.var')) + 
        scale_color_manual(values = color.mixo(1:20)) + theme_bw() + xlim(xlim[1], xlim[2]) + 
        ylim(ylim[1], ylim[2]) + labs(colour = batch.legend.title)
      
      pTop <- ggplot(data, aes(x = data[ ,1], fill = batch)) + 
        geom_density(size = density.lwd, alpha=0.5) + ylab('Density') + 
        theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.8)), 
              plot.title = element_text(hjust = 0.5, size = rel(title.cex)), 
              axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
              panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:20)) +
        xlim(xlim[1], xlim[2]) + labs(title = title)
      
      pRight <- ggplot(data, aes(x=data[ ,2], fill = batch)) + 
        geom_density(size = density.lwd, alpha = 0.5) +  coord_flip() + ylab('Density') +
        theme(axis.title.x = element_text(size = rel(0.8)), 
              axis.title.y = element_blank(), axis.line = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(),
              panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:20)) +
        xlim(ylim[1], ylim[2])
    }
    
    
    g <- ggplotGrob(pMain + theme(legend.position = 'right', legend.box = 'horizontal',
                                  legend.direction = 'vertical', 
                                  legend.key.height = unit(0.2, 'cm'),
                                  legend.key.width = unit(0.1, 'cm'),
                                  legend.title = element_text(size = rel(legend.title.cex)),
                                  legend.spacing.x = unit(0.1, 'cm'),
                                  legend.spacing.y = unit(0.1, 'cm'),
                                  legend.text = element_text(size = rel(legend.cex))))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    
    grid.arrange(pTop + theme(legend.position = 'none'), legend, pMain + 
                   theme(legend.position = 'none'), pRight + theme(legend.position = 'none'), 
                 ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
    
  }
  box_plot_fun = function(data = data, x = x, y = y, title = NULL, batch.legend.title = 'Batch', 
                          x.angle = 0, x.hjust = 0.5){
    ggplot(data = data, aes(x = x, y = y, fill = x)) + stat_boxplot(geom = "errorbar", width = 0.4) + 
      geom_boxplot() + scale_fill_manual(values = color.mixo(1:20)) + theme_bw() + 
      theme(axis.text.x = element_text(angle = x.angle, hjust = x.hjust), panel.grid = element_blank(),
            axis.title.x = element_blank(), axis.text = element_text(size = 10),
            axis.title = element_text(size = 12),
            plot.title = element_text(hjust = 0.5,size = rel(1))) + 
      labs(fill = batch.legend.title, y = 'value',title = title) 
  } 
  RleMicroRna2 <- function (object, maintitle = NULL, batch = batch, xlab = NA,
                            legend = TRUE, cex.lab = 1.2, cex.xaxis = 1, 
                            cex.yaxis = 1, abline.lwd=0.5, legend.cex = 0.8,
                            xaxis.dist.ratio = 0.1, outcex = 1, title.cex = 1.3) 
  {
    colorfill = color.mixo(batch)
    nARR = dim(object)[2]
    nGEN = dim(object)[1]
    y = apply(object, 1, median)
    mva = matrix(nrow = nGEN, ncol = nARR)
    for (i in 1:nARR) {
      x = object[, i]
      mva[ ,i] = (x - y)
    }
    med = apply(mva, 2, median)
    MIN = min(mva, na.rm = TRUE)
    MAX = max(mva, na.rm = TRUE)
    par(las = 3)
    plot(med, xlim = c(0, nARR + 1), ylim = c(MIN, MAX), axes = FALSE, 
         xlab = xlab, ylab = 'Deviations',cex.lab = cex.lab)
    colnames(mva) = colnames(object)
    res = boxplot(data.frame(mva), outline = TRUE, add = TRUE, col = colorfill,
                  xaxt = 'n', outcex = outcex, cex.axis = cex.yaxis) #outcex for outlier
    axis(1, cex.axis = cex.xaxis, at = 1:ncol(object), labels = NA)
    points(med, type = 'p', col = 'blue', cex = outcex)
    lines(med, type = 'l', col = 'blue', lty = 'dotted')
    title(main = maintitle, cex.main = title.cex)
    abline(0, 0, col = 'red', lwd = abline.lwd)
    par(las = 0)
    end_point = 0.5 + ncol(object)  # add degrees to the x axis
    box.max = max(max(res$stats), max(res$out))
    box.min = min(min(res$stats), min(res$out))
    box.range = box.max - box.min
    text(seq(1.2, end_point, by = 1), par("usr")[3] - xaxis.dist.ratio*box.range, 
         srt = 60, adj = 1, xpd = TRUE,
         labels = paste(colnames(object)), cex = cex.xaxis)
    if(legend == TRUE){
      legend('topright', legend = unique(batch), pch=15, col = unique(colorfill), cex = legend.cex)
    }
  }
  percentileofscore = function(df, control.index){
    df.percentile = df
    df.percentile[1:nrow(df), 1:ncol(df)] = NA
    for(i in 1:ncol(df)){
      control = sort(df[control.index, i])
      for(j in 1:nrow(df)){
        percentile.strick = sum(control < df[j, i])/length(control)
        percentile.weak = (length(control) - sum(control > df[j, i]))/length(control)
        percentile = (percentile.strick + percentile.weak)/2
        df.percentile[j, i] = percentile
        
      }
    }
    return(df.percentile)
  }
  ################
  percentile_norm = function(data = data, batch = batch, trt = trt){
    batch = as.factor(batch)
    trt = as.factor(trt)
    trt.list = list()
    data.pn.df = data.frame()
    for(i in 1:nlevels(batch)){
      trt.each.b = trt[batch == levels(batch)[i]]
      trt.list[[i]] = trt.each.b
      data.each.b.pn = percentileofscore(data[batch == levels(batch)[i],], 
                                         which(trt.each.b == levels(trt.each.b)[1]))
      data.pn.df = rbind(data.pn.df,data.each.b.pn)
    }
    names(trt.list) = levels(batch)
    data.pn.df.reorder = data.pn.df[rownames(data), ]
    return(data.pn.df.reorder)
  }
  #Silhouette coefficient
  
  calc.sil = function(
    x, # the PC variates
    y1, y2 = NULL, # factor of interest, e.g. known batch info or known treatment info
    name.y1, name.y2 = NULL # character of the factor of interest
  ){
    
    # calculate the distance, here euclidean is appropriate for PCA, NOT for t-SNE
    dist.res = daisy(x, metric = 'euclidean')
    # for factor 1
    sil.batch.res1 = silhouette(x = as.numeric(y1), dist = dist.res)
    # if factor 2 is provided
    if(!is.null(y2)) sil.batch.res2 = silhouette(x = as.numeric(y2), dist = dist.res)
    
    # extract average width silhouette per level
    res1 = c(summary(sil.batch.res1)['clus.avg.widths']$clus.avg.widths)
    names(res1) = levels(y1)
    
    
    if(!is.null(y2)){
      res2 = c(summary(sil.batch.res2)['clus.avg.widths']$clus.avg.widths)
      names(res2) = levels(y2)
    }
    
    # output data for plotting
    if(!is.null(y2)){
      silh.coeff = c(res1, res2)
      Cluster = c(levels(y1), levels (y2))
      Type = c(rep(name.y1, nlevels(y1)), rep(name.y2, nlevels(y2)))
      data.plot = data.frame(silh.coeff, Cluster, Type)
      
    }else{
      silh.coeff = c(res1)
      Cluster = c(levels(y1))
      Type = rep(name.y1, nlevels(y1))
      data.plot = data.frame(silh.coeff, Cluster, Type)
    }
    
    return(invisible(data.plot))
  }
  
  analyze_df <- function(df, plot_path, result_path, meta, batch_param) {
    
    color.mixo <- c('red', 'green', 'blue', 'yellow')
    #color.mixo <- color.mixo(1:20)
    
    # Check if the result path exists, create it if necessary
    if (!dir.exists(result_path)) {
      dir.create(result_path, recursive = TRUE)
    }
    
    for (col in 1:ncol(df)) {
      col_name <- colnames(df)[col]  # Get column name
      
      before.df <- data.frame(value = df[, col], batch = batch_param)
      
      # Box plot
      box_plot_fun(data = before.df, x = before.df$batch,
                   y = before.df$value, title = paste('B_P', col_name),
                   batch.legend.title = 'batch')
      ggsave(paste(plot_path, '/box_plot_', col, '.png', sep = ''), width = 6, height = 4)
      
      # Density plot
      density_plot <- ggplot(before.df, aes(x = value, fill = batch)) +
        geom_density(alpha = 0.5) +
        scale_fill_manual(values = color.mixo[1:20]) +
        labs(title = paste('D_P', col_name), x = 'Value', fill = 'experiment') +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid = element_blank())
      ggsave(paste(plot_path, '/density_plot_', col, '.png', sep = ''), width = 6, height = 4)
      
      # Linear regression summary
      
      
      data.lm <- lm(df[, col] ~ batch_param)
      summary_file <- paste(result_path, '/lm_summary_', col, '.txt', sep = '')
      sink(summary_file)
      summary(data.lm)
      sink(NULL)
    }
  }
  
  
  #log messages with different levels
  log_message <- function(message, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] [%s] %s\n", timestamp, level, message))
  }
  
  #main directory
  setwd("~")
  main_dir <- readline("Enter the main folder path (or press Enter to use the current directory): ")
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
  
  #prompt for results folder name 
  if (output_folder_prompt) {
    output_folder <- readline(prompt = "Enter the name for the results folder (default: 'batch_results_1'): ")
    if (output_folder == "") {
      output_folder <- "batch_results_1"
      log_message("No input provided. Using default results folder: 'batch_results_1'", "INFO")
    } else {
      log_message(sprintf("Results folder set to: '%s'", output_folder), "INFO")
    }
  } else {
    output_folder <- "batch_results_1"
    log_message("Using default results folder: 'batch_results_1'", "INFO")
  }
  
  #output directories with error handling
  create_directory <- function(path) {
    if (!dir.exists(path)) {
      tryCatch({
        dir.create(path, recursive = TRUE)
        log_message(sprintf("Created directory: '%s'", path), "INFO")
      }, error = function(e) {
        log_message(sprintf("Failed to create directory '%s': %s", path, e$message), "ERROR")
        stop(e)
      })
    } else {
      log_message(sprintf("Directory already exists: '%s'", path), "DEBUG")
    }
  }
  
  #subdirectories
  subdirs <- list(
    box_plots = "box_plots",
    lm_summary = "lm_summary",
    RLE_plots = "RLE_plots",
    heatmap = "heatmap",
    var_calc = "var_calc",
    batch_type = "batch_type"
  )
  
  #main output folder and subdirectories
  create_directory(output_folder)
  for (subdir in subdirs) {
    create_directory(file.path(output_folder, subdir))
  }
  
  #data and meta files read
  read_csv_safe <- function(filepath) {
    tryCatch({
      data <- read.csv(filepath, header = TRUE, stringsAsFactors = FALSE)
      log_message(sprintf("Successfully read file: '%s'", filepath), "DEBUG")
      return(data)
    }, error = function(e) {
      log_message(sprintf("Error reading file '%s': %s", filepath, e$message), "ERROR")
      stop(e)
    })
  }
  
  #initial folder inside dir
  data_folder <- readline("Enter the name of the folder containing your plate files: ")
  if (!dir.exists(data_folder)) {
    stop("Error: The specified folder does not exist within the main directory.")
  }
  
  data_file <- file.path(data_folder, 'feat_data.csv')
  meta_file <- file.path(data_folder, 'meta_data.csv')
  
  data <- read_csv_safe(data_file)
  meta <- read_csv_safe(meta_file)
  cli_alert_success("main folder, output folder and csv files were handled successfully!")
  cat("\033[1;31m==== Step 1 - processing data for the analysis... ====\033[0m\n")
  
  #prompt for row names
  if (rowname_prompt) {
    tryCatch({
      cat("Available columns for row names:\n")
      print(colnames(data))
      rowname_column <- readline(prompt = "Enter the column name for row names: ")
      
      if (!(rowname_column %in% colnames(data))) {
        stop(sprintf("Column '%s' does not exist in data.", rowname_column))
      }
      
      original_row_names <- data[[rowname_column]]
      rownames(data) <- data[[rowname_column]]
      names(data)[1] <- "rowname_column"
      data <- subset(data, select = -c(rowname_column))
      
      log_message(sprintf("Assigned row names using column: '%s'", rowname_column), "INFO")
    }, error = function(e) {
      log_message(sprintf("Error in row name assignment: %s", e$message), "ERROR")
      stop(e)
    })
  }
  
  #user to choose treatment and batch parameters
  get_meta_param <- function(meta, prompt_text) {
    tryCatch({
      cat(sprintf("Available columns in meta data for %s:\n", prompt_text))
      print(colnames(meta))
      param <- readline(prompt = sprintf("Enter the column name for %s: ", prompt_text))
      
      if (!(param %in% colnames(meta))) {
        stop(sprintf("Column '%s' does not exist in meta data.", param))
      }
      
      return(param)
    }, error = function(e) {
      log_message(sprintf("Error in selecting %s parameter: %s", prompt_text, e$message), "ERROR")
      stop(e)
    })
  }
  
  trt_param <- ""
  batch_param <- ""
  
  if (batch_param_prompt) {
    trt_param <- get_meta_param(meta, "treatment parameter( factor of interest )")
    batch_param <- get_meta_param(meta, "batch parameter")
  }
  
  #selected parameters
  log_message(sprintf("Treatment parameter selected: '%s'", trt_param), "INFO")
  log_message(sprintf("Batch parameter selected: '%s'", batch_param), "INFO")
  cli_alert_success("Data processing complete. Proceeding to analysis.")
  
  cat("\033[1;31m==== Step 2 - Batch effect analysis began... ====\033[0m\n")
  #log-ratio transformation
  cat("Normalizing data and saving normalized data frame to output directory.\n")
  tryCatch({
    data.clr <- logratio.transfo(data, logratio = 'CLR')
    class(data.clr) <- 'matrix'
    log_message("Performed CLR transformation on data.", "DEBUG")
  }, error = function(e) {
    log_message(sprintf("Error in CLR transformation: %s", e$message), "ERROR")
    stop(e)
  })
  
  # Save processed data.clr as CSV
  tryCatch({
    data_clr_df <- as.data.frame(data.clr)
    data_clr_df$original_row_names <- original_row_names
    clr_output_file <- file.path(output_folder, "data.clr_df.csv")
    write.csv(as.data.frame(data_clr_df), file = clr_output_file, row.names = FALSE)
    log_message(sprintf("Processed data saved to: '%s'", clr_output_file), "INFO")
  }, error = function(e) {
    log_message(sprintf("Error saving CLR data: %s", e$message), "ERROR")
    stop(e)
  })
  
  #PCA for batch effect detection
  cat("Generating PCA...\n")
  tryCatch({
    data.pca <- pca(data.clr, ncomp = 3)
    # Plot PCA before batch effect correction
    data.pca.plot <- Scatter_Density(data = data.pca$variates$X,
                                     batch = meta[[batch_param]],
                                     trt = meta[[trt_param]],
                                     expl.var = data.pca$explained_variance,
                                     xlim = c(-4.5,5), ylim = c(-3,4),
                                     batch.legend.title = 'Chosen Batch Parameter',
                                     trt.legend.title = 'Chosen Treatment Parameter',
                                     title = 'PCA Before Batch Effect Correction')
    
    if (!is.null(data.pca.plot)) {
      #print PCA plot
      print(data.pca.plot)
      pca_plot_file <- file.path(output_folder, "pca_plot.png")
      #save PCA plot
      ggsave(pca_plot_file, plot = data.pca.plot, width = 12, height = 8, dpi = 300)
      cat("PCA plot saved to:", pca_plot_file, "\n")
    } else {
      cat("PCA plot is empty. Skipping saving.\n")
    }
  })
  
  #analyze_df function with error handling
  cat("Generating BOX & Density plots for each feature in the analysis...\n")
  tryCatch({
    plot_subdir = file.path(output_folder, "box_plots")
    lm_subdir = file.path(output_folder, "lm_summary")
    analyze_df(data.clr, plot_subdir, lm_subdir, meta, meta[[batch_param]])
    log_message("Box plots generated successfully.", "INFO")
  }, error = function(e) {
    log_message(sprintf("Error in analyze_df function: %s", e$message), "ERROR")
  })
  
  #RLE plots
  cat("Generating RLE plots for treatment groups...\n")
  tryCatch({
    meta.trt <- meta[[trt_param]]
    meta.trt <- as.factor(meta.trt)
    meta.batch <- meta[[batch_param]]
    meta.batch <- as.factor(meta.batch)
    
    unique_treatment_factors <- unique(meta.trt)
    log_message(sprintf("Unique treatment factors: %s", paste(unique_treatment_factors, collapse = ", ")), "DEBUG")
    
    for (factor in unique_treatment_factors) {
      #subset data for the current treatment factor
      current_data <- data.clr[meta.trt == factor, ]
      batch_data <- meta.batch[meta.trt == factor]
      current_data_copy <- data_clr_df[meta.trt == factor, ]
      
      #saving subsetted data frame
      subset_data_file <- file.path(output_folder, subdirs$RLE_plots, paste0("subset_data_", factor, ".csv"))
      write.csv(current_data_copy, file = subset_data_file, row.names = TRUE)
      log_message(sprintf("Subset data for '%s' saved to: '%s'", factor, subset_data_file), "DEBUG")
      
      batch_data_file <- file.path(output_folder, subdirs$RLE_plots, paste0("batch_data_", factor, ".csv"))
      write.csv(batch_data, file = batch_data_file, row.names = TRUE)
      log_message(sprintf("Batch data for '%s' saved to: '%s'", factor, batch_data_file), "DEBUG")
      
      #RLE on subset of data
      rle_plot_file <- file.path(output_folder, subdirs$RLE_plots, paste0("RLE_plot_", factor, ".png"))
      png(rle_plot_file)
      tryCatch({
        #RleMicroRna2 function
        RleMicroRna2(object = t(current_data), batch = batch_data,
                     maintitle = paste('Batch Analysis (Type =', factor, ')'))
        dev.off()
        log_message(sprintf("RLE plot for '%s' saved to: '%s'", factor, rle_plot_file), "INFO")
      }, error = function(e) {
        dev.off()
        log_message(sprintf("Error in generating RLE plot for '%s': %s", factor, e$message), "ERROR")
      })
    }
    log_message("RLE plots generated successfully.", "INFO")
  }, error = function(e) {
    log_message(sprintf("Error in RLE plotting section: %s", e$message), "ERROR")
  })
  cat("Completed RLE analysis!\n")
  
  #Linear regression
  cat("Working on Linear Regression analysis...\n")
  tryCatch({
    data.anova_t <- apply(data.clr, 2, function(x){
      res.lm <- lm(x ~ meta.trt + meta.batch)
      a <- anova(res.lm)
      return(a)
    })
    
    anova_data_file <- file.path(output_folder, subdirs$lm_summary, "anova_data.csv")
    write.csv(data.anova_t, file = anova_data_file, row.names = TRUE)
    log_message(sprintf("ANOVA data saved to: '%s'", anova_data_file), "INFO")
    
    data.rsq <- apply(data.clr, 2, function(x){
      res.lm <- lm(x ~ meta.trt + meta.batch)
      summary.res <- summary(res.lm)
      r <- summary.res$adj.r.squared
      return(r)
    })
    
    rsq_data_file <- file.path(output_folder, subdirs$lm_summary, "rsq_data.csv")
    write.csv(data.rsq, file = rsq_data_file, row.names = TRUE)
    log_message(sprintf("R-squared data saved to: '%s'", rsq_data_file), "INFO")
  }, error = function(e) {
    log_message(sprintf("Error in linear regression analysis: %s", e$message), "ERROR")
  })
  cat("Completed Linear Regression analysis!\n")
  # Heatmap
  cat("Working on Heatmap for batch detection...\n")
  tryCatch({
    #scaling data
    data.clr.scale <- scale(data.clr, center = TRUE, scale = TRUE)
    data.clr.scale <- scale(t(data.clr.scale), center = TRUE, scale = TRUE)
    
    #annotation for columns
    data.anno_col <- data.frame(Batch = meta.batch, Tissue = meta.trt)
    rownames(data.anno_col) <- colnames(data.clr.scale)
    
    #colors using RColorBrewer
    batch_levels <- levels(meta.batch)
    trt_levels <- levels(meta.trt)
    batch_colors <- brewer.pal(min(length(batch_levels), 8), "Set1")
    trt_colors <- brewer.pal(min(length(trt_levels), 8), "Set2")
    data.anno_metabo_colors <- list(
      Batch = setNames(batch_colors, batch_levels),
      Tissue = setNames(trt_colors, trt_levels)
    )
    
    #Heatmap
    heatmap_file <- file.path(output_folder, subdirs$heatmap, "heatmap.png")
    pheatmap(data.clr.scale,
             scale = 'none',
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             fontsize_row = 2,
             fontsize_col = 2,
             fontsize = 3,
             clustering_distance_rows = 'euclidean',
             clustering_method = 'ward.D',
             treeheight_row = 30,
             annotation_col = data.anno_col,
             annotation_colors = data.anno_metabo_colors,
             border_color = NA,
             main = 'Heatmap - Scaled Data',
             filename = heatmap_file,
             width = 12, height = 8)
    log_message(sprintf("Heatmap saved to: '%s'", heatmap_file), "INFO")
  }, error = function(e) {
    log_message(sprintf("Error in generating heatmap: %s", e$message), "ERROR")
  })
  cat("Completed Heatmap on all features!\n")
  #variance calculation and plotting
  tryCatch({
    #variation explained for each feature
    variation_explained <- apply(t(data.clr), 1, function(x) var(x) / sum(var(x)))
    
    #sorting features based on variation explained
    sorted_features <- sort(variation_explained, decreasing = TRUE)
    
    #top features
    top_features <- names(sorted_features)[1:min(num_features, length(sorted_features))]
    
    #subset data for top features
    data_sub <- t(data.clr)[top_features, ]
    data.clr.scale <- scale(data.clr.scale, center = TRUE, scale = TRUE)
    data.clr.scale <- scale(t(data_sub), center = TRUE, scale = TRUE)
    data.clr.scale <- t(data.clr.scale)
    
    #annotation remains the same
    data.anno_col <- data.frame(Batch = meta.batch, Tissue = meta.trt)
    rownames(data.anno_col) <- colnames(data.clr.scale)
    batch_levels <- levels(meta.batch)
    trt_levels <- levels(meta.trt)
    batch_colors <- brewer.pal(min(length(batch_levels), 8), "Set1")
    trt_colors <- brewer.pal(min(length(trt_levels), 8), "Set2")
    data.anno_metabo_colors <- list(
      Batch = setNames(batch_colors, batch_levels),
      Tissue = setNames(trt_colors, trt_levels))
    
    #second Heatmap for top features
    heatmap_file2 <- file.path(output_folder, subdirs$heatmap, "heatmap_top_features.png")
    pheatmap(data.clr.scale,
             scale = 'none',
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             fontsize_row = 2,
             fontsize_col = 2,
             fontsize = 3,
             clustering_distance_rows = 'euclidean',
             clustering_method = 'ward.D',
             treeheight_row = 30,
             annotation_col = data.anno_col,
             annotation_colors = data.anno_metabo_colors,
             border_color = NA,
             main = 'Heatmap - Top Features Scaled',
             filename = heatmap_file2,
             width = 12, height = 8)
    log_message(sprintf("Top features heatmap saved to: '%s'", heatmap_file2), "INFO")
  }, error = function(e) {
    log_message(sprintf("Error in generating top features heatmap: %s", e$message), "ERROR")
  })
  cat("Completed Heatmap on selected features!\n")
  #variance analysis
  cat("Working on variance analysis...\n")
  tryCatch({
    #variance calculation
    data.form <- ~ meta.trt + meta.batch
    data.info <- as.data.frame(cbind(rownames(data.clr), meta.trt, meta.batch))
    rownames(data.info) <- rownames(data.clr)
    
    #variance partition model
    data.varPart.before <- fitExtractVarPartModel(exprObj = t(data.clr),
                                                  formula = data.form,
                                                  data = data.info)
    data.varmat.before <- as.matrix(data.varPart.before[, 1:2])
    
    #preparing data for plotting
    data.variance <- c(as.vector(data.varmat.before))
    data.variance <- cbind(variance = data.variance,
                           Type = rep(c('sample', 'Gender'), each = ncol(data.clr)),
                           method = rep(c('Before'), each = 2 * ncol(data.clr)))
    data.variance <- as.data.frame(data.variance)
    data.variance$method <- factor(data.variance$method,
                                   levels = unique(data.variance$method))
    data.variance$variance <- as.numeric(as.character(data.variance$variance))
    
    #plot variance
    var_plot <- ggplot(data.variance, aes(x = Type, y = variance, fill = Type)) +
      geom_boxplot() +
      facet_grid(cols = vars(method)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1),
            strip.text = element_text(size = 12),
            panel.grid = element_blank(),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12)) +
      labs(x = 'Type', y = 'Proportion Variance', fill = 'Type') +
      ylim(0, 1)
    
    var_plot_file <- file.path(output_folder, subdirs$var_calc, "varmat.pdf")
    ggsave(filename = var_plot_file, plot = var_plot, width = 10, height = 6, dpi = 300)
    log_message(sprintf("Variance plot saved to: '%s'", var_plot_file), "INFO")
  }, error = function(e) {
    log_message(sprintf("Error in variance plotting: %s", e$message), "ERROR")
  })
  cat("Completed Variance anaysis 1!\n")
  #variance calculation using RDA
  tryCatch({
    
    data.design <- numeric()
    data.design$group <- meta.trt
    data.design$batch <- meta.batch
    
    
    data.rda.before1 <- rda(data.clr ~ group + Condition(batch), 
                            data = data.design)
    data.rda.before2 <- rda(data.clr ~ batch + Condition(group), 
                            data = data.design)
    
    
    # amount of variance
    data.rda.bat_prop.before <- data.rda.before1$pCCA$tot.chi*100/data.rda.before1$tot.chi
    data.rda.trt_prop.before <- data.rda.before2$pCCA$tot.chi*100/data.rda.before2$tot.chi
    data.rda.prop.before <- c(data.rda.bat_prop.before, 
                              data.rda.trt_prop.before)
    
    # merge results
    data.rda.prop.val <- c(data.rda.prop.before)
    
    # add batch, trt and method info
    data.rda.prop <- data.frame(prop = data.rda.prop.val, 
                                prop.r = round(data.rda.prop.val, 2), 
                                Method = rep(c('Before'), each = 2), 
                                Type = rep(c('experiment', 'sample'), 6))
    
    # reorder levels
    data.rda.prop$Method <- factor(data.rda.prop$Method, 
                                   levels = unique(data.rda.prop$Method))
    
    
    #plot variance explained by RDA
    var_plot2 <- ggplot(data = data.rda.prop, aes(x = Method, y = prop, fill = Type)) +
      geom_bar(stat = "identity", position = 'dodge', colour = 'black') +
      geom_text(aes(label = prop.r),
                position = position_dodge(width = 0.9),
                vjust = -0.5, size = 3) +
      theme_bw() +
      labs(y = "Variance Explained (%)") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1),
            panel.grid = element_blank(),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12)) +
      ylim(0, 100)
    
    var_plot_file2 <- file.path(output_folder, subdirs$var_calc, "var_calc_rda.pdf")
    ggsave(filename = var_plot_file2, plot = var_plot2, width = 10, height = 6, dpi = 300)
    log_message(sprintf("Variance plot 2 saved to: '%s'", var_plot_file2), "INFO")
  }, error = function(e) {
    log_message(sprintf("Error in RDA variance plotting: %s", e$message), "ERROR")
  })
  cat("Completed variance calculation analysis!\n")
  
  #batch type
  cat("Working on batch type visualzation...\n")
  tryCatch({
    data.b.coeff <- c()
    for(i in 1:ncol(data.clr)){
      res <- lm(data.clr[ ,i] ~ meta.trt + meta.batch)
      sum.res <- summary(res)
      data.b.coeff[i] <- sum.res$coefficients[3,1]
    }
    batch_type_plot_file <- file.path(output_folder, "batch_type", "batch_type_plot.png")
    png(batch_type_plot_file, width = 1200, height = 800)
    par(mfrow = c(2,2))
    hist(data.b.coeff,col = 'gray')
    plot(density(data.b.coeff))
    qqnorm(data.b.coeff)
    qqline(data.b.coeff, col='red')
    par(mfrow = c(1,1))
    dev.off()
    
    log_message(sprintf("Batch type plot saved to: '%s'", batch_type_plot_file), "INFO")
  }, error = function(e) {
    log_message(sprintf("Error in batch type plotting: %s", e$message), "ERROR")
  })
  cat("Completed batch type visualzation!\n")
  cli_alert_success(" Completed batch effect analysis successfully!")
  log_message("Batch effect detection completed successfully.", "INFO")
  cat("\033[1;31m==== Analysis completed, review files and run again as needed:)  ====\033[0m\n")
  cat("\033[1;31m==== IMPORTANT! PLEASE REMOVE YOUR FILES FROM THE DIRECTORY ====\033[0m\n")
}

detect_batch_effects(output_folder_prompt = TRUE, rowname_prompt = TRUE, batch_param_prompt = TRUE, num_features = 20)
#############################################################################################################

