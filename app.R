# Load necessary libraries
library(tidyverse)
library(limma)
library(randomForest)
library(varImp)
library(caret)
library(FactoMineR)
library(ggfortify)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(shiny)
library(shinythemes)
library(DT)
library(umap)
library(doParallel)
library(foreach)
library(enrichR)
library(plotly)
library(shinyjs)  # Add this library at the beginning
library(randomForest)
library(pwr)
library(nnet)
library(pROC)

# Define EnrichR databases to use
enrichr_dbs <- c("KEGG_2021_Human", "GO Molecular Function 2023", "GO Biological Process 2023", "GWAS Catalog 2023", "UK Biobank GWAS v1")

# Remove Outlier Samples Function
remove_outliers <- function(df) {
  df$phenotype <- as.numeric(df$phenotype)  # Convert phenotype to numeric
  sample_sums <- rowSums(df[, -ncol(df)])  # Exclude the phenotype column
  
  Q1 <- quantile(sample_sums, 0.25)
  Q3 <- quantile(sample_sums, 0.75)
  IQR <- Q3 - Q1
  
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  non_outliers <- rownames(df)[sample_sums >= lower_bound & sample_sums <= upper_bound]
  
  df_clean <- df[rownames(df) %in% non_outliers, ]
  df_clean <- as.data.frame(df_clean)  # Ensure df_clean is a data frame
  
  df_clean$phenotype <- as.factor(df_clean$phenotype)  # Convert phenotype back to factor
  
  return(df_clean)
}
dif_expr <- function(expr_data, phenotype_data, chip) {
  set.seed(1)
  
  # Load expression data
  if (grepl("\\.csv$", expr_data)) {
    Exprmat <- read.csv(expr_data, row.names = 1)
  } else {
    Exprmat <- read.table(expr_data, header = TRUE, row.names = 1)
  }
  
  # Load phenotype data
  if (grepl("\\.csv$", phenotype_data)) {
    disease <- read.csv(phenotype_data, row.names = 1)
  } else {
    disease <- read.table(phenotype_data, header = TRUE, row.names = 1)
  }
  
  # Convert phenotype to factor
  disease_fact <- as.factor(disease[, 1])
  
  # Transpose expression matrix and add phenotype column
  df <- data.frame(t(Exprmat))
  df$phenotype <- disease_fact
  
  # Check for mismatched sample sizes
  if (nrow(df) != length(disease_fact)) {
    stop("The number of samples in expr_data and phenotype_data do not match.")
  }
  
  # Remove outliers
  df <- remove_outliers(df)
  
  # Separate phenotype column
  phenotype <- df$phenotype
  df <- df[, colnames(df) != "phenotype"]
  
  # Debugging prints
  print("Expression matrix before scaling:")
  print(head(df))
  
  # Scale the expression data
  df_scaled <- as.data.frame(scale(df))
  
  # Debugging prints
  print("Expression matrix after scaling:")
  print(head(df_scaled))
  
  # Reattach the phenotype column
  df_scaled$phenotype <- phenotype
  
  # Identify good genes
  good_stuff <- Exprmat[rowSums(Exprmat) > 10, ]
  good_stuff$goodgenes <- row.names(good_stuff)
  
  # Debugging prints
  print("Good genes identified:")
  print(good_stuff$goodgenes)
  
  print("Number of good genes identified:")
  print(length(good_stuff$goodgenes))
  
  # Remove 'X' prefix from column names in df_scaled
  colnames(df_scaled) <- gsub("^X", "", colnames(df_scaled))
  
  print("Column names in df_scaled before subsetting:")
  print(colnames(df_scaled))
  
  # Check for intersection of good genes with df_scaled columns
  intersect_genes <- intersect(good_stuff$goodgenes, colnames(df_scaled))
  
  print("Intersecting good genes with df_scaled columns:")
  print(intersect_genes)
  
  # Subset scaled data to include only good genes and phenotype
  cols_to_include <- c(intersect_genes, "phenotype")
  df_scaled <- df_scaled[, colnames(df_scaled) %in% cols_to_include]
  
  # Debugging prints
  print("Column names after subsetting:")
  print(colnames(df_scaled))
  
  # Check for missing phenotype column
  if (!"phenotype" %in% colnames(df_scaled)) {
    stop("Phenotype column is missing after subsetting.")
  }
  
  # Transpose scaled data (excluding phenotype)
  exprs <- t(df_scaled[, colnames(df_scaled) != "phenotype"])
  
  # Create design matrix
  design <- model.matrix(~0 + df_scaled$phenotype)
  colnames(design) <- levels(df_scaled$phenotype)
  
  # Perform differential expression analysis
  fit <- lmFit(exprs, design)
  fit <- eBayes(fit)
  
  # Extract top table results
  initial_results <- topTable(fit, coef =1, number = Inf, adjust.method = "BH",
                              sort.by = "B", p.value = 0.05, lfc = 0)
  
  initial_results$Probe_IDs <- row.names(initial_results)
  
  # Select relevant columns and order by adjusted p-value
  initial_results <- initial_results[, c("Probe_IDs", "adj.P.Val", "logFC")]
  initial_results <- initial_results[order(initial_results$adj.P.Val), ]
  top50_results <- head(initial_results, 50)
  
  # Subset scaled data to include only top 50 results and phenotype
  df_scaled <- df_scaled[, colnames(df_scaled) %in% c(top50_results$Probe_IDs, "phenotype")]
  
  # Debugging prints
  print("Column names after filtering for top 50 results:")
  print(colnames(df_scaled))
  
  # Train random forest model using data frame interface
  model <- randomForest(x = df_scaled[, -which(colnames(df_scaled) == "phenotype")], 
                        y = df_scaled$phenotype, importance = TRUE)
  
  # Extract important features
  importance_df <- data.frame(importance(model))
  importance_df$probenames <- row.names(importance_df)
  importance_df <- importance_df[order(importance_df$MeanDecreaseAccuracy, decreasing = TRUE), ]
  Top20_RFgenes <- importance_df$probenames[1:20]
  top20_results <- top50_results[top50_results$Probe_IDs %in% Top20_RFgenes, ]
  print(top20_results)
  # Function to retry getBM with multiple attempts
  retry_getBM <- function(genes, mart, chip, retries = 5) {
    for (i in 1:retries) {
      tryCatch({
        return(getBM(filters = chip, 
                     attributes = c(chip, "ensembl_gene_id", "hgnc_symbol"), 
                     values = genes, mart = mart, bmHeader = TRUE, useCache = FALSE))
      }, error = function(e) {
        if (i == retries) stop(e)
        Sys.sleep(10)
      })
    }
  }
  
  # Use Ensembl biomart to retrieve gene information
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- top20_results$Probe_IDs
  
  # Debugging: Print chip being used
  print("Using chip:")
  print(chip)
  
  G_list <- retry_getBM(genes, mart, chip)
  # Modify chip name to match the format used in Ensembl data
  chip_mod <- gsub("_", " ", toupper(chip))
  print(chip_mod)
  colnames(G_list)[1] <- chip_mod
  
  # Merge results with gene information
  final_results <- merge(top20_results, G_list, by.x = "Probe_IDs", by.y =  chip_mod,all.x=TRUE)
  #renaming columns
  colnames(final_results)[1]<-"Probe_ID"
  colnames(final_results)[4] <- "Ensembl_ID"
  
  return(final_results)
}
#

# Power Calculation Function for Two or More Classes
calculate_power <- function(phenotype_data, effect_size = 0.5, alpha = 0.05) {
  # Count the number of samples in each class
  class_counts <- table(phenotype_data)
  
  # Get the number of groups
  k <- length(class_counts)
  
  if (k == 2) {
    # For two classes, perform a two-sample t-test power calculation
    n_group1 <- class_counts[1]
    n_group2 <- class_counts[2]
    
    # Perform power calculation for a two-sample t-test
    result <- pwr.t.test(n = min(n_group1, n_group2), d = effect_size, sig.level = alpha, power = NULL, type = "two.sample")
    
  } else {
    # For more than two classes, perform a one-way ANOVA power calculation
    n_harmonic <- 1 / mean(1 / class_counts)
    
    # Perform power calculation for a one-way ANOVA
    result <- pwr.anova.test(k = k, n = n_harmonic, f = effect_size, sig.level = alpha, power = NULL)
  }
  
  return(result)
}

# Regularized Random Forest Classification Function
regularized_random_forest <- function(expr_data, phenotype_data, selected_genes) {
  set.seed(1)
  
  # Load expression data
  if (grepl("\\.csv$", expr_data)) {
    Exprmat <- read.csv(expr_data, row.names = 1)
  } else {
    Exprmat <- read.table(expr_data, header = TRUE, row.names = 1)
  }
  
  # Load phenotype data
  if (grepl("\\.csv$", phenotype_data)) {
    disease <- read.csv(phenotype_data, row.names = 1)
  } else {
    disease <- read.table(phenotype_data, header = TRUE, row.names = 1)
  }
  
  # Prepare data frame
  df <- data.frame(t(Exprmat))
  colnames(df) <- gsub("^X", "", colnames(df))  # Remove X prefix
  df$phenotype <- as.factor(disease[, 1])
  
  # Select only the genes present in the `selected_genes`
  df_subset <- df[, selected_genes, drop = FALSE]
  df_subset$phenotype <- df$phenotype
  
  # Split the data into training and testing sets
  set.seed(123)
  train_index <- createDataPartition(df_subset$phenotype, p = 0.7, list = FALSE)
  train_data <- df_subset[train_index, ]
  test_data <- df_subset[-train_index, ]
  
  # Train Random Forest with regularization parameters
  model <- randomForest(
    x = train_data[, -which(colnames(train_data) == "phenotype")], 
    y = train_data$phenotype,
    importance = TRUE,
    ntree = 500,          # Number of trees
    mtry = 3,             # Number of features considered at each split
    maxnodes = 30,        # Maximum number of nodes in each tree
    nodesize = 5          # Minimum size of terminal nodes
  )
  
  # Make predictions
  predictions <- predict(model, newdata = test_data[, -which(colnames(test_data) == "phenotype")])
  
  # Evaluate model performance
  confusion_matrix <- confusionMatrix(predictions, test_data$phenotype)
  
  accuracy <- confusion_matrix$overall['Accuracy']
  sensitivity <- confusion_matrix$byClass['Sensitivity']
  specificity <- confusion_matrix$byClass['Specificity']
  auc <- roc(response = test_data$phenotype, predictor = as.numeric(predictions))$auc
  prevalence <- mean(test_data$phenotype == levels(test_data$phenotype)[1])
  
  metrics <- data.frame(
    Metric = c("Accuracy", "Sensitivity", "Specificity", "AUC", "Prevalence"),
    Value = c(accuracy, sensitivity, specificity, auc, prevalence)
  )
  
  return(list(model = model, metrics = metrics))
}


# Pathway Analysis Functions with Combined Results
enrich_genes <- function(genes) {
  if (length(genes) == 0) {
    stop("No genes provided for enrichment analysis.")
  }
  enrichr(genes, enrichr_dbs)
}

enrich_all_genes <- function(final_results) {
  genes <- final_results[, 5]  
  print("Enriching all genes:")
  print(genes)
  enrich_res <- enrich_genes(genes)
  combined_res_all <<- combine_enrichr_results(enrich_res)
  return(combined_res_all)
}

enrich_up_genes <- function(final_results) {
  up_genes <- final_results[, 5][final_results$logFC > 0]  # Use Ensembl_IDs for enrichment
  enrich_res <- enrich_genes(up_genes)
  combined_res_up <<- combine_enrichr_results(enrich_res)
  return(combined_res_up)
}

enrich_down_genes <- function(final_results) {
  down_genes <- final_results[, 5][final_results$logFC < 0]  # Use Ensembl_IDs for enrichment
  enrich_res <- enrich_genes(down_genes)
  combined_res_down <<- combine_enrichr_results(enrich_res)
  return(combined_res_down)
}

combine_enrichr_results <- function(enrichr_results) {
  combined_results <- do.call(rbind, lapply(names(enrichr_results), function(db) {
    df <- enrichr_results[[db]]
    if (nrow(df) > 0) {
      df$Database <- db
      return(df)
    } else {
      return(NULL)
    }
  }))
  return(as.data.frame(combined_results))
}

Make_Volcano <- function(results, p_value_threshold, log_fc_threshold) {
  results$hgnc_symbol <- results[, 5]  # rename column for ease.
  
  # Convert the results data frame to a named list to avoid issues with jsonlite
  results_list <- lapply(seq_along(results), function(i) as.list(results[, i]))
  names(results_list) <- names(results)
  
  volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), text = results_list)) +
    geom_point(aes(color = adj.P.Val < p_value_threshold & abs(logFC) > log_fc_threshold)) +
    scale_color_manual(values = c("black", "red")) +
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
    theme_minimal()
  
  return(ggplotly(volcano_plot, tooltip = "text"))
}


Make_PCA <- function(expr_data, phenotype_data, n_components) {
  # Read data
  if (grepl("\\.csv$", expr_data)) {
    Exprmat <- read.csv(expr_data, row.names = 1)
  } else {
    Exprmat <- read.table(expr_data, header = TRUE, row.names = 1)
  }
  
  if (grepl("\\.csv$", phenotype_data)) {
    disease <- read.csv(phenotype_data, row.names = 1)
  } else {
    disease <- read.table(phenotype_data, header = TRUE, row.names = 1)
  }
  
  disease_fact <- as.factor(disease[, 1])
  
  # Transpose and prepare the data frame
  df <- data.frame(t(Exprmat))
  df$phenotype <- disease_fact
  df$Sample <- rownames(df)
  
  # Debugging prints
  print("Column names of the initial data frame:")
  print(colnames(df))
  
  # Keep only numeric columns for PCA
  numeric_cols <- sapply(df, is.numeric)
  df_numeric <- df[, numeric_cols]
  
  # Debugging prints
  print("Column names of the numeric data frame:")
  print(colnames(df_numeric))
  
  # Remove columns with zero variance
  zero_var_columns <- sapply(df_numeric, function(x) var(x, na.rm = TRUE) == 0)
  print("Columns with zero variance:")
  print(names(zero_var_columns)[zero_var_columns])
  
  df_numeric <- df_numeric[, !zero_var_columns]
  
  # Ensure 'phenotype' and 'Sample' columns are not removed
  df_numeric$phenotype <- disease_fact
  df_numeric$Sample <- rownames(df)
  
  # Debugging prints
  print("Column names after removing zero variance columns:")
  print(colnames(df_numeric))
  
  # Perform PCA
  pca <- prcomp(df_numeric[, !colnames(df_numeric) %in% c("phenotype", "Sample")], scale. = TRUE)
  pca_df <- data.frame(pca$x)
  pca_df$Sample <- rownames(pca_df)
  pca_df$Phenotype <- df_numeric$phenotype[match(rownames(pca_df), rownames(df_numeric))]
  
  # Debugging prints
  print("PCA Data Frame:")
  print(head(pca_df))
  
  # Plot PCA results with hover functionality
  pca_plot <- ggplot(pca_df, aes_string(x = "PC1", y = "PC2", color = "Phenotype", text = "Sample")) +
    geom_point() +
    labs(title = "PCA Plot", x = "PC1", y = "PC2") +
    scale_color_discrete(name = "Phenotype") +  # Change legend title
    theme_minimal()
  
  return(ggplotly(pca_plot, tooltip = "text"))
}



# UMAP Function with Parallel Processing, Sliders, Hover, and Scaling
Make_UMAP <- function(expr_data, phenotype_data, n_neighbors) {
  if (grepl("\\.csv$", expr_data)) {
    Exprmat <- read.csv(expr_data, row.names = 1)
  } else {
    Exprmat <- read.table(expr_data, header = TRUE, row.names = 1)
  }
  
  if (grepl("\\.csv$", phenotype_data)) {
    disease <- read.csv(phenotype_data, row.names = 1)
  } else {
    disease <- read.table(phenotype_data, header = TRUE, row.names = 1)
  }
  
  disease_fact <- as.factor(disease[, 1])
  
  df <- data.frame(t(Exprmat))
  df$phenotype <- disease_fact
  df$Sample <- rownames(df)
  
  if (nrow(df) != length(disease_fact)) {
    stop("The number of samples in expr_data and phenotype_data do not match.")
  }
  
  n_neighbors <- min(n_neighbors, nrow(df) - 1)
  
  # Scale the expression data
  df_scaled <- df %>% mutate(across(-c(phenotype, Sample), scale))
  
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  
  umap_res <- foreach(i = 1, .combine = rbind, .packages = 'umap') %dopar% {
    umap(df_scaled[, !(colnames(df_scaled) %in% c("phenotype", "Sample"))], n_neighbors = n_neighbors)
  }
  
  stopCluster(cl)
  
  umap_df <- data.frame(umap_res$layout)
  umap_df$Sample <- rownames(umap_df)
  umap_df$Phenotype <- df$phenotype
  
  if (nrow(umap_df) != nrow(df)) {
    stop("The number of rows in UMAP result and input data do not match.")
  }
  
  umap_plot <- ggplot(umap_df, aes(x = X1, y = X2, color = Phenotype, text = Sample)) +
    geom_point() +
    labs(title = "UMAP Plot", x = "UMAP1", y = "UMAP2") +
    theme_minimal()
  
  return(ggplotly(umap_plot, tooltip = "text"))
}
Make_Volcano <- function(results, p_value_threshold, log_fc_threshold) {
  results$hgnc_symbol <- results[, 5]  # rename column for ease.
  volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), text = hgnc_symbol)) +
    geom_point(aes(color = adj.P.Val < p_value_threshold & abs(logFC) > log_fc_threshold)) +
    scale_color_manual(values = c("black", "red")) +
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
    theme_minimal()
  
  return(ggplotly(volcano_plot, tooltip = "text"))
}



library(shinyjs)  # Add this library at the beginning# UI

      

# UI
ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  theme = shinytheme("cyborg"),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('notify', function(message) {
        alert(message);
      });
    ")),
    # Add CSS for scrolling
    tags$style(HTML("
      #scrollable-content {
        overflow-y: auto;
        max-height: 600px; /* Adjust this value based on your needs */
        padding: 10px;
      }
    "))
  ),
  
  titlePanel("JCAP Differential Expression and Feature Selection on Microarray Data"),
  
  # Wrap sidebarLayout in a div with the ID "scrollable-content" to make it scrollable
  div(
    id = "scrollable-content",
    sidebarLayout(
      sidebarPanel(
        fileInput("counts_input", "Input Expression Data", accept = c(".csv", ".txt")),
        fileInput("phenotype_input", "Input Phenotype Data", accept = c(".csv", ".txt")),
        uiOutput("chip_ui"),  # Dynamic UI for chip selection
        actionButton("running_Dif", "Perform Differential Expression Analysis"),
        actionButton("plot_pca", "Display PCA Plot"),
        actionButton("plot_umap", "Display UMAP Plot"),
        actionButton("pathway_all", "Enrich Pathways (All Genes)"),
        actionButton("pathway_up", "Enrich Pathways (Upregulated Genes)"),
        actionButton("pathway_down", "Enrich Pathways (Downregulated Genes)"),
        actionButton("plot_volcano", "Display Volcano Plot"),
        actionButton("run_rf", "Run Random Forest Analysis"),

        actionButton("calculate_power", "Calculate Statistical Power"),
        sliderInput("n_components_pca", "Number of PCA Components", min = 2, max = 10, value = 2),
        sliderInput("n_neighbors_umap", "Number of UMAP Neighbors", min = 5, max = 50, value = 15),
        sliderInput("p_value_threshold", "P-Value Threshold for Volcano Plot", min = 0, max = 0.1, value = 0.05, step = 0.001),
        sliderInput("log_fc_threshold", "Log Fold Change Threshold for Volcano Plot", min = 0, max = 5, value = 1, step = 0.1),
        downloadButton("output", "Export Results as CSV"),
        downloadButton("All_genes", "Export Enrichr results for all 20 genes as a CSV"),
        downloadButton("Up_genes", "Export Enrichr results for upregulated genes as a CSV"),
        downloadButton("Down_genes", "Export Enrichr results for downregulated genes as a CSV"),
      ),
      mainPanel(
        tabsetPanel(id = "main_tabset",
                    tabPanel("Differential Expression Results", DTOutput("Dif_expr_results")),
                    tabPanel("PCA Plot", plotlyOutput("pcaplot")),  
                    tabPanel("UMAP Plot", plotlyOutput("umapplot")),  
                    tabPanel("Volcano Plot", plotlyOutput("volcano_plot")),  
                    tabPanel("Pathway Enrichment", 
                             tabsetPanel(
                               tabPanel("All Genes", DTOutput("enrich_table_all")),
                               tabPanel("Upregulated Genes", DTOutput("enrich_table_up")),
                               tabPanel("Downregulated Genes", DTOutput("enrich_table_down"))
                             )
                    ),
                    tabPanel("Random Forest Performance Metrics", 
                             DTOutput("rf_metrics"),
                             downloadButton("download_metrics", "Download Performance Metrics as CSV")
                    ),
                   
                    
                    
                    tabPanel("Power Calculation", verbatimTextOutput("power_results")),
                    tabPanel("Read Me", verbatimTextOutput("readme"))
        )
      )
    )
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024^2)
  enrich_res_up <- reactiveVal()
  enrich_res_down <- reactiveVal()
  
  # Dynamic UI for Chip Selection
  observe({
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    filters <- listFilters(mart)
    affy_filters <- filters[grep("^affy", filters$name), "name"]
    
    output$chip_ui <- renderUI({
      selectInput("chip", "Input the chip that was used in your experiment", choices = affy_filters)
    })
  })
  
  final_results <- reactiveVal()
  
  observeEvent(input$running_Dif, {
    shinyjs::info("Starting Differential Expression Analysis...")
    res <- dif_expr(expr_data = input$counts_input$datapath, phenotype_data = input$phenotype_input$datapath, chip = input$chip)
    final_results(res)
    output$Dif_expr_results <- renderDT({
      datatable(res)
    })
    session$sendCustomMessage("notify", "Differential Expression Analysis Completed!")
  })
  
  output$output <- downloadHandler(
    filename = "differential_expression_results.csv",
    content = function(file) {
      write.csv(final_results(), file, row.names = FALSE)
    }
  )
  
  observeEvent(input$plot_pca, {
    req(input$counts_input, input$phenotype_input)
    shinyjs::info("Generating PCA Plot...")
    output$pcaplot <- renderPlotly({
      Make_PCA(expr_data = input$counts_input$datapath, phenotype_data = input$phenotype_input$datapath, n_components = input$n_components_pca)
    })
    session$sendCustomMessage("notify", "PCA Plot Generated!")
  })
  
  observeEvent(input$plot_umap, {
    req(input$counts_input, input$phenotype_input)
    shinyjs::info("Generating UMAP Plot...")
    output$umapplot <- renderPlotly({
      Make_UMAP(expr_data = input$counts_input$datapath, phenotype_data = input$phenotype_input$datapath, n_neighbors = input$n_neighbors_umap)
    })
    session$sendCustomMessage("notify", "UMAP Plot Generated!")
  })
  
  observeEvent(input$plot_volcano, {
    req(final_results())
    shinyjs::info("Generating Volcano Plot...")
    output$volcano_plot <- renderPlotly({
      Make_Volcano(final_results(), input$p_value_threshold, input$log_fc_threshold)
    })
    session$sendCustomMessage("notify", "Volcano Plot Generated!")
  })
  
  observeEvent(input$pathway_all, {
    req(final_results())
    shinyjs::info("Performing Pathway Enrichment (All Genes)...")
    res <- final_results()
    enrich_res <- enrich_all_genes(res)
    output$enrich_table_all <- renderDT({
      datatable(enrich_res)
    })
    session$sendCustomMessage("notify", "Pathway Enrichment (All Genes) Completed!")
  })
  
  output$All_genes <- downloadHandler(
    filename = "Enrichment_on_all_genes.csv",
    content = function(file) {
      write.csv(combined_res_all, file, row.names = FALSE)
    }
  )
  observeEvent(input$pathway_up, {
    req(final_results())
    shinyjs::info("Performing Pathway Enrichment (Upregulated Genes)...")
    res <- final_results()
    up_result <- enrich_up_genes(res)
    
    # Ensure up_result is a data frame
    if (!is.data.frame(up_result)) {
      up_result <- as.data.frame(up_result)
    }
    
    enrich_res_up(up_result)
    
    output$enrich_table_up <- renderDT({
      if ("Message" %in% colnames(up_result)) {
        datatable(data.frame(Message = "No upregulated genes found."), options = list(dom = 't'))
      } else {
        datatable(up_result)
      }
    })  
    session$sendCustomMessage("notify", "Pathway Enrichment (Upregulated Genes) Completed!")
  })
  
  output$Up_genes <- downloadHandler(
    filename = "Enrichment_on_upregulated_genes.csv",
    content = function(file) {
      up_data <- enrich_res_up()  # Get the reactive value
      if (!is.null(up_data) && !("Message" %in% colnames(up_data))) {
        # Ensure up_data is a data frame
        if (!is.data.frame(up_data)) {
          up_data <- as.data.frame(up_data)
        }
        write.csv(up_data, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No upregulated genes found."), file, row.names = FALSE)
      }
    }
  )
  
  observeEvent(input$pathway_down, {
    req(final_results())
    shinyjs::info("Performing Pathway Enrichment (Downregulated Genes)...")
    res <- final_results()
    down_result <- enrich_down_genes(res)
    
    # Ensure down_result is a data frame
    if (!is.data.frame(down_result)) {
      down_result <- as.data.frame(down_result)
    }
    
    enrich_res_down(down_result)
    
    output$enrich_table_down <- renderDT({
      if ("Message" %in% colnames(down_result)) {
        datatable(data.frame(Message = "No downregulated genes found."), options = list(dom = 't'))
      } else {
        datatable(down_result)
      }
    }) 
    session$sendCustomMessage("notify", "Pathway Enrichment (Downregulated Genes) Completed!")
  })
  
  output$Down_genes <- downloadHandler(
    filename = "Enrichment_on_downregulated_genes.csv",
    content = function(file) {
      down_data <- enrich_res_down()  # Get the reactive value
      if (!is.null(down_data) && !("Message" %in% colnames(down_data))) {
        # Ensure down_data is a data frame
        if (!is.data.frame(down_data)) {
          down_data <- as.data.frame(down_data)
        }
        write.csv(down_data, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No downregulated genes found."), file, row.names = FALSE)
      }
    }
  )
  

  rf_performance_metrics <- reactiveVal()
  
  
  observeEvent(input$run_rf, {
    shinyjs::info("Starting Regularized Random Forest...")
    
    try({
      req(input$counts_input, input$phenotype_input, final_results())
      
      # Extract selected genes from the first column of final_results
      selected_genes <- trimws(final_results()[, 1])
      
      # Run the regularized random forest
      rf_result <- regularized_random_forest(
        expr_data = input$counts_input$datapath, 
        phenotype_data = input$phenotype_input$datapath, 
        selected_genes = selected_genes
      )
      
      # Display the Random Forest metrics in a table
      output$rf_metrics <- renderDT({
        datatable(rf_result$metrics)
      })
      
      # Download handler for the metrics table
      output$download_metrics <- downloadHandler(
        filename = "regularized_rf_metrics.csv",
        content = function(file) {
          write.csv(rf_result$metrics, file, row.names = FALSE)
        }
      )
      
      # Notify the user that the Random Forest process is done
      session$sendCustomMessage("notify", "Regularized Random Forest Completed!")
      
    }, silent = FALSE)
  })
  
  

  # Render the README file content
  output$readme <- renderPrint({
    readme_path <- normalizePath("JCAP_Mircroarray_APP readme.txt", mustWork = TRUE)
    readLines(readme_path, warn = FALSE)
  })
  

  


    observeEvent(input$calculate_power, {
      # Ensure data is loaded
      print("Power analysis button clicked.")
      
      shinyjs::info("Starting Statistical Power Calculation...")
      
      req(input$counts_input, input$phenotype_input)
      
      # Load pheno data
      if (grepl("\\.csv$", input$phenotype_input$datapath)) {
        pheno <- read.csv(input$phenotype_input$datapath, row.names = 1)
      } else {
        Exprmat <- read.table(input$phenotype_input$datapath, header = TRUE, row.names = 1)
      }
      
      # Calculate power based on the number of samples in the expression data
      power_result <- calculate_power(phenotype_data = Exprmat, effect_size = 0.5)
      
      # Display the power result
      output$power_results <- renderPrint({
        print(power_result)
      })
      
      session$sendCustomMessage("notify", "Statistical Power Calculation Completed!")
    })
    


}

shinyApp(ui=ui,server = server)



