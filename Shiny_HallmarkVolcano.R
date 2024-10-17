# Interactive shiny app
##############################
## Prepairing the environment
##############################

# Don`t show code in the report set up 
knitr::opts_chunk$set(echo = FALSE)

# cleaning sessionInfo
rm(list = ls())
# Reset graphical parameters
graphics.off()
# Reset options to default (optional)
options(default = TRUE)

#### Library preparation
# Install if necessary 
if (!requireNamespace("edgeR", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("edgeR")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install('EnhancedVolcano')
}


# attaching libraries
library(ggplot2)
library(edgeR)
library(EnhancedVolcano)
library(limma)

library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(stringr)

library(org.Mm.eg.db)  # Mouse-specific annotation
library(msigdbr) # MSigDB
library(clusterProfiler)


##############################
## Prepairing the data
##############################


# Reading data
rawData_trimmed <- read.table('init_data/counts_matrix_trimmed.txt', 
                              header = TRUE, 
                              row.names = 1,
                              stringsAsFactors = FALSE)

# Assign sample names
colnames(rawData_trimmed) <- c('KO1', 'KO2', 'KO3', 'WT1', 'WT2', 'WT3')

# Strip daway metadata
# Check if there's any unwanted metadata row, e.g., 'Geneid', and remove it if necessary
# If the first row is metadata, you can drop it explicitly:
if ("Geneid" %in% rownames(rawData_trimmed)) {
  rawData_trimmed <- rawData_trimmed[-which(rownames(rawData_trimmed) == "Geneid"), ]
}

# Convert chr to int, preserving rownames - using [] to preserve row attributes
rawData_trimmed[] <- lapply(rawData_trimmed, as.integer)

# Check for any NA values (indicating conversion issues)
any_na <- any(is.na(rawData_trimmed))
# If any NA values are found, investigate
if (any_na) {
  print("Warning: Some non-numeric values were converted to NA. Please check your data.")
} else {
  print("No NAs introduced during values conversion")
}

# Create DGEList object
rawDGE_t <- DGEList(rawData_trimmed)

# Create a grouping factor 
group <- factor(
  c(
    rep("KO", 3),
    rep("WT", 3)
  )
)

# Add grouping to the DGEList object
rawDGE_t$samples$group <- group

# filtering genes by expression, getting rid of low count genes
keep <- filterByExpr(rawDGE_t)
rawDGE_t <- rawDGE_t[keep, , keep.lib.sizes=FALSE] # update lib sizes

# data normalization
rawDGE_t <- normLibSizes(rawDGE_t,
                         method = "TMM")


##############################
## Data modelling
##############################

# Releveling the factor to set WT as the reference level for further modelling
rawDGE_t$samples$group<- relevel(rawDGE_t$samples$group, ref = "WT")

# Setting the desing 
design <- model.matrix(~group, data = rawDGE_t$samples)
colnames(design) <- c("WT", "MecrKO")

## EdgeR`s voomLmFit does all this under the hood: logCPM, limmaWeighted, corfit, lmFit
fit <- edgeR::voomLmFit(counts = rawDGE_t,
                        design = design,
                        sample.weights = TRUE
)
# edgeR::voomLmFit performs limma`s voomWithQualityWeights with duplicateCorrelation
# If block is specified voom weights and intra-block correlation are each estimated twice

# Set the contrasts 
contrasts <- makeContrasts(
  # Overall treatment effect (across both batches)
  KOvsWT = MecrKO - WT,
  levels = design
)

# Compute contrasts from linear model
fit <- contrasts.fit(fit, contrasts)

# Compute moderated t-statistics
fit <- eBayes(fit)

# Extract contrast coefficients
## Overall treatment effect 
DE_KOvsWT <- topTable(fit, 
                      coef="KOvsWT", 
                      sort.by = "t",
                      n=Inf)


# Prepare the ranked gene list (from limma results)
# Rank by logFC or another relevant statistic

# Rank by t.value
ranked_genes <- DE_KOvsWT$t
names(ranked_genes) <- rownames(DE_KOvsWT) # Ensure Gene column contains gene symbols

# Sort genes by t (or other ranking criterion), just in case it hasn`t been sorted already
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Retrieve Mouse Hallmark gene sets from MSigDB
msigdb_H <- msigdbr(species = "Mus musculus", 
                    category = "H")  # "H" is for hallmark gene sets

# Perform GSEA using clusterProfiler
GSE_H <- GSEA(ranked_genes, 
              TERM2GENE = msigdb_H[, c("gs_name", "gene_symbol")], 
              pvalueCutoff = 1, 
              verbose = FALSE,
              pAdjustMethod = "fdr",
              eps = 0,
              by = "fgsea",
              seed = 123,
              nPermSimple = 100000)


# Define analyze_pathway_volcano_interactive

##################
### Combined Volcano Interactive: Claude with GPT Plotly Hover and p_method logic
##################

library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)
library(plotly)  # Add plotly for interactivity

# Modify the analyze_pathway_volcano function to return a plotly object with hover
analyze_pathway_volcano_interactive <- function(pathway_name, gsea_results, de_results, 
                                                p_cutoff = 0.05, fc_cutoff = 2.0, 
                                                label_method = "default", max_overlaps = 100, 
                                                style = 'clean', x_breaks = 2,
                                                p_method = "P.Value") {  # Add p_method parameter
  
  # Extract the relevant genes for the specified pathway
  message("Extracting genes for the pathway: ", pathway_name)
  pathway_genes <- as.data.frame(gsea_results) %>%
    filter(Description == pathway_name) %>%
    pull(core_enrichment) %>% 
    str_split("/", simplify = TRUE) %>% 
    as.vector()
  
  # Prepare DE results and identify pathway genes
  message("Preparing DE results...")
  de_results <- de_results %>%
    mutate(
      in_pathway = rownames(.) %in% pathway_genes,  # Mark pathway genes
      significant_fc = abs(logFC) > fc_cutoff,
      significant_p = get(p_method) < p_cutoff,  # Use the chosen p-value method dynamically
      highlight = significant_fc & significant_p,
      hover_info = paste("Gene:", rownames(.),
                         "<br>Log2FC:", round(logFC, 2),
                         "<br>", p_method, ":", round(get(p_method), 5))  # Dynamic hover text with gene info
    ) %>%
    # Add a color column based on the condition
    mutate(
      color = case_when(
        highlight ~ "p-value & Log2FC",   # Both significant fold change and p-value
        significant_fc ~ "Log2FC",        # Significant fold change only
        significant_p ~ "p-value",        # Significant p-value only
        TRUE ~ "NS"                       # Not significant
      )
    )
  
  # Determine labeling data based on method
  label_data <- switch(
    label_method,
    "default" = de_results %>% filter(in_pathway & highlight),  # Only label pathway genes with both p-value & Log2FC significant
    "fc" = de_results %>% filter(in_pathway & abs(logFC) > fc_cutoff),  # Label pathway genes crossing fold change threshold
    "p" = de_results %>% filter(in_pathway & get(p_method) < p_cutoff),  # Label pathway genes crossing p-value threshold
    "all" = de_results %>% filter(in_pathway & (abs(logFC) > fc_cutoff | get(p_method) < p_cutoff)),  # Label all significant pathway genes (either p-value or fold change)
    stop("Invalid label_method. Please use 'default', 'fc', 'p', or 'all'.")
  )
  
  # Determine the x-axis limits based on the user-defined x_breaks
  x_min <- floor(min(de_results$logFC) / x_breaks) * x_breaks  # Start from nearest multiple <= min(logFC)
  x_max <- ceiling(max(de_results$logFC) / x_breaks) * x_breaks  # End at nearest multiple >= max(logFC)
  
  # Define y-axis label as a string for plotly compatibility
  y_axis_label <- if (p_method == "P.Value") {
    "-Log10 P"
  } else {
    paste0("-Log10(", p_method, ")")
  }
  
  # Define base plot with hover information and color mapping
  base_plot <- ggplot(de_results, aes(x = logFC, y = -log10(get(p_method)), 
                                      text = hover_info, color = in_pathway)) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    labs(
      title = gsub("_", " ", pathway_name),
      subtitle = paste("p-value cutoff:", p_cutoff, "| FC cutoff:", fc_cutoff),
      x = "Log2FC",
      y = y_axis_label  # Use the dynamically created y-axis label as a string
    ) +
    # Set x-axis breaks based on user-defined x_breaks, ensuring alignment with the breaks pattern
    scale_x_continuous(breaks = seq(x_min, x_max, by = x_breaks), limits = c(x_min, x_max)) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +  # Pathway genes in red, others in grey
    custom_minimal_theme_with_grid()
  
  # Choose plot style based on user input
  ggplot_obj <- switch(
    style,
    'clean' = plot_style_clean(base_plot, de_results, pathway_genes, label_data, p_cutoff, fc_cutoff, max_overlaps),
    'claude' = plot_style_claude(base_plot, de_results, pathway_genes, label_data, p_cutoff, fc_cutoff, max_overlaps),
    'gpt' = plot_style_gpt(base_plot, de_results, pathway_genes, label_data, p_cutoff, fc_cutoff, max_overlaps),
    stop("Invalid style. Please use 'clean', 'claude', or 'gpt'.")
  )
  
  # Convert ggplot object to plotly for interactivity
  plotly_obj <- ggplotly(ggplot_obj, tooltip = c("text"))  # Use the custom hover text
  
  return(plotly_obj)
}




##################
### Define the Shiny app
##################

ui <- fluidPage(
  titlePanel("Interactive Volcano Plot with Pathway Selection"),
  
  # Create a dropdown menu for pathway selection
  selectInput(
    inputId = "pathway", 
    label = "Select a pathway to highlight:",
    choices = NULL,  # Choices will be dynamically updated
    selected = NULL
  ),
  
  # Create a space for the plot
  plotlyOutput("volcanoPlot")
)

server <- function(input, output, session) {
  
  # Load your data (GSE_H and DE_KOvsWT)
  # GSE_H and DE_KOvsWT should already be available in your environment
  gsea_results <- as.data.frame(GSE_H@result)  # Convert gseaResult to a data frame
  de_results <- DE_KOvsWT
  
  # Filter pathways based on q_value < 0.05 and dynamically update dropdown choices
  pathways_filtered <- reactive({
    gsea_results %>%
      filter(qvalue < 0.05) %>%
      pull(Description)  # Extract the pathway descriptions
  })
  
  # Update the pathway dropdown choices dynamically
  observe({
    updateSelectInput(session, "pathway", choices = pathways_filtered(), selected = pathways_filtered()[1])
  })
  
  # Generate the interactive volcano plot based on the selected pathway
  output$volcanoPlot <- renderPlotly({
    
    # Call the analyze_pathway_volcano_interactive function to create the plot
    interactive_plot <- analyze_pathway_volcano_interactive(
      pathway_name = input$pathway,
      gsea_results = gsea_results,
      de_results = de_results,
      p_cutoff = 0.05,
      fc_cutoff = 2.0,
      label_method = "default",
      max_overlaps = 50,
      style = 'clean',
      x_breaks = 2  # Set x-axis breaks to every 2 logFC units
    )
    
    interactive_plot  # Return the plotly object
  })
}

# Run the application
shinyApp(ui = ui, server = server)
