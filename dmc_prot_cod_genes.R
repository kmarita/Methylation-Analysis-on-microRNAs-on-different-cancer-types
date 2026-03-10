library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tidyr)

main_dir <- "~/methylation_project"
results_dir <- file.path(main_dir, "methylation_analysis_results")

cancer_types <- c("TCGA-BRCA","TCGA-BLCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", 
                  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
                  "TCGA-LIHC","TCGA-LUSC", "TCGA-LUAD", "TCGA-PAAD", "TCGA-PCPG", 
                  "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", 
                  "TCGA-THYM", "TCGA-UCEC")

# Load protein-coding gene coordinates
protein_coding_bed_file <- file.path(main_dir, "prot_cod_genes.bed")

protein_coding_genes <- read.table(protein_coding_bed_file, header = FALSE, sep = "\t", 
                                   stringsAsFactors = FALSE, comment.char = "")
colnames(protein_coding_genes) <- c("Chromosome", "Start", "End", "Gene", "Score", "Strand")

cat(sprintf(" Loaded %d protein-coding gene coordinates from BED file\n", nrow(protein_coding_genes)))

# Function to read all CpGs in protein-coding genes
read_all_gene_cpgs <- function(cancer_type, results_dir, protein_coding_genes) {
  file_path <- file.path(results_dir, cancer_type, 
                         paste0(cancer_type, "_DMCs_in_protein_coding_genes.tsv"))
  
  if (!file.exists(file_path)) {
    cat(sprintf("File not found for %s\n", cancer_type))
    return(NULL)
  }
  
  all_cpgs <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  required_cols <- c("Chromosome", "Start", "End")
  missing_cols <- setdiff(required_cols, colnames(all_cpgs))
  
  if (length(missing_cols) > 0) {
    cat(sprintf("Missing columns in %s: %s\n", cancer_type, paste(missing_cols, collapse=", ")))
    return(NULL)
  }
  
  all_cpgs <- all_cpgs[!is.na(all_cpgs$Chromosome) & 
                         !is.na(all_cpgs$Start) & 
                         !is.na(all_cpgs$End), ]
  
  all_cpgs <- all_cpgs[all_cpgs$End >= all_cpgs$Start, ]
  all_cpgs <- all_cpgs[all_cpgs$Start >= 0 & all_cpgs$End >= 0, ]
  
  if (nrow(all_cpgs) == 0) {
    cat(sprintf("No valid CpGs for %s after QC\n", cancer_type))
    return(NULL)
  }
  
  all_cpgs$Gene <- NA_character_
  
  for (i in 1:nrow(protein_coding_genes)) {
    chr <- protein_coding_genes$Chromosome[i]
    start <- protein_coding_genes$Start[i]
    end <- protein_coding_genes$End[i]
    gene_name <- protein_coding_genes$Gene[i]
    
    overlaps <- all_cpgs$Chromosome == chr & 
      all_cpgs$Start >= start & 
      all_cpgs$End <= end
    
    all_cpgs$Gene[overlaps] <- gene_name
  }
  
  all_cpgs <- all_cpgs[!is.na(all_cpgs$Gene), ]
  
  cat(sprintf("Loaded %d CpGs in protein-coding gene bodies from %s\n", nrow(all_cpgs), cancer_type))
  return(all_cpgs)
}

# Function to read significant DMCs in protein-coding genes
read_gene_dmcs <- function(cancer_type, results_dir, protein_coding_genes) {
  file_path <- file.path(results_dir, cancer_type, 
                         paste0(cancer_type, "_DMCs_in_protein_coding_genes.tsv"))
  
  if (!file.exists(file_path)) {
    cat(sprintf("File not found for %s\n", cancer_type))
    return(NULL)
  }
  
  dmcs <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  if (!"Methylation_Status" %in% colnames(dmcs)) {
    cat(sprintf("Column 'Methylation_Status' not found in %s\n", cancer_type))
    return(NULL)
  }
  
  dmcs_sig <- dmcs[dmcs$Methylation_Status %in% c("Hypermethylated", "Hypomethylated"), ]
  
  if (nrow(dmcs_sig) == 0) {
    cat(sprintf("No significant DMCs for %s\n", cancer_type))
    return(NULL)
  }
  
  dmcs_sig$DMC_Status <- dmcs_sig$Methylation_Status
  
  required_cols <- c("Chromosome", "Start", "End")
  missing_cols <- setdiff(required_cols, colnames(dmcs_sig))
  
  if (length(missing_cols) > 0) {
    cat(sprintf("Missing columns in %s: %s\n", cancer_type, paste(missing_cols, collapse=", ")))
    return(NULL)
  }
  
  if (!"ProbeID" %in% colnames(dmcs_sig)) {
    dmcs_sig$DMC_idx <- paste0("DMC_", 1:nrow(dmcs_sig))
  } else {
    dmcs_sig$DMC_idx <- dmcs_sig$ProbeID
  }
  
  dmcs_sig <- dmcs_sig[!is.na(dmcs_sig$Chromosome) & 
                         !is.na(dmcs_sig$Start) & 
                         !is.na(dmcs_sig$End), ]
  
  dmcs_sig <- dmcs_sig[dmcs_sig$End >= dmcs_sig$Start, ]
  dmcs_sig <- dmcs_sig[dmcs_sig$Start >= 0 & dmcs_sig$End >= 0, ]
  
  if (nrow(dmcs_sig) == 0) {
    cat(sprintf("No valid DMCs remaining for %s after QC\n", cancer_type))
    return(NULL)
  }
  
  dmcs_sig$Gene <- NA_character_
  
  for (i in 1:nrow(protein_coding_genes)) {
    chr <- protein_coding_genes$Chromosome[i]
    start <- protein_coding_genes$Start[i]
    end <- protein_coding_genes$End[i]
    gene_name <- protein_coding_genes$Gene[i]
    
    overlaps <- dmcs_sig$Chromosome == chr & 
      dmcs_sig$Start >= start & 
      dmcs_sig$End <= end
    
    dmcs_sig$Gene[overlaps] <- gene_name
  }
  
  dmcs_sig <- dmcs_sig[!is.na(dmcs_sig$Gene), ]
  
  cat(sprintf("Loaded %d significant protein-coding gene DMCs from %s\n", nrow(dmcs_sig), cancer_type))
  return(dmcs_sig)
}

all_cpg_list <- list()
all_dmcs_list <- list()

for (cancer in cancer_types) {
  all_cpgs <- read_all_gene_cpgs(cancer, results_dir, protein_coding_genes)
  if (!is.null(all_cpgs)) {
    all_cpg_list[[cancer]] <- all_cpgs
  }
  
  dmcs_sig <- read_gene_dmcs(cancer, results_dir, protein_coding_genes)
  if (!is.null(dmcs_sig)) {
    all_dmcs_list[[cancer]] <- dmcs_sig
  }
}

if (length(all_dmcs_list) == 0) {
  stop("No protein-coding gene DMC files found with significant results.")
}

cat(sprintf("Total cancer types with CpG data: %d\n", length(all_cpg_list)))
cat(sprintf("Total cancer types with significant protein-coding gene DMCs: %d\n", length(all_dmcs_list)))


cpg_counts_per_gene_cancer <- data.frame()

for (cancer in names(all_cpg_list)) {
  cpg_data <- all_cpg_list[[cancer]]
  
  gene_counts <- cpg_data %>%
    group_by(Gene) %>%
    summarise(Total_CpGs = n(), .groups = "drop") %>%
    mutate(Cancer = cancer)
  
  cpg_counts_per_gene_cancer <- rbind(cpg_counts_per_gene_cancer, gene_counts)
}

cat(sprintf(" Calculated CpG counts for %d gene-cancer combinations\n", 
            nrow(cpg_counts_per_gene_cancer)))

cpg_counts_summary <- cpg_counts_per_gene_cancer %>%
  group_by(Gene) %>%
  summarise(
    Mean_CpGs = mean(Total_CpGs),
    Min_CpGs = min(Total_CpGs),
    Max_CpGs = max(Total_CpGs),
    Num_Cancers = n(),
    .groups = "drop"
  )

cat("\nCpG count summary (top 10 genes by mean CpGs):\n")
print(head(cpg_counts_summary[order(-cpg_counts_summary$Mean_CpGs), ], 10))

cpg_counts_file <- file.path(results_dir, "protein_coding_gene_CpG_counts_per_cancer.tsv")
write.table(cpg_counts_per_gene_cancer, file = cpg_counts_file, 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n Saved CpG counts: %s\n", cpg_counts_file))

all_dmcs <- data.frame()

cancer_names <- names(all_dmcs_list)
for (cancer in cancer_names) {
  dmcs_sig <- all_dmcs_list[[cancer]]
  
  if (nrow(dmcs_sig) > 0) {
    dmc_info <- dmcs_sig[, c("DMC_idx", "Chromosome", "Start", "End", "Gene", "DMC_Status")]
    dmc_info$Cancer <- cancer
    dmc_info$Coord <- paste(dmc_info$Chromosome, dmc_info$Start, sep = "_")
    all_dmcs <- rbind(all_dmcs, dmc_info)
  }
}

cat(sprintf("\nTotal DMCs across all cancers: %d\n\n", nrow(all_dmcs)))

cat("\n NORMALIZATION \n")
cat("Formula: [Sum(DMCs in all genes)] / [Sum(CpGs in all genes)] per cancer\n\n")

dmc_counts_per_gene_cancer <- all_dmcs %>%
  group_by(Cancer, Gene, DMC_Status) %>%
  summarise(DMC_Count = n(), .groups = "drop")

total_dmcs_per_gene_cancer <- dmc_counts_per_gene_cancer %>%
  group_by(Cancer, Gene) %>%
  summarise(Total_DMCs = sum(DMC_Count), .groups = "drop")

gene_dmc_cpg <- total_dmcs_per_gene_cancer %>%
  left_join(cpg_counts_per_gene_cancer, by = c("Cancer", "Gene"))

gene_dmc_cpg <- gene_dmc_cpg %>%
  filter(!is.na(Total_CpGs))

if (nrow(gene_dmc_cpg) == 0) {
  stop("No matching gene names between DMC data and CpG counts!")
}

cat(sprintf("Successfully matched %d gene-cancer pairs with CpG counts\n", 
            nrow(gene_dmc_cpg)))

normalization_per_cancer <- gene_dmc_cpg %>%
  group_by(Cancer) %>%
  summarise(
    Sum_DMCs = sum(Total_DMCs),
    Sum_CpGs = sum(Total_CpGs),
    Normalization_Ratio = Sum_DMCs / Sum_CpGs,
    .groups = "drop"
  )

cat("\nNormalization ratios per cancer:\n")
print(normalization_per_cancer)

dmc_status_with_cpg <- dmc_counts_per_gene_cancer %>%
  left_join(cpg_counts_per_gene_cancer, by = c("Cancer", "Gene")) %>%
  filter(!is.na(Total_CpGs))

dmc_status_sums <- dmc_status_with_cpg %>%
  group_by(Cancer, DMC_Status) %>%
  summarise(
    Sum_DMCs_Status = sum(DMC_Count),
    Sum_CpGs = sum(Total_CpGs),
    .groups = "drop"
  )

dmc_status_sums <- dmc_status_sums %>%
  mutate(Normalized_Ratio = (Sum_DMCs_Status / Sum_CpGs) * 100)

dmc_status_sums$Cancer <- factor(dmc_status_sums$Cancer, levels = cancer_types)

normalization_stats_file <- file.path(results_dir, "dmc_normalization_statistics_protein_coding.tsv")
write.table(dmc_status_sums, file = normalization_stats_file, 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n Saved normalization statistics: %s\n", normalization_stats_file))

# PLOT 1: Hypermethylated vs Hypomethylated per cancer
p_status_norm <- ggplot(dmc_status_sums, aes(x = Cancer, y = Normalized_Ratio, fill = DMC_Status)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sprintf("%.2f%%", Normalized_Ratio)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Hypermethylated" = "steelblue", "Hypomethylated" = "salmon")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Hypermethylated vs Hypomethylated Protein-Coding Gene Body DMCs per Cancer Type",
    subtitle = "Normalized as: (DMC count / Total CpG count) × 100",
    x = "Cancer Type", 
    y = "Normalized DMC Percentage (%)",
    fill = "DMC Status"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "top"
  ) +
  ylim(0, max(dmc_status_sums$Normalized_Ratio) * 1.15)

plot_status_norm_file <- file.path(results_dir, "dmc_status_hyper_hypo_protein_coding_per_cancer_NORMALIZED.pdf")
ggsave(plot_status_norm_file, p_status_norm, width = 14, height = 6)
cat(sprintf("Saved Plot 1: %s\n", plot_status_norm_file))

# Continue with concordance analysis
dmc_coords <- unique(all_dmcs$Coord)
dmc_concordance <- data.frame()

cat(sprintf("\nAnalyzing concordance across %d unique CpG positions...\n", length(dmc_coords)))

for (coord in dmc_coords) {
  coord_data <- all_dmcs[all_dmcs$Coord == coord, ]
  num_cancers <- length(unique(coord_data$Cancer))
  
  if (num_cancers >= 2) {
    all_hyper <- all(coord_data$DMC_Status == "Hypermethylated")
    all_hypo <- all(coord_data$DMC_Status == "Hypomethylated")
    
    if (all_hyper) {
      concordance <- "All_Hypermethylated"
    } else if (all_hypo) {
      concordance <- "All_Hypomethylated"
    } else {
      concordance <- "Discordant"
    }
    
    unique_genes <- unique(coord_data$Gene)
    unique_genes <- unique_genes[unique_genes != "" & !is.na(unique_genes)]
    
    dmc_concordance <- rbind(dmc_concordance, data.frame(
      Coord = coord,
      DMC_idx = paste(unique(coord_data$DMC_idx), collapse = ";"),
      Chromosome = coord_data$Chromosome[1],
      Start = coord_data$Start[1],
      End = coord_data$End[1],
      Gene = paste(sort(unique_genes), collapse = ";"),
      Num_Cancers = num_cancers,
      Cancers = paste(sort(unique(coord_data$Cancer)), collapse = ";"),
      Concordance = concordance,
      stringsAsFactors = FALSE
    ))
  }
}

concordant_dmcs <- dmc_concordance[dmc_concordance$Concordance != "Discordant", ]

cat(sprintf("\nTotal gene body DMC positions in 2+ cancers (exact coordinates): %d\n", nrow(dmc_concordance)))
cat(sprintf("Concordant gene body DMCs (all same direction): %d\n", nrow(concordant_dmcs)))
cat(sprintf("  All Hypermethylated: %d\n", sum(concordant_dmcs$Concordance == "All_Hypermethylated")))
cat(sprintf("  All Hypomethylated: %d\n", sum(concordant_dmcs$Concordance == "All_Hypomethylated")))
cat(sprintf("Discordant gene body DMCs (mixed directions): %d\n\n", sum(dmc_concordance$Concordance == "Discordant")))

cat("Distribution of gene body DMCs by number of cancers:\n")
cancer_count_table <- table(dmc_concordance$Num_Cancers)
print(cancer_count_table)
cat("\n")

if (nrow(concordant_dmcs) > 0) {
  
  gene_cancer_dmc <- list()
  
  for (i in 1:nrow(concordant_dmcs)) {
    genes <- unlist(strsplit(concordant_dmcs$Gene[i], ";"))
    genes <- genes[genes != "" & !is.na(genes)]
    genes <- trimws(genes)
    
    cancers <- unlist(strsplit(concordant_dmcs$Cancers[i], ";"))
    cancers <- trimws(cancers)
    
    dmc_coord <- concordant_dmcs$Coord[i]
    
    for (gene in genes) {
      for (cancer in cancers) {
        key <- paste(gene, cancer, sep = "||")
        if (is.null(gene_cancer_dmc[[key]])) {
          gene_cancer_dmc[[key]] <- character()
        }
        gene_cancer_dmc[[key]] <- c(gene_cancer_dmc[[key]], dmc_coord)
      }
    }
  }
  
  all_genes <- unique(unlist(strsplit(concordant_dmcs$Gene, ";")))
  all_genes <- all_genes[all_genes != "" & !is.na(all_genes)]
  all_genes <- trimws(all_genes)
  all_genes <- sort(all_genes)
  
  all_cancers_with_dmcs <- unique(unlist(strsplit(concordant_dmcs$Cancers, ";")))
  all_cancers_with_dmcs <- trimws(all_cancers_with_dmcs)
  all_cancers_with_dmcs <- sort(all_cancers_with_dmcs)
  
  dmc_matrix <- matrix(0, nrow = length(all_genes), ncol = length(all_cancers_with_dmcs),
                       dimnames = list(all_genes, all_cancers_with_dmcs))
  
  for (gene in all_genes) {
    for (cancer in all_cancers_with_dmcs) {
      key <- paste(gene, cancer, sep = "||")
      if (!is.null(gene_cancer_dmc[[key]])) {
        dmc_matrix[gene, cancer] <- length(unique(gene_cancer_dmc[[key]]))
      }
    }
  }
  
  cat("\nMATRIX NORMALIZATION \n")
  cat("Formula: (DMC count / CpG count) × 100 for each gene-cancer pair\n\n")
  
  dmc_matrix_norm <- dmc_matrix
  
  for (gene in rownames(dmc_matrix_norm)) {
    for (cancer in colnames(dmc_matrix_norm)) {
      cpg_count <- cpg_counts_per_gene_cancer %>%
        filter(Gene == gene, Cancer == cancer) %>%
        pull(Total_CpGs)
      
      if (length(cpg_count) > 0 && cpg_count > 0) {
        dmc_matrix_norm[gene, cancer] <- (dmc_matrix[gene, cancer] / cpg_count) * 100
      } else {
        dmc_matrix_norm[gene, cancer] <- 0
      }
    }
  }
  
  dmc_matrix_df <- as.data.frame(dmc_matrix)
  dmc_matrix_df <- cbind(Gene = rownames(dmc_matrix_df), dmc_matrix_df)
  rownames(dmc_matrix_df) <- NULL
  
  dmc_matrix_norm_df <- as.data.frame(dmc_matrix_norm)
  dmc_matrix_norm_df <- cbind(Gene = rownames(dmc_matrix_norm_df), dmc_matrix_norm_df)
  rownames(dmc_matrix_norm_df) <- NULL
  
  output_matrix_file <- file.path(results_dir, "protein_coding_gene_cancer_dmc_matrix.tsv")
  write.table(dmc_matrix_df, file = output_matrix_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  output_matrix_norm_file <- file.path(results_dir, "protein_coding_gene_cancer_dmc_matrix_NORMALIZED.tsv")
  write.table(dmc_matrix_norm_df, file = output_matrix_norm_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat(sprintf("\nProtein-Coding Gene x Cancer DMC matrix created with %d genes and %d cancer types\n", 
              length(all_genes), length(all_cancers_with_dmcs)))
  cat(sprintf("Raw matrix saved to: %s\n", output_matrix_file))
  cat(sprintf("Normalized matrix saved to: %s\n\n", output_matrix_norm_file))
  
  cat("Matrix summary:\n")
  cat(sprintf("  Total genes: %d\n", nrow(dmc_matrix)))
  cat(sprintf("  Total cancer types: %d\n", ncol(dmc_matrix)))
  cat(sprintf("  Total non-zero entries: %d\n", sum(dmc_matrix > 0)))
  cat(sprintf("  Max DMCs for any gene-cancer pair: %d\n", max(dmc_matrix)))
  
  gene_totals <- rowSums(dmc_matrix)
  top_genes <- head(sort(gene_totals, decreasing = TRUE), 20)
  cat("\nTop 20 genes by total exact coordinate gene body DMCs across all cancers:\n")
  print(top_genes)
  
  gene_concordant_counts <- data.frame()
  
  for (i in 1:nrow(concordant_dmcs)) {
    genes <- unlist(strsplit(concordant_dmcs$Gene[i], ";"))
    genes <- genes[genes != "" & !is.na(genes)]
    genes <- trimws(genes)
    
    for (gene in genes) {
      concordance_type <- concordant_dmcs$Concordance[i]
      num_cancers <- concordant_dmcs$Num_Cancers[i]
      
      gene_concordant_counts <- rbind(gene_concordant_counts,
                                      data.frame(Gene = gene,
                                                 Concordance = concordance_type,
                                                 Num_Cancers = num_cancers,
                                                 stringsAsFactors = FALSE))
    }
  }
  
  gene_summary <- data.frame()
  for (gene in unique(gene_concordant_counts$Gene)) {
    gene_data <- gene_concordant_counts[gene_concordant_counts$Gene == gene, ]
    
    total_concordant <- nrow(gene_data)
    all_hyper <- sum(gene_data$Concordance == "All_Hypermethylated")
    all_hypo <- sum(gene_data$Concordance == "All_Hypomethylated")
    
    gene_rows <- grepl(paste0("(^|;)", gene, "(;|$)"), concordant_dmcs$Gene)
    num_cancers <- length(unique(unlist(strsplit(
      paste(concordant_dmcs$Cancers[gene_rows], collapse = ";"), ";"))))
    
    gene_summary <- rbind(gene_summary,
                          data.frame(Gene = gene,
                                     Total_Concordant_DMCs = total_concordant,
                                     All_Hypermethylated_DMCs = all_hyper,
                                     All_Hypomethylated_DMCs = all_hypo,
                                     Num_Cancers = num_cancers,
                                     stringsAsFactors = FALSE))
  }
  
  gene_summary <- gene_summary[order(-gene_summary$Total_Concordant_DMCs), ]
  
  output_file1 <- file.path(results_dir, "concordant_protein_coding_gene_DMCs_summary.tsv")
  write.table(gene_summary, file = output_file1, sep = "\t", quote = FALSE, row.names = FALSE)
  
  output_file2 <- file.path(results_dir, "concordant_protein_coding_DMCs_detailed.tsv")
  write.table(concordant_dmcs, file = output_file2, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("\nAnalysis complete.\n")
  cat(sprintf("Gene-Cancer DMC matrix: %s\n", output_matrix_file))
  cat(sprintf("Gene-Cancer NORMALIZED DMC matrix: %s\n", output_matrix_norm_file))
  cat(sprintf("Concordant gene DMCs summary: %s\n", output_file1))
  cat(sprintf("Concordant DMCs detailed: %s\n", output_file2))
  
  # PLOT 2: Distribution by number of cancers
  cancer_count_df <- as.data.frame(table(dmc_concordance$Num_Cancers))
  colnames(cancer_count_df) <- c("Num_Cancers", "Count")
  cancer_count_df$Num_Cancers <- as.numeric(as.character(cancer_count_df$Num_Cancers))
  
  total_concordant_dmcs <- sum(cancer_count_df$Count)
  cancer_count_df$Normalized_Percent <- (cancer_count_df$Count / total_concordant_dmcs) * 100
  
  p2 <- ggplot(cancer_count_df, aes(x = Num_Cancers, y = Normalized_Percent)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black") +
    geom_text(aes(label = sprintf("%.1f%%", Normalized_Percent)), 
              vjust = -0.5, size = 3.5) +
    theme_minimal(base_size = 12) +
    labs(title = "Distribution of Protein-Coding Gene Body DMCs by Number of Cancers",
         x = "Number of Cancers",
         y = "Percentage of Total Concordant DMCs (%)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ylim(0, max(cancer_count_df$Normalized_Percent) * 1.15)
  
  plot2_file <- file.path(results_dir, "protein_coding_dmc_cancer_distribution_NORMALIZED.pdf")
  ggsave(plot2_file, p2, width = 10, height = 6)
  cat(sprintf(" Saved Plot 2: %s\n", plot2_file))
  
  # PLOT 3: Top 20 genes
  top20_summary <- head(gene_summary, 20)
  
  top20_cpg_summary <- data.frame()
  for (gene in top20_summary$Gene) {
    cpg_data <- cpg_counts_per_gene_cancer %>%
      filter(Gene == gene)
    
    if (nrow(cpg_data) > 0) {
      top20_cpg_summary <- rbind(top20_cpg_summary, data.frame(
        Gene = gene,
        Mean_CpGs = mean(cpg_data$Total_CpGs),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  top20_summary <- top20_summary %>%
    left_join(top20_cpg_summary, by = "Gene")
  
  top20_summary$Hyper_Percent <- (top20_summary$All_Hypermethylated_DMCs / top20_summary$Mean_CpGs) * 100
  top20_summary$Hypo_Percent <- (top20_summary$All_Hypomethylated_DMCs / top20_summary$Mean_CpGs) * 100
  
  top20_long <- top20_summary %>%
    select(Gene, Hyper_Percent, Hypo_Percent, Mean_CpGs, Total_Concordant_DMCs) %>%
    pivot_longer(
      cols = c(Hyper_Percent, Hypo_Percent),
      names_to = "DMC_Status",
      values_to = "Percent"
    ) %>%
    mutate(
      DMC_Status = ifelse(DMC_Status == "Hyper_Percent", 
                          "Hypermethylated", 
                          "Hypomethylated")
    )
  
  top20_long$Gene <- factor(top20_long$Gene, levels = rev(top20_summary$Gene))
  top20_long$DMC_Status <- factor(top20_long$DMC_Status, 
                                  levels = c("Hypomethylated", "Hypermethylated"))
  
  p3 <- ggplot(top20_long, aes(x = Gene, y = Percent, fill = DMC_Status)) +
    geom_bar(stat = "identity", color = "black", position = "stack") +
    geom_text(aes(label = ifelse(Percent > 2, sprintf("%.1f%%", Percent), "")),
              position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(values = c("Hypomethylated" = "salmon",
                                 "Hypermethylated" = "steelblue")) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(
      title = "Top 20 Protein-Coding Genes by Concordant Gene Body DMC Count",
      subtitle = "Normalized as: (DMC count / Mean CpG count) × 100",
      x = "Gene",
      y = "Normalized DMC Percentage (%)",
      fill = "DMC Status"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 10)
    )
  
  plot3_file <- file.path(results_dir, "top20_protein_coding_genes_dmcs_NORMALIZED.pdf")
  ggsave(plot3_file, p3, width = 11, height = 8)
  cat(sprintf(" Saved Plot 3: %s\n", plot3_file))
  
  # PLOT 4: Chromosome distribution
  chr_counts <- table(concordant_dmcs$Chromosome)
  chr_df <- data.frame(
    Chromosome = names(chr_counts),
    Count = as.numeric(chr_counts)
  )
  
  chr_df <- chr_df[!chr_df$Chromosome %in% c("chrX", "chrY", "chrM"), ]
  
  if (nrow(chr_df) > 0) {
    total_chr_dmcs <- sum(chr_df$Count)
    chr_df$Normalized_Percent <- (chr_df$Count / total_chr_dmcs) * 100
    
    chr_order <- paste0("chr", 1:22)
    chr_df$Chromosome <- factor(chr_df$Chromosome, 
                                levels = chr_order[chr_order %in% chr_df$Chromosome])
    chr_df <- chr_df[order(chr_df$Chromosome), ]
    
    p4 <- ggplot(chr_df, aes(x = Chromosome, y = Normalized_Percent)) +
      geom_bar(stat = "identity", fill = "mediumpurple", color = "black") +
      geom_text(aes(label = sprintf("%.1f%%", Normalized_Percent)), 
                vjust = -0.5, size = 3) +
      theme_minimal(base_size = 11) +
      labs(title = "Protein-Coding Gene Body DMCs Distribution Across Chromosomes",
           x = "Chromosome",
           y = "Percentage of Total DMCs (%)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      ylim(0, max(chr_df$Normalized_Percent) * 1.15)
    
    plot4_file <- file.path(results_dir, "protein_coding_dmcs_by_chromosome_NORMALIZED.pdf")
    ggsave(plot4_file, p4, width = 12, height = 6)
    cat(sprintf(" Saved Plot 4: %s\n", plot4_file))
  }
  
  # PLOT 5: DMCs per cancer
  cancer_norm_data <- normalization_per_cancer %>%
    mutate(Normalized_Percent = Normalization_Ratio * 100) %>%
    arrange(desc(Normalized_Percent))
  
  cancer_norm_data$Cancer <- factor(cancer_norm_data$Cancer, levels = cancer_norm_data$Cancer)
  
  p5 <- ggplot(cancer_norm_data, aes(x = Cancer, y = Normalized_Percent)) +
    geom_bar(stat = "identity", fill = "salmon", color = "black") +
    geom_text(aes(label = sprintf("%.2f%%", Normalized_Percent)), 
              vjust = -0.5, size = 3, angle = 0) +
    theme_minimal(base_size = 11) +
    labs(
      title = "Concordant Protein-Coding Gene Body DMCs per Cancer Type",
      subtitle = "Normalized as: (Sum DMCs / Sum CpGs) × 100",
      x = "Cancer Type",
      y = "Normalized DMC Percentage (%)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    ) +
    ylim(0, max(cancer_norm_data$Normalized_Percent) * 1.15)
  
  plot5_file <- file.path(results_dir, "protein_coding_dmcs_per_cancer_NORMALIZED.pdf")
  ggsave(plot5_file, p5, width = 12, height = 6)
  cat(sprintf(" Saved Plot 5: %s\n", plot5_file))
  
  # PLOT 6: Heatmap
  min_cancer_types <- 2
  top_n_genes <- 20
  
  if (exists("dmc_matrix_norm")) {
    cat(sprintf("Original matrix dimensions: %d genes x %d cancers\n", 
                nrow(dmc_matrix_norm), ncol(dmc_matrix_norm)))
    cat(sprintf("Total non-zero entries: %d\n", sum(dmc_matrix_norm > 0)))
    
    multi_cancer_genes <- rowSums(dmc_matrix_norm > 0) >= min_cancer_types
    cat(sprintf("Genes present in >= %d cancers: %d\n", 
                min_cancer_types, sum(multi_cancer_genes)))
    
    dmc_matrix_filtered <- dmc_matrix_norm[multi_cancer_genes, , drop = FALSE]
    
    if (nrow(dmc_matrix_filtered) == 0) {
      cat("WARNING: No genes found in >= 2 cancers. Lowering threshold to 1.\n")
      min_cancer_types <- 1
      multi_cancer_genes <- rowSums(dmc_matrix_norm > 0) >= min_cancer_types
      dmc_matrix_filtered <- dmc_matrix_norm[multi_cancer_genes, , drop = FALSE]
    }
  } else {
    stop("dmc_matrix_norm does not exist in environment!")
  }
  
  if (nrow(dmc_matrix_filtered) > top_n_genes) {
    gene_sums <- rowSums(dmc_matrix_filtered)
    top_genes_idx <- order(gene_sums, decreasing=TRUE)[1:top_n_genes]
    dmc_matrix_top <- dmc_matrix_filtered[top_genes_idx, , drop=FALSE]
    cat(sprintf("Selected top %d genes by normalized DMC percentage\n", top_n_genes))
  } else {
    dmc_matrix_top <- dmc_matrix_filtered
    cat(sprintf("Using all %d filtered genes (less than top_n threshold)\n", 
                nrow(dmc_matrix_top)))
  }
  
  dmc_matrix_top <- dmc_matrix_top[rowSums(dmc_matrix_top) > 0, , drop = FALSE]
  dmc_matrix_top <- dmc_matrix_top[, colSums(dmc_matrix_top) > 0, drop = FALSE]
  
  cat(sprintf("After removing zero rows/cols: %d rows, %d cols\n", 
              nrow(dmc_matrix_top), ncol(dmc_matrix_top)))
  
  if (nrow(dmc_matrix_top) == 0 || ncol(dmc_matrix_top) == 0) {
    cat("\nERROR: No data remaining for heatmap after filtering !!!\n")
    cat("Skipping heatmap generation.\n")
  } else {
    cat(sprintf("\n=== GENERATING HEATMAP ===\n"))
    cat(sprintf("Matrix shape: %d genes x %d cancers\n", 
                nrow(dmc_matrix_top), ncol(dmc_matrix_top)))
    cat(sprintf("Number of nonzero cells: %d\n", sum(dmc_matrix_top > 0)))
    cat(sprintf("Value range: %.2f%% to %.2f%%\n", min(dmc_matrix_top), max(dmc_matrix_top)))
    
    plot6_file <- file.path(results_dir, "protein_coding_gene_cancer_heatmap_DMC_NORMALIZED.pdf")
    
    tryCatch({
      pdf(plot6_file, width=14, height=12)
      
      pheatmap(
        dmc_matrix_top,
        color = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize_row = 8,
        fontsize_col = 10,
        main = paste0("Normalized Concordant Gene Body DMCs: Top ", nrow(dmc_matrix_top),
                      " Protein-Coding Genes\n(DMC/CpG × 100%)"),
        border_color = "grey60"
      )
      
      dev.off()
      cat(sprintf("Saved Plot 6 (Heatmap): %s\n", plot6_file))
      
    }, error = function(e) {
      dev.off()
      cat(sprintf("\nERROR generating heatmap: %s\n", e$message))
    })
  }
  
} else {
  cat("No concordant gene body DMCs found.\n")
}

cat("\n=== ALL PLOTS GENERATED (6 TOTAL - ALL NORMALIZED BY CPG COUNT) ===\n")
cat("1. dmc_status_hyper_hypo_protein_coding_per_cancer_NORMALIZED.pdf\n")
cat("2. protein_coding_dmc_cancer_distribution_NORMALIZED.pdf\n")
cat("3. top20_protein_coding_genes_dmcs_NORMALIZED.pdf\n")
cat("4. protein_coding_dmcs_by_chromosome_NORMALIZED.pdf\n")
cat("5. protein_coding_dmcs_per_cancer_NORMALIZED.pdf\n")
cat("6. protein_coding_gene_cancer_heatmap_DMC_NORMALIZED.pdf\n")
cat("\nNormalization formula: (DMC count / CpG count) × 100\n")
cat("Per-cancer normalization: (Sum DMCs / Sum CpGs) × 100\n")
cat("\nDone!\n")