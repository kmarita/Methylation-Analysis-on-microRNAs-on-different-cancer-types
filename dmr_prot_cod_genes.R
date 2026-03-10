library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(GenomicRanges)

main_dir <- "~/methylation_project"
results_dir <- file.path(main_dir, "methylation_analysis_results")

cancer_types <- c("TCGA-BRCA","TCGA-BLCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", 
                  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
                  "TCGA-LIHC","TCGA-LUSC", "TCGA-LUAD", "TCGA-PAAD", "TCGA-PCPG", 
                  "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", 
                  "TCGA-THYM", "TCGA-UCEC")


protein_coding_bed_file <- file.path(main_dir, "prot_cod_genes.bed")

protein_coding_genes <- read.table(protein_coding_bed_file, header = FALSE, sep = "\t", 
                                   stringsAsFactors = FALSE, comment.char = "")
colnames(protein_coding_genes) <- c("Chromosome", "Start", "End", "Gene", "Score", "Strand")

cat(sprintf("Loaded %d protein-coding gene coordinates from BED file\n", nrow(protein_coding_genes)))

genes_gr <- GRanges(
  seqnames = protein_coding_genes$Chromosome,
  ranges = IRanges(start = protein_coding_genes$Start, end = protein_coding_genes$End),
  Gene = protein_coding_genes$Gene
)

read_gene_dmrs_with_counts <- function(cancer_type, results_dir, genes_gr) {
  file_path <- file.path(results_dir, cancer_type, 
                         paste0(cancer_type, "_DMRs_in_protein_coding_genes.tsv"))
  
  if (!file.exists(file_path)) {
    cat(sprintf("DMR file not found for %s\n", cancer_type))
    return(NULL)
  }
  
  dmrs <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  required_cols <- c("Chromosome", "Start", "End")
  missing_cols <- setdiff(required_cols, colnames(dmrs))
  
  if (length(missing_cols) > 0) {
    cat(sprintf("Missing columns in %s: %s\n", cancer_type, paste(missing_cols, collapse=", ")))
    return(NULL)
  }
  
  
  dmrs <- dmrs[!is.na(dmrs$Chromosome) & 
                 !is.na(dmrs$Start) & 
                 !is.na(dmrs$End), ]
  
  dmrs <- dmrs[dmrs$End >= dmrs$Start, ]
  dmrs <- dmrs[dmrs$Start >= 0 & dmrs$End >= 0, ]
  
  if (nrow(dmrs) == 0) {
    cat(sprintf("No valid DMRs for %s after QC\n", cancer_type))
    return(NULL)
  }
  
  if (!"DMR_ID" %in% colnames(dmrs)) {
    dmrs$DMR_ID <- paste0("DMR_", 1:nrow(dmrs))
  }
  
  dmrs$DMR_Length <- dmrs$End - dmrs$Start
  
  dmrs_gr <- GRanges(
    seqnames = dmrs$Chromosome,
    ranges = IRanges(start = dmrs$Start, end = dmrs$End)
  )
  
  overlaps <- findOverlaps(dmrs_gr, genes_gr)
  
  dmrs$Gene <- NA_character_
  if (length(overlaps) > 0) {
    dmrs$Gene[queryHits(overlaps)] <- genes_gr$Gene[subjectHits(overlaps)]
  }
  
  dmrs <- dmrs[!is.na(dmrs$Gene), ]
  
  if (nrow(dmrs) == 0) {
    cat(sprintf("No DMRs overlapping genes for %s\n", cancer_type))
    return(NULL)
  }
  
  result <- list()
  result$all_dmrs <- dmrs
  
  if ("DMR_Status" %in% colnames(dmrs)) {
    result$sig_dmrs <- dmrs[dmrs$DMR_Status %in% c("Hypermethylated", "Hypomethylated"), ]
    cat(sprintf("✓ Loaded %d total DMRs (%d significant) from %s\n", 
                nrow(dmrs), nrow(result$sig_dmrs), cancer_type))
  } else {
    result$sig_dmrs <- dmrs
    cat(sprintf("✓ Loaded %d DMRs from %s (all treated as significant)\n", 
                nrow(dmrs), cancer_type))
  }
  
  return(result)
}


all_dmrs_list <- list()
all_dmrs_total_list <- list()

for (cancer in cancer_types) {
  dmr_data <- read_gene_dmrs_with_counts(cancer, results_dir, genes_gr)
  
  if (!is.null(dmr_data)) {
    all_dmrs_list[[cancer]] <- dmr_data$sig_dmrs
    all_dmrs_total_list[[cancer]] <- dmr_data$all_dmrs
  }
}

if (length(all_dmrs_list) == 0) {
  stop("No protein-coding gene DMR files found with significant results.")
}

cat(sprintf("Total cancer types with significant DMRs: %d\n", length(all_dmrs_list)))
cat(sprintf("Total cancer types with ALL DMRs data: %d\n", length(all_dmrs_total_list)))


dmr_counts_per_gene_cancer <- data.frame()

for (cancer in names(all_dmrs_total_list)) {
  dmr_data <- all_dmrs_total_list[[cancer]]
  
  gene_counts <- dmr_data %>%
    group_by(Gene) %>%
    summarise(Total_DMRs = n(), .groups = "drop") %>%
    mutate(Cancer = cancer)
  
  dmr_counts_per_gene_cancer <- rbind(dmr_counts_per_gene_cancer, gene_counts)
}

cat(sprintf(" Calculated total DMR counts for %d gene-cancer combinations\n", 
            nrow(dmr_counts_per_gene_cancer)))

if (nrow(dmr_counts_per_gene_cancer) == 0) {
  stop("No DMR counts calculated! Check your input files.")
}

dmr_counts_summary <- dmr_counts_per_gene_cancer %>%
  group_by(Gene) %>%
  summarise(
    Mean_DMRs = mean(Total_DMRs),
    Min_DMRs = min(Total_DMRs),
    Max_DMRs = max(Total_DMRs),
    Num_Cancers = n(),
    .groups = "drop"
  )

cat("\nTotal DMR count summary (top 10 genes by mean DMRs):\n")
print(head(dmr_counts_summary[order(-dmr_counts_summary$Mean_DMRs), ], 10))

dmr_counts_file <- file.path(results_dir, "protein_coding_gene_total_DMR_counts_per_cancer.tsv")
write.table(dmr_counts_per_gene_cancer, file = dmr_counts_file, 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n saved total DMR counts: %s\n", dmr_counts_file))

all_dmrs <- data.frame()

for (cancer in names(all_dmrs_list)) {
  dmrs_sig <- all_dmrs_list[[cancer]]
  
  if (nrow(dmrs_sig) > 0) {

        if (!"DMR_Status" %in% colnames(dmrs_sig)) {
      dmrs_sig$DMR_Status <- "Significant"
    }
    
    dmr_info <- dmrs_sig[, c("DMR_ID", "Chromosome", "Start", "End", "Gene", "DMR_Status", "DMR_Length")]
    dmr_info$Cancer <- cancer
    dmr_info$Region <- paste(dmr_info$Chromosome, dmr_info$Start, dmr_info$End, sep = "_")
    all_dmrs <- rbind(all_dmrs, dmr_info)
  }
}

cat(sprintf("\nTotal SIGNIFICANT DMRs across all cancers: %d\n\n", nrow(all_dmrs)))

cat("\n NORMALIZATION \n")
cat("Formula: [Sum(Significant DMRs)] / [Sum(Total DMRs)] per cancer\n\n")

sig_dmr_counts_per_gene_cancer <- all_dmrs %>%
  group_by(Cancer, Gene, DMR_Status) %>%
  summarise(Sig_DMR_Count = n(), .groups = "drop")

total_sig_dmrs_per_gene_cancer <- sig_dmr_counts_per_gene_cancer %>%
  group_by(Cancer, Gene) %>%
  summarise(Total_Sig_DMRs = sum(Sig_DMR_Count), .groups = "drop")

gene_dmr_ratio <- total_sig_dmrs_per_gene_cancer %>%
  left_join(dmr_counts_per_gene_cancer, by = c("Cancer", "Gene"))

gene_dmr_ratio <- gene_dmr_ratio %>%
  filter(!is.na(Total_DMRs))

if (nrow(gene_dmr_ratio) == 0) {
  stop("No matching gene names between significant DMR data and total DMR counts!")
}

cat(sprintf("✓ Successfully matched %d gene-cancer pairs\n", nrow(gene_dmr_ratio)))

normalization_per_cancer <- gene_dmr_ratio %>%
  group_by(Cancer) %>%
  summarise(
    Sum_Sig_DMRs = sum(Total_Sig_DMRs),
    Sum_Total_DMRs = sum(Total_DMRs),
    Normalization_Ratio = Sum_Sig_DMRs / Sum_Total_DMRs,
    .groups = "drop"
  )

sig_dmr_status_with_total <- sig_dmr_counts_per_gene_cancer %>%
  left_join(dmr_counts_per_gene_cancer, by = c("Cancer", "Gene")) %>%
  filter(!is.na(Total_DMRs))

dmr_status_sums <- sig_dmr_status_with_total %>%
  group_by(Cancer, DMR_Status) %>%
  summarise(
    Sum_Sig_DMRs_Status = sum(Sig_DMR_Count),
    Sum_Total_DMRs = sum(Total_DMRs),
    .groups = "drop"
  )

dmr_status_sums <- dmr_status_sums %>%
  mutate(Normalized_Ratio = (Sum_Sig_DMRs_Status / Sum_Total_DMRs) * 100)

dmr_status_sums$Cancer <- factor(dmr_status_sums$Cancer, levels = cancer_types)

normalization_stats_file <- file.path(results_dir, "dmr_normalization_statistics_protein_coding.tsv")
write.table(dmr_status_sums, file = normalization_stats_file, 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n✓ Saved normalization statistics: %s\n", normalization_stats_file))

p_status_norm <- ggplot(dmr_status_sums, aes(x = Cancer, y = Normalized_Ratio, fill = DMR_Status)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sprintf("%.2f%%", Normalized_Ratio)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Hypermethylated" = "steelblue", "Hypomethylated" = "salmon")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Hypermethylated vs Hypomethylated Protein-Coding Gene Body DMRs per Cancer Type",
    subtitle = "Normalized as: (Significant DMR count / Total DMR count) × 100",
    x = "Cancer Type", 
    y = "Normalized DMR Percentage (%)",
    fill = "DMR Status"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "top"
  ) +
  ylim(0, max(dmr_status_sums$Normalized_Ratio) * 1.15)

plot_status_norm_file <- file.path(results_dir, "dmr_status_hyper_hypo_protein_coding_per_cancer_NORMALIZED.pdf")
ggsave(plot_status_norm_file, p_status_norm, width = 14, height = 6)
cat(sprintf("Saved Plot 1: %s\n", plot_status_norm_file))


all_dmrs_info <- data.frame()

for (cancer in names(all_dmrs_list)) {
  dmrs <- all_dmrs_list[[cancer]]
  
  dmrs_df <- data.frame(
    Cancer = cancer,
    Chromosome = dmrs$Chromosome,
    Start = dmrs$Start,
    End = dmrs$End,
    Gene = dmrs$Gene,
    DMR_Status = dmrs$DMR_Status,
    DMR_ID = dmrs$DMR_ID,
    DMR_Length = dmrs$DMR_Length,
    Index = 1:nrow(dmrs),
    stringsAsFactors = FALSE
  )
  
  all_dmrs_info <- rbind(all_dmrs_info, dmrs_df)
}

all_dmrs_gr <- GRanges(
  seqnames = all_dmrs_info$Chromosome,
  ranges = IRanges(start = all_dmrs_info$Start, end = all_dmrs_info$End)
)

cat(sprintf("Total DMRs to analyze: %d\n", length(all_dmrs_gr)))

overlaps_self <- findOverlaps(all_dmrs_gr, all_dmrs_gr, 
                              type = "any", ignore.strand = TRUE)

overlaps_self <- overlaps_self[queryHits(overlaps_self) != subjectHits(overlaps_self)]

cat(sprintf("Found %d pairwise overlaps between DMRs\n", length(overlaps_self)))

dmr_concordance <- data.frame()
processed_dmrs <- rep(FALSE, length(all_dmrs_gr))
cluster_id <- 0

cat("Clustering overlapping DMRs...\n")

for (i in seq_along(all_dmrs_gr)) {
  if (processed_dmrs[i]) next
  
  if (i %% 1000 == 0) {
    cat(sprintf("  Processing DMR %d / %d\n", i, length(all_dmrs_gr)))
  }
  
  cluster_id <- cluster_id + 1
  
  current_cluster <- i
  to_check <- i
  
  while (length(to_check) > 0) {
    current_idx <- to_check[1]
    to_check <- to_check[-1]
    
    new_overlaps <- unique(c(
      subjectHits(overlaps_self[queryHits(overlaps_self) == current_idx]),
      queryHits(overlaps_self[subjectHits(overlaps_self) == current_idx])
    ))
    
    new_ones <- setdiff(new_overlaps, current_cluster)
    if (length(new_ones) > 0) {
      current_cluster <- c(current_cluster, new_ones)
      to_check <- c(to_check, new_ones)
    }
  }
  
  processed_dmrs[current_cluster] <- TRUE
  
  cluster_dmrs <- all_dmrs_info[current_cluster, ]
  
  unique_cancers <- unique(cluster_dmrs$Cancer)
  num_cancers <- length(unique_cancers)
  
  if (num_cancers >= 2) {
    all_hyper <- all(cluster_dmrs$DMR_Status == "Hypermethylated")
    all_hypo <- all(cluster_dmrs$DMR_Status == "Hypomethylated")
    
    if (all_hyper) {
      concordance <- "All_Hypermethylated"
    } else if (all_hypo) {
      concordance <- "All_Hypomethylated"
    } else {
      concordance <- "Discordant"
    }
    
    unique_genes <- unique(cluster_dmrs$Gene)
    unique_genes <- unique_genes[!is.na(unique_genes) & unique_genes != ""]
    
    chr <- cluster_dmrs$Chromosome[1]
    start_coord <- min(cluster_dmrs$Start)
    end_coord <- max(cluster_dmrs$End)
    
    dmr_concordance <- rbind(dmr_concordance, data.frame(
      Cluster_ID = cluster_id,
      Region = paste(chr, start_coord, end_coord, sep = "_"),
      Chromosome = chr,
      Start = start_coord,
      End = end_coord,
      Gene = paste(sort(unique_genes), collapse = ";"),
      Num_Cancers = num_cancers,
      Num_DMRs = nrow(cluster_dmrs),
      Cancers = paste(sort(unique_cancers), collapse = ";"),
      Concordance = concordance,
      stringsAsFactors = FALSE
    ))
  }
}

concordant_dmrs <- dmr_concordance[dmr_concordance$Concordance != "Discordant", ]

cat(sprintf("\nTotal overlapping DMR clusters in 2+ cancers: %d\n", nrow(dmr_concordance)))
cat(sprintf("Concordant DMR clusters (all same direction): %d\n", nrow(concordant_dmrs)))
cat(sprintf("  All Hypermethylated: %d\n", sum(concordant_dmrs$Concordance == "All_Hypermethylated")))
cat(sprintf("  All Hypomethylated: %d\n", sum(concordant_dmrs$Concordance == "All_Hypomethylated")))
cat(sprintf("Discordant DMR clusters (mixed directions): %d\n\n", sum(dmr_concordance$Concordance == "Discordant")))

cat("Distribution of concordant DMRs by number of cancers:\n")
cancer_count_table <- table(dmr_concordance$Num_Cancers)
print(cancer_count_table)
cat("\n")

if (nrow(concordant_dmrs) > 0) {
  
  gene_cancer_dmr <- list()
  
  for (i in 1:nrow(concordant_dmrs)) {
    genes <- unlist(strsplit(concordant_dmrs$Gene[i], ";"))
    genes <- genes[genes != "" & !is.na(genes)]
    genes <- trimws(genes)
    
    cancers <- unlist(strsplit(concordant_dmrs$Cancers[i], ";"))
    cancers <- trimws(cancers)
    
    dmr_region <- concordant_dmrs$Region[i]
    
    for (gene in genes) {
      for (cancer in cancers) {
        key <- paste(gene, cancer, sep = "||")
        if (is.null(gene_cancer_dmr[[key]])) {
          gene_cancer_dmr[[key]] <- character()
        }
        gene_cancer_dmr[[key]] <- c(gene_cancer_dmr[[key]], dmr_region)
      }
    }
  }
  
  all_genes <- unique(unlist(strsplit(concordant_dmrs$Gene, ";")))
  all_genes <- all_genes[all_genes != "" & !is.na(all_genes)]
  all_genes <- trimws(all_genes)
  all_genes <- sort(all_genes)
  
  all_cancers_with_dmrs <- unique(unlist(strsplit(concordant_dmrs$Cancers, ";")))
  all_cancers_with_dmrs <- trimws(all_cancers_with_dmrs)
  all_cancers_with_dmrs <- sort(all_cancers_with_dmrs)
  
  dmr_matrix <- matrix(0, nrow = length(all_genes), ncol = length(all_cancers_with_dmrs),
                       dimnames = list(all_genes, all_cancers_with_dmrs))
  
  for (gene in all_genes) {
    for (cancer in all_cancers_with_dmrs) {
      key <- paste(gene, cancer, sep = "||")
      if (!is.null(gene_cancer_dmr[[key]])) {
        dmr_matrix[gene, cancer] <- length(unique(gene_cancer_dmr[[key]]))
      }
    }
  }
  
  cat("\n MATRIX NORMALIZATION \n")
  cat("Formula: (Significant DMR count / Total DMR count) × 100 for each gene-cancer pair\n\n")
  
  dmr_matrix_norm <- dmr_matrix
  
  for (gene in rownames(dmr_matrix_norm)) {
    for (cancer in colnames(dmr_matrix_norm)) {
      total_dmr_count <- dmr_counts_per_gene_cancer %>%
        filter(Gene == gene, Cancer == cancer) %>%
        pull(Total_DMRs)
      
      if (length(total_dmr_count) > 0 && total_dmr_count > 0) {
        dmr_matrix_norm[gene, cancer] <- (dmr_matrix[gene, cancer] / total_dmr_count) * 100
      } else {
        dmr_matrix_norm[gene, cancer] <- 0
      }
    }
  }
  
  dmr_matrix_df <- as.data.frame(dmr_matrix)
  dmr_matrix_df <- cbind(Gene = rownames(dmr_matrix_df), dmr_matrix_df)
  rownames(dmr_matrix_df) <- NULL
  
  dmr_matrix_norm_df <- as.data.frame(dmr_matrix_norm)
  dmr_matrix_norm_df <- cbind(Gene = rownames(dmr_matrix_norm_df), dmr_matrix_norm_df)
  rownames(dmr_matrix_norm_df) <- NULL
  
  output_matrix_file <- file.path(results_dir, "protein_coding_gene_cancer_dmr_matrix.tsv")
  write.table(dmr_matrix_df, file = output_matrix_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  output_matrix_norm_file <- file.path(results_dir, "protein_coding_gene_cancer_dmr_matrix_NORMALIZED.tsv")
  write.table(dmr_matrix_norm_df, file = output_matrix_norm_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat(sprintf("\nProtein-Coding Gene x Cancer DMR matrix created with %d genes and %d cancer types\n", 
              length(all_genes), length(all_cancers_with_dmrs)))
  cat(sprintf(" matrix saved to: %s\n", output_matrix_file))
  cat(sprintf("Normalized matrix saved to: %s\n\n", output_matrix_norm_file))
  
  cat("Matrix summary:\n")
  cat(sprintf("  Total genes: %d\n", nrow(dmr_matrix)))
  cat(sprintf("  Total cancer types: %d\n", ncol(dmr_matrix)))
  cat(sprintf("  Total non-zero entries: %d\n", sum(dmr_matrix > 0)))
  cat(sprintf("  Max DMRs for any gene-cancer pair: %d\n", max(dmr_matrix)))
  
  gene_totals <- rowSums(dmr_matrix)
  top_genes <- head(sort(gene_totals, decreasing = TRUE), 20)
  cat("\nTop 20 genes by total concordant DMRs across all cancers:\n")
  print(top_genes)
  
  gene_concordant_counts <- data.frame()
  
  for (i in 1:nrow(concordant_dmrs)) {
    genes <- unlist(strsplit(concordant_dmrs$Gene[i], ";"))
    genes <- genes[genes != "" & !is.na(genes)]
    genes <- trimws(genes)
    
    for (gene in genes) {
      concordance_type <- concordant_dmrs$Concordance[i]
      num_cancers <- concordant_dmrs$Num_Cancers[i]
      
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
    
    gene_rows <- grepl(paste0("(^|;)", gene, "(;|$)"), concordant_dmrs$Gene)
    num_cancers <- length(unique(unlist(strsplit(
      paste(concordant_dmrs$Cancers[gene_rows], collapse = ";"), ";"))))
    
    gene_summary <- rbind(gene_summary,
                          data.frame(Gene = gene,
                                     Total_Concordant_DMRs = total_concordant,
                                     All_Hypermethylated_DMRs = all_hyper,
                                     All_Hypomethylated_DMRs = all_hypo,
                                     Num_Cancers = num_cancers,
                                     stringsAsFactors = FALSE))
  }
  
  gene_summary <- gene_summary[order(-gene_summary$Total_Concordant_DMRs), ]
  
  output_file1 <- file.path(results_dir, "concordant_protein_coding_gene_DMRs_summary.tsv")
  write.table(gene_summary, file = output_file1, sep = "\t", quote = FALSE, row.names = FALSE)
  
  output_file2 <- file.path(results_dir, "concordant_protein_coding_DMRs_detailed.tsv")
  write.table(concordant_dmrs, file = output_file2, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("\nAnalysis complete.\n")
  cat(sprintf("Gene-Cancer DMR matrix: %s\n", output_matrix_file))
  cat(sprintf("Gene-Cancer NORMALIZED DMR matrix: %s\n", output_matrix_norm_file))
  cat(sprintf("Concordant gene DMRs summary: %s\n", output_file1))
  cat(sprintf("Concordant DMRs detailed: %s\n", output_file2))
  
  cancer_count_df <- as.data.frame(table(dmr_concordance$Num_Cancers))
  colnames(cancer_count_df) <- c("Num_Cancers", "Count")
  cancer_count_df$Num_Cancers <- as.numeric(as.character(cancer_count_df$Num_Cancers))
  
  total_concordant_dmrs <- sum(cancer_count_df$Count)
  cancer_count_df$Normalized_Percent <- (cancer_count_df$Count / total_concordant_dmrs) * 100
  
  p2 <- ggplot(cancer_count_df, aes(x = Num_Cancers, y = Normalized_Percent)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black") +
    geom_text(aes(label = sprintf("%.1f%%", Normalized_Percent)), 
              vjust = -0.5, size = 3.5) +
    theme_minimal(base_size = 12) +
    labs(title = "Distribution of Concordant Protein-Coding Gene Body DMRs by Number of Cancers",
         x = "Number of Cancers",
         y = "Percentage of Total Concordant DMRs (%)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ylim(0, max(cancer_count_df$Normalized_Percent) * 1.15)
  
  plot2_file <- file.path(results_dir, "protein_coding_dmr_cancer_distribution_NORMALIZED.pdf")
  ggsave(plot2_file, p2, width = 10, height = 6)
  cat(sprintf("Saved Plot 2: %s\n", plot2_file))
  
  top20_summary <- head(gene_summary, 20)
  
  top20_dmr_summary <- data.frame()
  for (gene in top20_summary$Gene) {
    dmr_data <- dmr_counts_per_gene_cancer %>%
      filter(Gene == gene)
    
    if (nrow(dmr_data) > 0) {
      top20_dmr_summary <- rbind(top20_dmr_summary, data.frame(
        Gene = gene,
        Mean_Total_DMRs = mean(dmr_data$Total_DMRs),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  top20_summary <- top20_summary %>%
    left_join(top20_dmr_summary, by = "Gene")
  
  top20_summary$Hyper_Percent <- (top20_summary$All_Hypermethylated_DMRs / top20_summary$Mean_Total_DMRs) * 100
  top20_summary$Hypo_Percent <- (top20_summary$All_Hypomethylated_DMRs / top20_summary$Mean_Total_DMRs) * 100
  
  top20_long <- top20_summary %>%
    select(Gene, Hyper_Percent, Hypo_Percent, Mean_Total_DMRs, Total_Concordant_DMRs) %>%
    pivot_longer(
      cols = c(Hyper_Percent, Hypo_Percent),
      names_to = "DMR_Status",
      values_to = "Percent"
    ) %>%
    mutate(
      DMR_Status = ifelse(DMR_Status == "Hyper_Percent", 
                          "Hypermethylated", 
                          "Hypomethylated")
    )
  
  top20_long$Gene <- factor(top20_long$Gene, levels = rev(top20_summary$Gene))
  top20_long$DMR_Status <- factor(top20_long$DMR_Status, 
                                  levels = c("Hypomethylated", "Hypermethylated"))
  
  p3 <- ggplot(top20_long, aes(x = Gene, y = Percent, fill = DMR_Status)) +
    geom_bar(stat = "identity", color = "black", position = "stack") +
    geom_text(aes(label = ifelse(Percent > 2, sprintf("%.1f%%", Percent), "")),
              position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(values = c("Hypomethylated" = "salmon",
                                 "Hypermethylated" = "steelblue")) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(
      title = "Top 20 Protein-Coding Genes by Concordant DMR Count",
      subtitle = "Normalized as: (Significant DMR count / Mean Total DMR count) × 100",
      x = "Gene",
      y = "Normalized DMR Percentage (%)",
      fill = "DMR Status"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 10)
    )
  
  plot3_file <- file.path(results_dir, "top20_protein_coding_genes_dmrs_NORMALIZED.pdf")
  ggsave(plot3_file, p3, width = 11, height = 8)
  cat(sprintf("Saved Plot 3: %s\n", plot3_file))
  
  # PLOT 4: Chromosome distribution
  chr_counts <- table(concordant_dmrs$Chromosome)
  chr_df <- data.frame(
    Chromosome = names(chr_counts),
    Count = as.numeric(chr_counts)
  )
  
  chr_df <- chr_df[!chr_df$Chromosome %in% c("chrX", "chrY", "chrM"), ]
  
  if (nrow(chr_df) > 0) {
    total_chr_dmrs <- sum(chr_df$Count)
    chr_df$Normalized_Percent <- (chr_df$Count / total_chr_dmrs) * 100
    
    chr_order <- paste0("chr", 1:22)
    chr_df$Chromosome <- factor(chr_df$Chromosome, 
                                levels = chr_order[chr_order %in% chr_df$Chromosome])
    chr_df <- chr_df[order(chr_df$Chromosome), ]
    
    p4 <- ggplot(chr_df, aes(x = Chromosome, y = Normalized_Percent)) +
      geom_bar(stat = "identity", fill = "mediumpurple", color = "black") +
      geom_text(aes(label = sprintf("%.1f%%", Normalized_Percent)), 
                vjust = -0.5, size = 3) +
      theme_minimal(base_size = 11) +
      labs(title = "Concordant Protein-Coding Gene Body DMRs Distribution Across Chromosomes",
           x = "Chromosome",
           y = "Percentage of Total Concordant DMRs (%)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      ylim(0, max(chr_df$Normalized_Percent) * 1.15)
    
    plot4_file <- file.path(results_dir, "protein_coding_dmrs_by_chromosome_NORMALIZED.pdf")
    ggsave(plot4_file, p4, width = 12, height = 6)
    cat(sprintf(" Saved Plot 4: %s\n", plot4_file))
  }
  
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
      title = "Concordant Protein-Coding Gene Body DMRs per Cancer Type",
      subtitle = "Normalized as: (Sum Significant DMRs / Sum Total DMRs) × 100",
      x = "Cancer Type",
      y = "Normalized DMR Percentage (%)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    ) +
    ylim(0, max(cancer_norm_data$Normalized_Percent) * 1.15)
  
  plot5_file <- file.path(results_dir, "protein_coding_dmrs_per_cancer_NORMALIZED.pdf")
  ggsave(plot5_file, p5, width = 12, height = 6)
  cat(sprintf("Saved Plot 5: %s\n", plot5_file))
  
  min_cancer_types <- 2
  top_n_genes <- 20
  
  if (exists("dmr_matrix_norm")) {
    cat(sprintf("Original matrix dimensions: %d genes x %d cancers\n", 
                nrow(dmr_matrix_norm), ncol(dmr_matrix_norm)))
    cat(sprintf("Total non-zero entries: %d\n", sum(dmr_matrix_norm > 0)))
    
    multi_cancer_genes <- rowSums(dmr_matrix_norm > 0) >= min_cancer_types
    cat(sprintf("Genes present in >= %d cancers: %d\n", 
                min_cancer_types, sum(multi_cancer_genes)))
    
    dmr_matrix_filtered <- dmr_matrix_norm[multi_cancer_genes, , drop = FALSE]
    
    if (nrow(dmr_matrix_filtered) == 0) {
      cat("WARNING: No genes found in >= 2 cancers. Lowering threshold to 1.\n")
      min_cancer_types <- 1
      multi_cancer_genes <- rowSums(dmr_matrix_norm > 0) >= min_cancer_types
      dmr_matrix_filtered <- dmr_matrix_norm[multi_cancer_genes, , drop = FALSE]
    }
  } else {
    stop("dmr_matrix_norm does not exist in environment!")
  }
  
  if (nrow(dmr_matrix_filtered) > top_n_genes) {
    gene_sums <- rowSums(dmr_matrix_filtered)
    top_genes_idx <- order(gene_sums, decreasing=TRUE)[1:top_n_genes]
    dmr_matrix_top <- dmr_matrix_filtered[top_genes_idx, , drop=FALSE]
    cat(sprintf("Selected top %d genes by normalized DMR percentage\n", top_n_genes))
  } else {
    dmr_matrix_top <- dmr_matrix_filtered
    cat(sprintf("Using all %d filtered genes (less than top_n threshold)\n", 
                nrow(dmr_matrix_top)))
  }
  
  dmr_matrix_top <- dmr_matrix_top[rowSums(dmr_matrix_top) > 0, , drop = FALSE]
  dmr_matrix_top <- dmr_matrix_top[, colSums(dmr_matrix_top) > 0, drop = FALSE]
  
  cat(sprintf("After removing zero rows/cols: %d rows, %d cols\n", 
              nrow(dmr_matrix_top), ncol(dmr_matrix_top)))
  
  if (nrow(dmr_matrix_top) == 0 || ncol(dmr_matrix_top) == 0) {
    cat("\n!!! ERROR: No data remaining for heatmap after filtering !!!\n")
    cat("Skipping heatmap generation.\n")
  } else {
    cat(sprintf("\n=== GENERATING HEATMAP ===\n"))
    cat(sprintf("Matrix shape: %d genes x %d cancers\n", 
                nrow(dmr_matrix_top), ncol(dmr_matrix_top)))
    cat(sprintf("Number of nonzero cells: %d\n", sum(dmr_matrix_top > 0)))
    cat(sprintf("Value range: %.2f%% to %.2f%%\n", min(dmr_matrix_top[dmr_matrix_top > 0]), max(dmr_matrix_top)))
    
    plot6_file <- file.path(results_dir, "protein_coding_gene_cancer_heatmap_DMR_NORMALIZED.pdf")
    
    tryCatch({
      pdf(plot6_file, width=14, height=12)
      
      pheatmap(
        dmr_matrix_top,
        color = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize_row = 8,
        fontsize_col = 10,
        main = paste0("Normalized Concordant Gene Body DMRs: Top ", nrow(dmr_matrix_top),
                      " Protein-Coding Genes\n(Significant DMR / Total DMR × 100%)"),
        border_color = "grey60"
      )
      
      dev.off()
      cat(sprintf(" Saved Plot 6 (Heatmap): %s\n", plot6_file))
      
    }, error = function(e) {
      dev.off()
      cat(sprintf("\n ERROR generating heatmap: %s\n", e$message))
    })
  }
  
} else {
  cat("No concordant gene body DMRs found.\n")
}

cat("\n=== ALL PLOTS GENERATED (6 TOTAL - ALL NORMALIZED BY TOTAL DMR COUNT) ===\n")
cat("1. dmr_status_hyper_hypo_protein_coding_per_cancer_NORMALIZED.pdf\n")
cat("2. protein_coding_dmr_cancer_distribution_NORMALIZED.pdf\n")
cat("3. top20_protein_coding_genes_dmrs_NORMALIZED.pdf\n")
cat("4. protein_coding_dmrs_by_chromosome_NORMALIZED.pdf\n")
cat("5. protein_coding_dmrs_per_cancer_NORMALIZED.pdf\n")
cat("6. protein_coding_gene_cancer_heatmap_DMR_NORMALIZED.pdf\n")
cat("\nNormalization formula: (Significant DMR count / Total DMR count) × 100\n")
cat("Per-cancer normalization: (Sum Significant DMRs / Sum Total DMRs) × 100\n")
cat("\nDone!\n")