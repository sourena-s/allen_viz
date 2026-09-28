#module load 2025;module load R/4.5.1-gfbf-2025a


library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)
library(fastICA)

donor_ids <- c(178236545, 178238266, 178238316, 178238359, 178238373, 178238387)

# Probe annotation (identical across donors)
probe_annotation <- read.csv(
  file.path("raw", donor_ids[1], "Probes.csv"),
  header = TRUE,
  quote = "\"",
  stringsAsFactors = FALSE
)

# Helper: load one donor's expression matrix + sample annotation
# and return a probes x samples expression matrix with proper row/col names
load_donor_data <- function(donor_id) {

  # MicroarrayExpression.csv has NO header row.
  # Column 1 = probe_id, remaining columns = one per sample (well_id order
  # matches the row order of SampleAnnot.csv for that donor)
  donor_expression <- read.csv(
    file.path("raw", donor_id, "MicroarrayExpression.csv"),
    header = FALSE,
    quote = "\"",
    stringsAsFactors = FALSE
  )

  probe_ids <- donor_expression[[1]]
  donor_expression <- as.matrix(donor_expression[, -1, drop = FALSE])
  rownames(donor_expression) <- probe_ids

  # --- Sample (brain region) annotation ---
  sample_annotation <- read.csv(
    file.path("raw", donor_id, "SampleAnnot.csv"),
    header = TRUE,
    quote = "\"",
    stringsAsFactors = FALSE
  )

  # Sanity check: number of samples must match number of expression columns
  if (nrow(sample_annotation) != ncol(donor_expression)) {
    stop(sprintf(
      "Donor %s: mismatch between SampleAnnot rows (%d) and expression columns (%d)",
      donor_id, nrow(sample_annotation), ncol(donor_expression)
    ))
  }

  # Sample names = slab type, structure acronym, donor id and full structure name
  colnames(donor_expression) <- paste(sample_annotation$slab_type, sample_annotation$structure_acronym, donor_id, "____", sample_annotation$structure_name, sep = "_")

  list(expression = donor_expression, probe_ids = probe_ids)
}

# -----------------------------------------------------------------
# Load all donors
# -----------------------------------------------------------------
donor_data <- lapply(donor_ids, load_donor_data)
names(donor_data) <- donor_ids

# -----------------------------------------------------------------
# Sanity check: all donors must share the same probe_id row order
# (they normally do, since Probes.csv is identical across donors)
# -----------------------------------------------------------------
reference_probe_ids <- donor_data[[1]]$probe_ids
for (donor in names(donor_data)) {
  if (!identical(donor_data[[donor]]$probe_ids, reference_probe_ids)) {
    stop(sprintf("Donor %s has a different probe_id order than the reference donor", donor))
  }
}

# Concatenate all donors column-wise (cbind), since rows (probes) are aligned across donors
all_donors_expression <- do.call(cbind, lapply(donor_data, function(x) x$expression))


rownames(all_donors_expression) <- reference_probe_ids

dim(all_donors_expression)
head(colnames(all_donors_expression))

probe_gene_info <- read.table(file.path("raw", "probes_updated_geneids_entrez.csv"), header = TRUE, sep = ",")
gene_length <- probe_gene_info[, 8]
gene_symbol <- probe_gene_info[, 7]
probe_name <- probe_gene_info[, 3]

#gene_length_threshold <- 250000
gene_length_threshold <- 500000
#gene_length_threshold <- 1000000

is_long_gene <- gene_length > gene_length_threshold & !is.na(gene_length)

# Probes x samples expression for long genes only
long_gene_expression <- all_donors_expression[is_long_gene, ]
rownames(long_gene_expression) <- paste(gene_symbol[is_long_gene], probe_name[is_long_gene], sep = " ")
stopifnot(!any(duplicated(rownames(long_gene_expression))))

colnames(long_gene_expression) <- make.unique(colnames(long_gene_expression))

dim(long_gene_expression)

ica_result <- fastICA(long_gene_expression, n.comp = 20, alg.typ = "parallel", fun = "logcosh", method = "C", row.norm = FALSE)

# Mixing matrix: independent components x samples
ic_sample_loadings <- ica_result$A
rownames(ic_sample_loadings) <- paste0("IC", seq_len(nrow(ic_sample_loadings)))
colnames(ic_sample_loadings) <- colnames(long_gene_expression)

#ic_distance <- dist(long_gene_expression)
#sample_distance <- dist(t(long_gene_expression))

ic_distance <- dist(ic_sample_loadings)
sample_distance <- dist(t(ic_sample_loadings))

ic_clustering <- hclust(ic_distance, method = "average")
sample_clustering <- hclust(sample_distance, method = "average")

gene_probe_order <- rownames(long_gene_expression)[ic_clustering$order]
sample_order <- colnames(long_gene_expression)[sample_clustering$order]
#sample_order <- sort(colnames(long_gene_expression))

expression_ordered <- long_gene_expression[gene_probe_order, sample_order]
ic_loadings_ordered <- ic_sample_loadings[, sample_order]

# Melt to long format for ggplot2, preserving cluster order via factors
# -----------------------------------------------------------------
expression_long <- melt(expression_ordered, varnames = c("gene_probe", "sample"), value.name = "expression")
expression_long$gene_probe <- factor(expression_long$gene_probe, levels = gene_probe_order)
expression_long$sample     <- factor(expression_long$sample, levels = sample_order)


ic_loadings_long <- melt(ic_loadings_ordered, varnames = c("IC", "sample"), value.name = "IC_loading")
ic_loadings_long$sample <- factor(ic_loadings_long$sample, levels = sample_order)



# geom_jitter(width = 0.15, size = 0.4, alpha = 1, color = "steelblue") +

expression_boxplot <- ggplot(expression_long, aes(x = sample, y = expression)) +
  stat_summary(geom="boxplot",fun.data=function(x)setNames(quantile(x,c(.05,.25,.5,.75,.95),na.rm=TRUE),c("ymin","lower","middle","upper","ymax")),fill="grey90",color="grey30",width=.6) +
  stat_summary(fun=median,geom="errorbar",aes(ymin=after_stat(y),ymax=after_stat(y)),width=1,color="red") +
  geom_hline(yintercept=8, color="blue", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=9, color="red", linetype="dashed", linewidth=0.5)+
  theme_minimal(base_size = 8) +
  theme(legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(y = "Expression")


#  scale_fill_gradientn(colours = c("white", "yellow", "orange", "red", "darkred"), values = scales::rescale(c(0, 1, 3, 5, 15)), limits = c(0, 15), oob = scales::squish) +

ic_loading_heatmap <- ggplot(ic_loadings_long, aes(x = sample, y = IC, fill = IC_loading)) +
  geom_tile() + scale_fill_viridis_c(option = "viridis") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "none",
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4),
    axis.text.y = element_text(size = 4),
    panel.grid = element_blank()
  ) +
  labs(x = "Sample (region)", y = "Independent component", fill = "IC loading")

combined_plot <- ic_loading_heatmap / expression_boxplot + plot_layout(heights = c(1, 1))

ggsave("heatmap.png", plot = combined_plot, width = 500, height = 200, units = "cm", limitsize = FALSE, dpi=200)
