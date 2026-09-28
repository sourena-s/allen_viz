#module load 2025;module load R/4.5.1-gfbf-2025a


library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)
library(fastICA)

donor_ids <- c(178236545, 178238266, 178238316, 178238359, 178238373, 178238387)

# Probe annotation (identical across donors)
probe_annot <- read.csv(
  file.path("raw", donor_ids[1], "Probes.csv"),
  header = TRUE,
  quote = "\"",
  stringsAsFactors = FALSE
)

# Helper: load one donor's expression matrix + region annotation
# and return an expression matrix with proper row/col names
load_donor_data <- function(donor_id) {

  # MicroarrayExpression.csv has NO header row.
  # Column 1 = probe_id, remaining columns = one per sample (well_id order
  # matches the row order of SampleAnnot.csv for that donor)
  exp_mat <- read.csv(
    file.path("raw",donor_id, "MicroarrayExpression.csv"),
    header = FALSE,
    quote = "\"",
    stringsAsFactors = FALSE
  )

  probe_ids <- exp_mat[[1]]
  exp_mat <- as.matrix(exp_mat[, -1, drop = FALSE])
  rownames(exp_mat) <- probe_ids

  # --- Region / sample annotation ---
  region_annot <- read.csv(
    file.path("raw",donor_id, "SampleAnnot.csv"),
    header = TRUE,
    quote = "\"",
    stringsAsFactors = FALSE
  )

  # Sanity check: number of samples must match number of expression columns
  if (nrow(region_annot) != ncol(exp_mat)) {
    stop(sprintf(
      "Donor %s: mismatch between SampleAnnot rows (%d) and expression columns (%d)",
      donor_id, nrow(region_annot), ncol(exp_mat)
    ))
  }

  # Column names = structure_name augmented with donor_id
colnames(exp_mat) <- paste(region_annot$slab_type,region_annot$structure_acronym,donor_id, "____", region_annot$structure_name, sep = "_")

  list(exp_mat = exp_mat, probe_ids = probe_ids)
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
ref_probe_ids <- donor_data[[1]]$probe_ids
for (d in names(donor_data)) {
  if (!identical(donor_data[[d]]$probe_ids, ref_probe_ids)) {
    stop(sprintf("Donor %s has a different probe_id order than the reference donor", d))
  }
}

# Concatenate all donors column-wise (cbind), since rows (probes) are aligned across donors
concatenated_matrix <- do.call(cbind, lapply(donor_data, function(x) x$exp_mat))


rownames(concatenated_matrix) <- ref_probe_ids

dim(concatenated_matrix)
head(colnames(concatenated_matrix))

probe_extra <- read.table(file.path("raw", "probes_updated_geneids_entrez.csv"), header = TRUE, sep = ",")
gene_length <- probe_extra[, 8]
gene_symbol <- probe_extra[, 7]
probe_name <- probe_extra[, 3]

#gene_length_threshold <- 250000
gene_length_threshold <- 500000
#gene_length_threshold <- 1000000

gene_symbol_filtered <- gene_symbol[gene_length > gene_length_threshold & !is.na(gene_length) ]
probe_name_filtered <- probe_name[gene_length > gene_length_threshold & !is.na(gene_length)  ]

mat <- concatenated_matrix[gene_length > gene_length_threshold & !is.na(gene_length)  , ]
rownames(mat) <- paste(gene_symbol_filtered, probe_name_filtered, sep = " ")
stopifnot(!any(duplicated(rownames(mat))))

colnames(mat) <- make.unique(colnames(mat))

dim(mat)

ica <- fastICA(mat, n.comp = 20, alg.typ = "parallel", fun = "logcosh", method = "C", row.norm=F, )

rnames<-rownames(mat)
cnames<-colnames(mat)



mat2<-ica$A

rownames(mat) <-rnames
rownames(mat2) <- paste0("IC", seq_len(nrow(mat2)))

colnames(mat) <-cnames
colnames(mat2) <-cnames

#row_dist <- dist(mat)
#col_dist <- dist(t(mat))

row_dist <- dist(mat2)
col_dist <- dist(t(mat2))

row_hclust <- hclust(row_dist, method = "average")
col_hclust <- hclust(col_dist, method = "average")

row_order <- rownames(mat)[row_hclust$order]
col_order <- colnames(mat)[col_hclust$order]
#col_order <- sort(colnames(mat))

mat_ordered <- mat[row_order, col_order]
mat2_ordered <- mat2[,col_order]

# Melt to long format for ggplot2, preserving cluster order via factors
# -----------------------------------------------------------------
mat_long <- melt(mat_ordered, varnames = c("gene_probe", "region"), value.name = "expression")
mat_long$gene_probe <- factor(mat_long$gene_probe, levels = row_order)
mat_long$region     <- factor(mat_long$region, levels = col_order)


mat2_long <- melt(mat2_ordered, varnames = c("IC", "region"), value.name = "IC_loading")
mat2_long$region     <- factor(mat2_long$region, levels = col_order)



# geom_jitter(width = 0.15, size = 0.4, alpha = 1, color = "steelblue") +

p_bottom <- ggplot(mat_long, aes(x = region, y = expression)) +
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

p_top <- ggplot(mat2_long, aes(x = region, y = IC, fill = IC_loading)) +
  geom_tile() + scale_fill_viridis_c(option = "viridis") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "none",
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4),
    axis.text.y = element_text(size = 4),
    panel.grid = element_blank()
  ) +
  labs(x = "Region", y = "Gene_Probe", fill = "Expression")

combined_plot <- p_top / p_bottom + plot_layout(heights = c(1, 1))

ggsave("heatmap.png", plot = combined_plot, width = 500, height = 200, units = "cm", limitsize = FALSE, dpi=200)
