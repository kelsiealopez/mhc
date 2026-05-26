###############################################################################
# MHC-IIB exon-level expression (TPM) for AW and AC
#
# first read in mhc_counts_matrices.RData (matrix of all counts from all *_q3_primary_counts.tsv)
#   - turns raw counts from STAR whole genome transcript mapping into MHC exon-level TPM per individual and tissue for AW and AC
#   - saves results tables and main plots
#
###############################################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpattern)

proj_dir <- "/n/netscratch/edwards_lab/Lab/kelsielopez/MHC_test_scott"  
setwd(proj_dir)

# read in raw counts of transcripts
counts_rdata <- file.path(proj_dir, "mhc_counts_matrices.RData") 

out_dir <- file.path(proj_dir, "analysis_out_cleaned")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

aw_color <- "#8AA5EE"
ac_color <- "#A9ADB8"

theme_set(
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90")
    )
)

## load the combined counts matrix for all gene expression individuals
load(counts_rdata)

###############################################################################
# AW exon-level TPM per individual
###############################################################################

# counts are length‑normalized by calculating reads per kilobase (RPK),
# then scaled so that the sum of TPMs over all exons in a given individual × tissue is 1,000,000.

aw_annot <- exon_annot[exon_annot$SPECIES == "AW", ]

aw_annot_small <- aw_annot[, c("exon_id",
                               "CONTIG",
                               "START",
                               "STOP",
                               "EXON",
                               "CONTIG_TYPE",
                               "GENE_TYPE",
                               "FRAMESHIFTS_BLAST",
                               "INFRAME_STOPS_BLAST",
                               "FRAMESHIFTS_MINIPROT",
                               "INFRAME_STOPS_MINIPROT",
                               "HAP",
                               "SPECIES")]

aw_annot_small$exon_len_bp <- aw_annot_small$STOP - aw_annot_small$START + 1L
aw_annot_small$exon_len_kb <- aw_annot_small$exon_len_bp / 1000

aw_counts <- counts_long[counts_long$SPECIES == "AW", ]
aw_counts <- aw_counts[, c("exon_id", "count", "sample_id", "tissue")]

aw_counts_df <- as.data.frame(aw_counts)
aw_annot_df  <- as.data.frame(aw_annot_small[, c("exon_id",
                                                 "EXON",
                                                 "CONTIG",
                                                 "CONTIG_TYPE",
                                                 "GENE_TYPE",
                                                 "FRAMESHIFTS_BLAST",
                                                 "INFRAME_STOPS_BLAST",
                                                 "FRAMESHIFTS_MINIPROT",
                                                 "INFRAME_STOPS_MINIPROT",
                                                 "exon_len_bp",
                                                 "exon_len_kb")])

aw_counts_annot <- merge(
  aw_counts_df,
  aw_annot_df,
  by = "exon_id",
  all.x = TRUE
)

aw_counts_annot$tissue <- as.character(aw_counts_annot$tissue)
aw_counts_annot$individual_id <- sub("\\.hap[12]$", "", aw_counts_annot$sample_id)

keep_tissues <- c("Eye", "Heart", "Liver", "Gonad", "Brain", "Blood")
aw_counts_annot <- aw_counts_annot[!is.na(aw_counts_annot$tissue) &
                                     aw_counts_annot$tissue %in% keep_tissues, ]

aw_counts_annot$tissue <- factor(aw_counts_annot$tissue)
aw_counts_annot$EXON   <- factor(aw_counts_annot$EXON)

aw_counts_indiv <- aw_counts_annot %>%
  group_by(
    individual_id, tissue, exon_id, EXON,
    CONTIG_TYPE, GENE_TYPE,
    FRAMESHIFTS_BLAST, INFRAME_STOPS_BLAST,
    FRAMESHIFTS_MINIPROT, INFRAME_STOPS_MINIPROT,
    exon_len_bp, exon_len_kb
  ) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

aw_tpm <- aw_counts_indiv

bad_len_aw <- is.na(aw_tpm$exon_len_kb) | aw_tpm$exon_len_kb <= 0 # sanity check
aw_tpm <- aw_tpm[!bad_len_aw, ]

aw_tpm$RPK <- aw_tpm$count / aw_tpm$exon_len_kb

aw_sf <- aw_tpm %>%
  group_by(individual_id, tissue) %>%
  summarise(scaling_factor = sum(RPK, na.rm = TRUE) / 1e6, .groups = "drop")

aw_tpm <- aw_tpm %>%
  left_join(aw_sf, by = c("individual_id", "tissue"))

aw_tpm$TPM <- ifelse(aw_tpm$scaling_factor > 0,
                     aw_tpm$RPK / aw_tpm$scaling_factor, 0)
aw_tpm$logTPM1 <- log1p(aw_tpm$TPM)

###############################################################################
# AC exon-level TPM per individual
###############################################################################

ac_annot <- exon_annot[exon_annot$SPECIES == "AC", ]

ac_annot_small <- ac_annot[, c("exon_id",
                               "CONTIG",
                               "START",
                               "STOP",
                               "EXON",
                               "CONTIG_TYPE",
                               "GENE_TYPE",
                               "FRAMESHIFTS_BLAST",
                               "INFRAME_STOPS_BLAST",
                               "FRAMESHIFTS_MINIPROT",
                               "INFRAME_STOPS_MINIPROT",
                               "HAP",
                               "SPECIES")]

ac_annot_small$exon_len_bp <- ac_annot_small$STOP - ac_annot_small$START + 1L
ac_annot_small$exon_len_kb <- ac_annot_small$exon_len_bp / 1000

ac_counts <- counts_long[counts_long$SPECIES == "AC", ]
ac_counts <- ac_counts[, c("exon_id", "count", "sample_id", "tissue")]

ac_counts_df <- as.data.frame(ac_counts)
ac_annot_df  <- as.data.frame(ac_annot_small[, c("exon_id",
                                                 "EXON",
                                                 "CONTIG",
                                                 "CONTIG_TYPE",
                                                 "GENE_TYPE",
                                                 "FRAMESHIFTS_BLAST",
                                                 "INFRAME_STOPS_BLAST",
                                                 "FRAMESHIFTS_MINIPROT",
                                                 "INFRAME_STOPS_MINIPROT",
                                                 "exon_len_bp",
                                                 "exon_len_kb")])

ac_counts_annot <- merge(
  ac_counts_df,
  ac_annot_df,
  by = "exon_id",
  all.x = TRUE
)

ac_counts_annot$tissue <- as.character(ac_counts_annot$tissue)
ac_counts_annot$individual_id <- sub("\\.hap[12]$", "", ac_counts_annot$sample_id)

ac_counts_annot <- ac_counts_annot[!is.na(ac_counts_annot$tissue) &
                                     ac_counts_annot$tissue %in% keep_tissues, ]

ac_counts_annot$tissue <- factor(ac_counts_annot$tissue)
ac_counts_annot$EXON   <- factor(ac_counts_annot$EXON)

ac_counts_indiv <- ac_counts_annot %>%
  group_by(
    individual_id, tissue, exon_id, EXON,
    CONTIG_TYPE, GENE_TYPE,
    FRAMESHIFTS_BLAST, INFRAME_STOPS_BLAST,
    FRAMESHIFTS_MINIPROT, INFRAME_STOPS_MINIPROT,
    exon_len_bp, exon_len_kb
  ) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

ac_tpm <- ac_counts_indiv

bad_len_ac <- is.na(ac_tpm$exon_len_kb) | ac_tpm$exon_len_kb <= 0
ac_tpm <- ac_tpm[!bad_len_ac, ]

ac_tpm$RPK <- ac_tpm$count / ac_tpm$exon_len_kb

ac_sf <- ac_tpm %>%
  group_by(individual_id, tissue) %>%
  summarise(scaling_factor = sum(RPK, na.rm = TRUE) / 1e6, .groups = "drop")

ac_tpm <- ac_tpm %>%
  left_join(ac_sf, by = c("individual_id", "tissue"))

ac_tpm$TPM <- ifelse(ac_tpm$scaling_factor > 0,
                     ac_tpm$RPK / ac_tpm$scaling_factor, 0)
ac_tpm$logTPM1 <- log1p(ac_tpm$TPM)

###############################################################################
# combine AW + AC exon-level TPM
###############################################################################

aw_tpm$SPECIES <- "AW"
ac_tpm$SPECIES <- "AC"

tpm_both_exon <- rbind(aw_tpm, ac_tpm)

tpm_both_exon$SPECIES <- factor(tpm_both_exon$SPECIES, levels = c("AW", "AC"))
tpm_both_exon$EXON    <- factor(tpm_both_exon$EXON, levels = c("EXON2", "EXON3"))
tpm_both_exon$tissue  <- droplevels(as.factor(tpm_both_exon$tissue))

###############################################################################
# total TPM vs copy number (core functional EXON2/3)
###############################################################################

aw_core <- aw_tpm[
  aw_tpm$CONTIG_TYPE == "Core" &
    aw_tpm$GENE_TYPE == "Potential_coding_gene" &
    aw_tpm$EXON %in% c("EXON2", "EXON3"),
]

aw_cn <- aw_core %>%
  group_by(individual_id, EXON) %>%
  summarise(copy_number = n_distinct(exon_id), .groups = "drop")

aw_expr <- aw_core %>%
  group_by(individual_id, tissue, EXON) %>%
  summarise(total_TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

aw_expr <- aw_expr %>%
  left_join(aw_cn, by = c("individual_id", "EXON"))

aw_expr$species <- "AW"
aw_expr$copy_number <- as.integer(aw_expr$copy_number)
aw_expr$log2_total_TPM1 <- log2(aw_expr$total_TPM + 1)
aw_expr <- aw_expr[!is.na(aw_expr$copy_number), ]

ac_core <- ac_tpm[
  ac_tpm$CONTIG_TYPE == "Core" &
    ac_tpm$GENE_TYPE == "Potential_coding_gene" &
    ac_tpm$EXON %in% c("EXON2", "EXON3"),
]

ac_cn <- ac_core %>%
  group_by(individual_id, EXON) %>%
  summarise(copy_number = n_distinct(exon_id), .groups = "drop")

ac_expr <- ac_core %>%
  group_by(individual_id, tissue, EXON) %>%
  summarise(total_TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

ac_expr <- ac_expr %>%
  left_join(ac_cn, by = c("individual_id", "EXON"))

ac_expr$species <- "AC"
ac_expr$copy_number <- as.integer(ac_expr$copy_number)
ac_expr$log2_total_TPM1 <- log2(ac_expr$total_TPM + 1)
ac_expr <- ac_expr[!is.na(ac_expr$copy_number), ]

expr_both <- rbind(aw_expr, ac_expr)
expr_both$SPECIES <- factor(expr_both$species, levels = c("AW", "AC"))
expr_both$EXON    <- factor(expr_both$EXON, levels = c("EXON2", "EXON3"))
expr_both$tissue  <- droplevels(as.factor(expr_both$tissue))

###############################################################################
# save files for results
###############################################################################

save(
  aw_tpm, ac_tpm,
  tpm_both_exon,
  aw_expr, ac_expr, expr_both,
  file = file.path(out_dir, "mhc_expression_objects.RData")
)

write.csv(aw_tpm, file.path(out_dir, "tpm_AW_exonlevel.csv"), row.names = FALSE)
write.csv(ac_tpm, file.path(out_dir, "tpm_AC_exonlevel.csv"), row.names = FALSE)

###############################################################################
# PLOTS !
###############################################################################

expr_AW <- aw_expr
expr_AC <- ac_expr

###############################################
# Violin plot TPM by tissue & exon, AW + AC
###############################################

p_violin_exon <- ggplot(
  tpm_both_exon,
  aes(x = tissue, y = TPM, fill = SPECIES)
) +
  geom_violin(
    position = position_dodge(width = 0.8),
    alpha    = 0.5,
    trim     = FALSE
  ) +
  geom_jitter(
    color       = "black",
    position    = position_jitterdodge(jitter.width = 0.15,
                                       dodge.width  = 0.8),
    size        = 0.8,
    alpha       = 0.4,
    show.legend = FALSE
  ) +
  scale_y_continuous(trans = "log1p") +
  scale_fill_manual(values = c(AW = aw_color, AC = ac_color)) +
  facet_wrap(~ EXON, nrow = 1) +
  labs(
    title = "AW + AC: exon-level TPM by tissue and exon type",
    x     = "Tissue",
    y     = "TPM (log1p)",
    fill  = "Species"
  )

print(p_violin_exon)
ggsave(file.path(out_dir, "fig_01_violin_exon_by_tissue_species.png"),
       p_violin_exon, width = 9, height = 5)

##########################################################
# Violin plot EXON2 + EXON3 combined, TPM by tissue
##########################################################

p_violin_combined <- ggplot(
  tpm_both_exon,
  aes(x = tissue, y = TPM, fill = SPECIES)
) +
  geom_violin(
    position = position_dodge(width = 0.8),
    alpha    = 0.5,
    trim     = FALSE
  ) +
  geom_jitter(
    color       = "black",
    position    = position_jitterdodge(jitter.width = 0.15,
                                       dodge.width  = 0.8),
    size        = 0.8,
    alpha       = 0.4,
    show.legend = FALSE
  ) +
  scale_y_continuous(trans = "log1p") +
  scale_fill_manual(values = c(AW = aw_color, AC = ac_color)) +
  labs(
    title = "AW + AC: exon-level TPM by tissue (EXON2 + EXON3 combined)",
    x     = "Tissue",
    y     = "TPM (log1p)",
    fill  = "Species"
  )

print(p_violin_combined)
ggsave(file.path(out_dir, "fig_02_violin_exons_combined_by_tissue.png"),
       p_violin_combined, width = 7, height = 5)

#############################################################
# Boxplot Core vs Not_core by tissue & species
#############################################################

tpm_both_core_tissue <- tpm_both_exon %>%
  filter(CONTIG_TYPE %in% c("Core", "Not_core")) %>%
  mutate(CONTIG_TYPE = factor(CONTIG_TYPE, levels = c("Core", "Not_core")))

p_core_not_by_tissue <- ggplot(
  tpm_both_core_tissue,
  aes(x = CONTIG_TYPE, y = TPM, fill = SPECIES)
) +
  geom_boxplot(
    alpha         = 0.7,
    outlier.alpha = 0.2,
    position      = position_dodge(width = 0.8)
  ) +
  scale_y_continuous(trans = "log1p") +
  scale_fill_manual(values = c(AW = aw_color, AC = ac_color)) +
  facet_wrap(~ tissue, nrow = 1) +
  labs(
    title = "AW + AC: exon-level TPM, Core vs Not_core by tissue",
    x     = "CONTIG_TYPE",
    y     = "TPM (log1p)",
    fill  = "Species"
  )

print(p_core_not_by_tissue)
ggsave(file.path(out_dir, "fig_03_core_vs_not_core_by_tissue.png"),
       p_core_not_by_tissue, width = 10, height = 4)

#################################################################
# Boxplot (pattern): Core vs Not_core by species
#################################################################

tpm_both_core_log2 <- tpm_both_core_tissue %>%
  mutate(log2_TPM1 = log2(TPM + 1))

p_core_not_species_pattern <- ggplot(
  tpm_both_core_log2,
  aes(
    x       = SPECIES,
    y       = log2_TPM1,
    fill    = SPECIES,
    pattern = CONTIG_TYPE,
    group   = interaction(SPECIES, CONTIG_TYPE)
  )
) +
  geom_boxplot_pattern(
    alpha           = 0.8,
    color           = "black",
    position        = position_dodge(width = 0.6),
    outlier.alpha   = 0.3,
    pattern_fill    = "black",
    pattern_colour  = "black",
    pattern_angle   = 45,
    pattern_spacing = 0.03,
    pattern_density = 0.05,
    show.legend     = TRUE
  ) +
  scale_x_discrete(
    labels = c(
      "AW" = "Woodhouse’s\nScrub-Jay",
      "AC" = "Florida\nScrub-Jay"
    )
  ) +
  scale_fill_manual(values = c(AW = aw_color, AC = ac_color), guide = "none") +
  scale_pattern_manual(
    values = c(Core = "none", Not_core = "crosshatch"),
    labels = c(
      Core     = expression(Core~italic("MHC-IIB")),
      Not_core = expression(Outlying~italic("MHC-IIB"))
    ),
    name = NULL
  ) +
  labs(
    y = "Expression level per individual (log2(TPM + 1))"
  ) +
  guides(
    pattern = guide_legend(
      title        = NULL,
      keywidth     = unit(1.2, "lines"),
      keyheight    = unit(0.8, "lines"),
      override.aes = list(
        size         = 0.4,
        linetype     = 1,
        alpha        = 1,
        fill         = "white",
        pattern_fill = "black",
        colour       = "black"
      )
    )
  ) +
  theme(
    legend.position = "top",
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(vjust = 1)
  )

print(p_core_not_species_pattern)
ggsave(file.path(out_dir, "fig_04_core_vs_not_core_pattern_species.png"),
       p_core_not_species_pattern, width = 6, height = 4)

##################################################################
# Boxplot coding vs pseudogene by tissue & species
##################################################################

tpm_both_pseudo_tissue <- tpm_both_exon %>%
  mutate(
    IS_PSEUDOGENE = ifelse(GENE_TYPE == "Potential_pseudogene",
                           "Pseudogene", "Non_pseudogene"),
    IS_PSEUDOGENE = factor(IS_PSEUDOGENE,
                           levels = c("Non_pseudogene", "Pseudogene"))
  )

p_pseudo_by_tissue <- ggplot(
  tpm_both_pseudo_tissue,
  aes(x = IS_PSEUDOGENE, y = TPM, fill = SPECIES)
) +
  geom_boxplot(
    alpha         = 0.7,
    outlier.alpha = 0.2,
    position      = position_dodge(width = 0.8)
  ) +
  scale_y_continuous(trans = "log1p") +
  scale_fill_manual(values = c(AW = aw_color, AC = ac_color)) +
  facet_wrap(~ tissue, nrow = 1) +
  labs(
    title = "AW + AC: exon-level TPM, coding vs pseudogene by tissue",
    x     = "Pseudogene status",
    y     = "TPM (log1p)",
    fill  = "Species"
  )

print(p_pseudo_by_tissue)
ggsave(file.path(out_dir, "fig_05_pseudogene_vs_non_by_tissue.png"),
       p_pseudo_by_tissue, width = 10, height = 4)

#######################################################################
# Boxplot exon-level TPM by tissue & exon, AW + AC (box style)
#######################################################################

p_box_exon_by_tissue <- ggplot(
  tpm_both_exon,
  aes(x = tissue, y = TPM, fill = SPECIES)
) +
  geom_boxplot(
    alpha         = 0.7,
    outlier.alpha = 0.2,
    position      = position_dodge(width = 0.8)
  ) +
  scale_y_continuous(trans = "log1p") +
  scale_fill_manual(values = c(AW = aw_color, AC = ac_color)) +
  facet_wrap(~ EXON, nrow = 1) +
  labs(
    title = "AW + AC: exon-level TPM by tissue and exon type",
    x     = "Tissue",
    y     = "TPM (log1p)",
    fill  = "Species"
  )

print(p_box_exon_by_tissue)
ggsave(file.path(out_dir, "fig_06_box_exon_by_tissue_species.png"),
       p_box_exon_by_tissue, width = 9, height = 4)

###############################################################
# Scatter plot of total TPM vs copy number, AW + AC together
###############################################################

tissue_colors <- c(
  Brain = "#1f78b4",
  Eye   = "#33a02c",
  Gonad = "#6a3d9a",
  Heart = "#e31a1c",
  Liver = "#ff7f00",
  Blood = "gray40"
)

expr_both2 <- expr_both %>%
  mutate(
    log2_total_TPM1 = log2(total_TPM + 1),
    SPECIES         = factor(SPECIES, levels = c("AW", "AC"))
  )

p_scatter_both_exon_tissue_species <- ggplot(
  expr_both2,
  aes(x = copy_number,
      y = log2_total_TPM1,
      color = tissue,
      shape = SPECIES)
) +
  geom_point(
    size     = 2.5,
    alpha    = 0.6,
    position = position_jitter(width = 0.1, height = 0)
  ) +
  geom_smooth(
    method = "lm",
    se     = FALSE,
    aes(group = interaction(tissue, SPECIES))
  ) +
  scale_shape_manual(values = c(AW = 16, AC = 17)) +
  scale_color_manual(values = tissue_colors) +
  facet_wrap(~ EXON, nrow = 1) +
  labs(
    title = "AW + AC: total TPM vs copy number by exon, tissue, and species",
    x     = "Core functional MHC-IIB copies per individual (per exon type)",
    y     = "Total expression per individual (log2(TPM + 1))",
    color = "Tissue",
    shape = "Species"
  )

print(p_scatter_both_exon_tissue_species)
ggsave(file.path(out_dir, "fig_07_scatter_totalTPM_copy_species_tissue.png"),
       p_scatter_both_exon_tissue_species, width = 9, height = 4)

########################################################
# AW only total TPM vs copy number by tissue (Brain, Heart, Liver)
########################################################

aw_expr_sel <- expr_AW %>%
  filter(tissue %in% c("Heart", "Liver", "Brain")) %>%
  mutate(
    tissue          = droplevels(tissue),
    EXON            = factor(EXON, levels = c("EXON2", "EXON3")),
    log2_total_TPM1 = log2(total_TPM + 1)
  )

p_aw_scatter_log2_by_tissue_exon <- ggplot(
  aw_expr_sel,
  aes(x = copy_number,
      y = log2_total_TPM1,
      color = EXON)
) +
  geom_point(
    size     = 2.5,
    alpha    = 0.7,
    position = position_jitter(width = 0.1, height = 0)
  ) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ tissue, nrow = 1) +
  labs(
    title = "AW: total TPM vs copy number, by tissue",
    x     = "Core functional MHC-IIB copies per individual (per exon type)",
    y     = "Total expression per individual (log2(TPM + 1))",
    color = "Exon"
  )

print(p_aw_scatter_log2_by_tissue_exon)
ggsave(file.path(out_dir, "fig_08_AW_scatter_log2_by_tissue_exon.png"),
       p_aw_scatter_log2_by_tissue_exon, width = 9, height = 4)

#############################################################
# AW only total TPM vs copy number (Heart Liver Brain combined)
#############################################################

p_aw_scatter_HLB_log2 <- ggplot(
  aw_expr_sel,
  aes(x = copy_number,
      y = log2_total_TPM1)
) +
  geom_point(
    size     = 3,
    alpha    = 0.7,
    position = position_jitter(width = 0.1, height = 0)
  ) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ EXON, nrow = 1) +
  labs(
    title = "AW: total TPM vs copy number (Heart, Liver, Brain)",
    x     = "Core functional MHC-IIB copies per individual (per exon type)",
    y     = "Total expression level per individual (log2(TPM + 1))"
  )

print(p_aw_scatter_HLB_log2)
ggsave(file.path(out_dir, "fig_09_AW_scatter_HLB_log2_exon.png"),
       p_aw_scatter_HLB_log2, width = 7, height = 4)

###############################################################################
# final numbers for results (mean, standard deviation (SD), standard error (SE) of TPM)
###############################################################################

summary_by_species_tissue_exon <- tpm_both_exon %>%
  group_by(SPECIES, tissue, EXON) %>%
  summarise(
    n        = n(),
    mean_TPM = mean(TPM, na.rm = TRUE),
    sd_TPM   = sd(TPM,   na.rm = TRUE),
    se_TPM   = sd_TPM / sqrt(n),
    .groups  = "drop"
  )

print(summary_by_species_tissue_exon)

write.csv(
  summary_by_species_tissue_exon,
  file = file.path(out_dir, "TPM_summary_bySpecies_tissue_exon_forScott.csv"),
  row.names = FALSE
)

summary_by_species_tissue_combined <- tpm_both_exon %>%
  group_by(SPECIES, tissue) %>%
  summarise(
    n        = n(),
    mean_TPM = mean(TPM, na.rm = TRUE),
    sd_TPM   = sd(TPM,   na.rm = TRUE),
    se_TPM   = sd_TPM / sqrt(n),
    .groups  = "drop"
  )

print(summary_by_species_tissue_combined)

write.csv(
  summary_by_species_tissue_combined,
  file = file.path(out_dir, "TPM_summary_bySpecies_tissue_exonsCombined_forScott.csv"),
  row.names = FALSE
)

summary_by_tissue <- tpm_both_exon %>%
  group_by(tissue) %>%
  summarise(
    n        = n(),
    mean_TPM = mean(TPM, na.rm = TRUE),
    sd_TPM   = sd(TPM,   na.rm = TRUE),
    se_TPM   = sd_TPM / sqrt(n),
    .groups  = "drop"
  )

print(summary_by_tissue)

summary_by_species <- tpm_both_exon %>%
  group_by(SPECIES) %>%
  summarise(
    n        = n(),
    mean_TPM = mean(TPM, na.rm = TRUE),
    sd_TPM   = sd(TPM,   na.rm = TRUE),
    se_TPM   = sd_TPM / sqrt(n),
    .groups  = "drop"
  )

print(summary_by_species)

summary_by_species_exon <- tpm_both_exon %>%
  group_by(SPECIES, EXON) %>%
  summarise(
    n        = n(),
    mean_TPM = mean(TPM, na.rm = TRUE),
    sd_TPM   = sd(TPM,   na.rm = TRUE),
    se_TPM   = sd_TPM / sqrt(n),
    .groups  = "drop"
  )

print(summary_by_species_exon)
