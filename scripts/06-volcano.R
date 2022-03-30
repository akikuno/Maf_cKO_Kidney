###############################################################################
# Initialization
###############################################################################

options(repos = "http://cran.us.r-project.org")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, ggrepel)

###############################################################################
# Input arguments
###############################################################################

file <- "reports/degs.csv"

###############################################################################
# Import data
###############################################################################

tcc_result <- read_csv(file)

###############################################################################
# Volucano and MA plots
###############################################################################

TARGETS <- c("SLC5A2", "SLC2A2", "NOX4", "MAF")
TARGETS <- str_c(TARGETS, "$")
TARGETS <- str_c(TARGETS, collapse = "|")
# qval_min10 <- sort(tcc_result$q.value)[10]
# mval_min10 <- sort(tcc_result$m.value)[10]
# mval_max10 <- sort(tcc_result$m.value, decreasing = TRUE)[10]

tcc_volcano <-
    tcc_result %>%
    mutate(log10qval = -log10(q.value)) %>%
    mutate(label = case_when(
        estimatedDEG == 1 & expression == "cko-high" ~ "cKO-UP",
        estimatedDEG == 1 & expression == "cko-low" ~ "cKO-DOWN",
        TRUE ~ "NO"
    )) %>%
    mutate(annotate = case_when(
        estimatedDEG == 1 & str_detect(gene_id, TARGETS) ~ gene_id,
        TRUE ~ as.character(NA)
    ))


mycolors <- c("#5588BB", "#F06060", "#E5E5E5")
names(mycolors) <- c("cKO-DOWN", "cKO-UP", "NO")

g_volcano <-
    ggplot(tcc_volcano, aes(x = m.value, y = log10qval, color = label, label = annotate)) +
    geom_point() +
    scale_colour_manual(values = mycolors) +
    geom_label_repel(box.padding = 1, max.overlaps = Inf, show.legend = FALSE) +
    labs(x = "log2 Fold Change", y = "-log10 q-value", colour = "") +
    theme_bw() +
    theme(
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.background = element_rect(fill = "white")
    )

###############################################################################
# Export results
###############################################################################

ggsave("reports/volcano.png", g_volcano, width = 10, height = 5)
ggsave("reports/volcano.pdf", g_volcano, width = 10, height = 5)