###############################################################################
# Initialization
###############################################################################

options(repos = "http://cran.us.r-project.org")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, enrichR)

###############################################################################
# Input arguments
###############################################################################

file <- "reports/degs.csv"

pathway <- c("WikiPathway_2021")

###############################################################################
# Import data
###############################################################################

tcc_result <- read_csv(file)

###############################################################################
# Enrichr
###############################################################################

dbs <- listEnrichrDbs() %>%
    as_tibble() %>%
    filter(str_detect(libraryName, pathway))

deg_up <- tcc_result %>%
    filter(estimatedDEG == 1 & str_detect(expression, "high")) %>%
    pull(gene_id)
deg_down <- tcc_result %>%
    filter(estimatedDEG == 1 & str_detect(expression, "low")) %>%
    pull(gene_id)

enrichr_deg_up <- enrichr(deg_up, dbs)
enrichr_deg_down <- enrichr(deg_down, dbs)

df_deg_up <- enrichr_deg_up[[dbs$libraryName]] %>%
    as_tibble() %>%
    filter(Adjusted.P.value < 0.05) %>%
    head(10) %>%
    select(Term, Adjusted.P.value) %>%
    mutate(Term = str_remove(Term, " WP[0-9].*$")) %>%
    mutate(log10AdjPval = -log10(Adjusted.P.value))

df_deg_down <- enrichr_deg_down[[dbs$libraryName]] %>%
    as_tibble() %>%
    filter(Adjusted.P.value < 0.05) %>%
    head(10) %>%
    select(Term, Adjusted.P.value) %>%
    mutate(Term = str_remove(Term, " WP[0-9].*$")) %>%
    mutate(log10AdjPval = -log10(Adjusted.P.value))

###############################################################################
# Bar plot
###############################################################################

g_deg_up <-
    df_deg_up %>%
    mutate(color = "#F06060") %>%
    ggplot(aes(x = fct_reorder(Term, log10AdjPval), y = log10AdjPval, fill = color)) +
    geom_col() +
    scale_fill_identity() +
    scale_y_reverse() +
    coord_flip() +
    labs(x = "-log10 Adjusted P-value", y = "", colour = "") +
    theme_bw() +
    theme(
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "none"
    )

g_deg_down <-
    df_deg_down %>%
    mutate(color = "#5588BB") %>%
    ggplot(aes(x = fct_reorder(Term, log10AdjPval), y = log10AdjPval, fill = color)) +
    geom_col() +
    scale_fill_identity() +
    scale_y_reverse() +
    coord_flip() +
    labs(x = "-log10 Adjusted P-value", y = "", colour = "") +
    theme_bw() +
    theme(
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "none"
    )

###############################################################################
# Export results
###############################################################################

ggsave("reports/enrchr_deg_up.png", g_deg_up, width = 10, height = 6)
ggsave("reports/enrich_deg_up.pdf", g_deg_up, width = 10, height = 6)
ggsave("reports/enrchr_deg_down.png", g_deg_down, width = 10, height = 6)
ggsave("reports/enrich_deg_down.pdf", g_deg_down, width = 10, height = 6)