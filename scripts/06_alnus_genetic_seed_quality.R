# Alnus genetic diversity and seed quality --------------------------------
# Run from the project root (open seed_itv.Rproj first).

dir.create("output/imgs", recursive = TRUE, showWarnings = FALSE)
dir.create("output/csv", recursive = TRUE, showWarnings = FALSE)


# Packages --------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(patchwork)

# Import and prepare Alnus glutinosa data ------------------------------------

gen_alnus = read.csv("data/genetics_alnus.csv")

embryoless_df <- read.csv("data/embryoless_seeds.csv")

embryoless_df = embryoless_df %>%
  select(-starts_with("CHELSA_")) %>%
  dplyr::filter(Species == "Alnus glutinosa")

sum_emb = embryoless_df %>%
  ungroup() %>%
  mutate(empty_per = Empty/Sown) %>%
  group_by(Species, SBD) %>%
  summarise(n_seedtraits = length(unique(test_tree_number)),
            empty_med = round(median(empty_per),2),
            empty_sd = sd(empty_per)
  ) %>%
  dplyr::filter(n_seedtraits > 7L) %>%
  mutate(SBD = as.integer(as.character(SBD)))

data = sum_emb %>%
  inner_join(., gen_alnus, by = c("Species","SBD"))

names(data)

##--------------------------------------------------------
# Correlations and plots ------------------------------------------------------

seed_vars <- c("empty_med", "empty_sd")
genetic_vars <- c("empty_med","P", "Obs_Het", "Pi", "Fis")

model_results <- list()
plots <- list()
results_table <- data.frame()  # Table to store model results


for (genetic in genetic_vars) {
  for (seed in seed_vars) {

    x <- data[[genetic]]
    y <- data[[seed]]

    mod <- cor.test(x, y, method = 'pearson')
    model_results[[paste(genetic, seed, sep = "_vs_")]] <- list(
      genetic = genetic,
      seed = seed,
      estimate = mod$estimate,
      p.value = mod$p.value
    )

    estimate <- round(mod$estimate,2)
    p_value <- round(mod$p.value,3)
    ci_low <- round(mod$conf.int[1],2)
    ci_high <- round(mod$conf.int[2],2)
    method <- mod$method

    results_table <- rbind(results_table, data.frame(
      Seed_var = seed,
      Genetic_var = genetic,
      Estimate = estimate,
      P_Value = p_value,
      ci_low = ci_low,
      ci_high = ci_high,
      Method = method
    ))

    y_label <- if (seed == "empty_med") "Median embryoless seed fraction" else "SD embryoless seed fraction"
    y_limits <- if (seed == "empty_med") c(-0.02, 1.02) else c(-0.02, 0.5)

    x_label <- if (genetic == "Obs_Het") {
      "Observed heterozygosity"
    } else if (genetic == "Pi") {
      "Nucleotide diversity"
    } else if (genetic == "P") {
      "Mean frequency of the most\ncommon allele at each locus "
    } else if (genetic == "Fis") {
      "Inbreeding coefficient"
    } else if (genetic == "empty_med") {
      "Median empty seed fraction"
    } else {
      "ERROR!"
    }

    # Label every correlation; line style indicates significance
    plot_label <- sprintf("atop(r == %.2f, italic(P) == %.3f)", estimate, p_value)

    p <- ggplot(data = data) +
      geom_point(aes(x = .data[[genetic]], y = .data[[seed]]), size = 2, stroke = 0.8, alpha = 0.5) +
      scale_y_continuous(name = y_label, limits = y_limits) +
      scale_x_continuous(name = x_label) +
      annotate("text", x = max(data[[genetic]], na.rm = TRUE),
               y = min(y_limits) + 0.92 * (max(y_limits) - min(y_limits)),
               label = plot_label,
               hjust = 1,
               size = 4,
               parse = TRUE) +
      theme_minimal()

    if (p_value < 0.05) {
      p <- p +
        geom_smooth(aes(x = .data[[genetic]], y = .data[[seed]]),
                    method = "lm", color = "black", fill = "black", alpha = 0) +
        theme_minimal()

    } else {
      p <- p +
        geom_smooth(aes(x = .data[[genetic]], y = .data[[seed]]),
                    method = "lm", color = "black", fill = "black", alpha = 0, linetype = "dashed", linewidth = 0.5) +
        theme_minimal()
    }
    plots[[paste(seed, genetic, sep = "_")]] <- p

  }
}

##--------------------------------------------------------
combined_plotA <- wrap_plots(plots[-c(1,2)], ncol = 2)

combined_plotA = combined_plotA +
  plot_layout(axes = "collect_y")  +
  plot_annotation(tag_levels = 'a') &
  theme(
    axis.title.y = element_text(margin = margin(r = 20), size = 12),
    axis.text = element_text(size = 9),
    plot.tag.position = c(0.1, 0.95),
    plot.tag = element_text(face = "bold"),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
  )

ggsave("output/imgs/figure_5_empty_genetic_diversity_alnus.pdf", plot = combined_plotA, width = 7, height = 9.2, dpi = 300)
write.csv(results_table, "output/csv/empty_genetic_diversity_table_alnus.csv", row.names = FALSE)


ps <- plots[["empty_sd_empty_med"]]

ggsave("output/imgs/figure_s3_alnus_embryoless_median_sd.png", plot = ps, width = 3, height = 3, dpi = 300, bg = "white")
