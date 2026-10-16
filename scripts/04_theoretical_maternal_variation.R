# Theoretical maternal variation ------------------------------------------
# Run from the project root (open seed_itv.Rproj first).

dir.create("output/imgs", recursive = TRUE, showWarnings = FALSE)


# Packages --------------------------------------------------------------------
library(ggplot2)
library(patchwork)

# Shared theme and population colours
my_theme <- theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, vjust = 0, size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )
pop_colors <- c("A" = "purple4", "B" = "aquamarine3")

# ---------------------------------------------------------------------------
# Substantial maternal-level variation
ranges <- list(
  A = list(c(1,2.5), c(3,4), c(0,3), c(0.5,3.5)),    # Varied ranges for A
  B = list(c(1.5,2.5), c(3,4.33), c(1.33,4), c(2,5.33))  # Different ranges for B
)
A_vals_maternal <- list(
  A = c(15, 15, 1, 1.8),
  B = c(10, 12, 0.9, 1.5)
)

generate_maternal <- function(pop) {
  do.call(rbind, lapply(1:4, function(i) {
    x_seq <- seq(ranges[[pop]][[i]][1], ranges[[pop]][[i]][2], len = 100)
    data.frame(
      x = x_seq,
      y = A_vals_maternal[[pop]][i] * (x_seq - ranges[[pop]][[i]][1]) *
        (ranges[[pop]][[i]][2] - x_seq),
      Population = paste0(pop),
      Replicate = paste0("Rep", i)
    )
  }))
}
df_right <- rbind(generate_maternal("A"), generate_maternal("B"))

p_right <- ggplot(df_right, aes(x, y, color = Population, group = interaction(Population, Replicate))) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = pop_colors) +
  labs(title = "Substantial Maternal-level Variation", x = "Trait Axis", y = "Frequency") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  my_theme

# ---------------------------------------------------------------------------
# p_top: Population-level curves (for two populations)
x_A_top <- seq(0, 4, length.out = 200)
y_A_top <- x_A_top * (4 - x_A_top)
df_top_A <- data.frame(x = x_A_top, y = y_A_top, Population = "A")
x_B_top <- seq(1.33, 5.33, length.out = 200)
y_B_top <- (x_B_top - 1.33) * (5.33 - x_B_top)
df_top_B <- data.frame(x = x_B_top, y = y_B_top, Population = "B")
df_top <- rbind(df_top_A, df_top_B)

p_top <- ggplot(df_top, aes(x = x, y = y, color = Population)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = pop_colors,
                     guide = guide_legend(override.aes = list(linewidth = 0.8))) +
  labs(title = "Population-level Variation", x = "", y = "Frequency") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  my_theme

# p_left: Without maternal-level variation (using fixed endpoints for each population)
A_vals_left_A <- c(0.8, 1.0, 1.2, 0.9)
df_left_A <- do.call(rbind, lapply(1:4, function(i) {
  x_vals <- seq(0, 4, length.out = 200)
  y_vals <- A_vals_left_A[i] * (x_vals - 0) * (4 - x_vals)
  data.frame(x = x_vals, y = y_vals, Population = "A", Replicate = paste0("Rep", i))
}))
A_vals_left_B <- c(1.0, 0.85, 1.1, 0.95)
df_left_B <- do.call(rbind, lapply(1:4, function(i) {
  x_vals <- seq(1.33, 5.33, length.out = 200)
  y_vals <- A_vals_left_B[i] * (x_vals - 1.33) * (5.33 - x_vals)
  data.frame(x = x_vals, y = y_vals, Population = "B", Replicate = paste0("Rep", i))
}))
df_left <- rbind(df_left_A, df_left_B)
p_left <- ggplot(df_left, aes(x = x, y = y, color = Population, group = interaction(Population, Replicate))) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = pop_colors) +
  labs(title = "Maternal-level Variation Neglected", x = "", y = "Frequency") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  my_theme

# ---------------------------------------------------------------------------
# Minimal maternal-level variation
ranges_centre <- list(
  A = list(c(0,3.7), c(0.1,4), c(0.2,3.9), c(0.3,3.8)),  # A: small differences
  B = list(c(1.33,5.2),   c(1.5,5.33), c(1.45,5.25),   c(1.4,5.15))  # B: small differences
)
A_vals_centre <- list(
  A = c(4, 5, 4.5, 5.5),
  B = c(3.9, 5.1, 4.7, 5.5)
)

generate_centre <- function(pop) {
  do.call(rbind, lapply(1:4, function(i) {
    x_seq <- seq(ranges_centre[[pop]][[i]][1], ranges_centre[[pop]][[i]][2], len = 100)
    data.frame(
      x = x_seq,
      y = A_vals_centre[[pop]][i] * (x_seq - ranges_centre[[pop]][[i]][1]) *
        (ranges_centre[[pop]][[i]][2] - x_seq),
      Population = paste0(pop),
      Replicate = paste0("Rep", i)
    )
  }))
}
df_centre <- rbind(generate_centre("A"), generate_centre("B"))

p_centre <- ggplot(df_centre, aes(x, y, color = Population, group = interaction(Population, Replicate))) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = pop_colors) +
  labs(title = "Minimal Maternal-level Variation", x = "", y = "Frequency") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  my_theme

# ---------------------------------------------------------------------------
# Final Composition:
(final_layout <- p_top / p_left / p_centre / p_right  +
   plot_annotation(tag_levels = 'a') +
   plot_layout(guides ='collect') &
   theme(legend.position='bottom')
)


ggsave("output/imgs/figure_3_theoretical_model.pdf", final_layout, width = 2.8, height = 8.5, units = "in", bg = "white")
