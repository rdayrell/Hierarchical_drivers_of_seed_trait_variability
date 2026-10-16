# Assemble Figure 2 -------------------------------------------------------
# Run from the project root (open seed_itv.Rproj first).

dir.create("output/imgs", recursive = TRUE, showWarnings = FALSE)

# Packages --------------------------------------------------------------------
library(patchwork)
library(ggplot2)
library(forcats)

# Run scripts/01_seed_quality.R and scripts/02_seed_morphology.R 
# before assembling Figure 2.

p1 <- readRDS("output/rds/plot_map_seed_empty.rds")
p2 <- readRDS("output/rds/plot_fixed_effects_embryoless_seeds.rds")
p3 <- readRDS("output/rds/plot_variance_partition_seed_quality_random_only.rds")
p4 <- readRDS("output/rds/plot_variance_partition_seed_traits.rds")
p4 <- p4 +  theme(panel.spacing = unit(0.5, "cm"))

p1 = p1 +ggtitle("a") +
  theme(plot.title = element_text(hjust = 0, vjust = -9, size = 16),
        plot.margin = margin(t = -20, r = 0, b = -15, l = 0)
        )

p1 <- wrap_elements(p1, clip = FALSE)

p2[[1]] = p2[[1]] + ggtitle("b") +
  theme(plot.title = element_text(hjust = -0.52, vjust = -5, size = 16))

p2 = p2 + theme(
  plot.margin = margin(t = 5, r = 0, b = 20, l = 0)
)

p2 <- wrap_elements(p2, clip = FALSE)

p3 = p3 +ggtitle("c") +
  theme(plot.title = element_text(hjust = -0.9, vjust = -6),

        axis.text.y  = element_text(margin = margin(l = 0)),
        axis.title.y = element_text(margin = margin(l = 0)),

        legend.spacing.y = unit(2, "pt"),
        legend.margin    = margin(t = 0, b = 0),

        plot.margin = margin(t = 0, r = 5, b = 5, l = 0)  # lines
  )

p4 = p4 + ggtitle("d") +
  theme(plot.title = element_text(hjust = -0.3, vjust = -6),

        legend.spacing.y = unit(2, "pt"),
        legend.margin    = margin(t = 0, b = 0),

        plot.margin = margin(t = 0, r = 5, b = 5, l = 0)  # lines)
  )

p34 = p3 + p4 +
    plot_layout(widths = c(5, 9))

p34 <- wrap_elements(p34, clip = FALSE)

layout <- c(
  patchwork::area(t = 0, l = 0,  b = 9,  r = 16),  # p1 (top)
  patchwork::area(t = 10, l = 2,  b = 18,  r = 14),  # p2 (middle, narrower)
  patchwork::area(t = 18, l = 1,  b = 29, r = 16)   # p34 (bottom, full width)
)

fig2 <- p1 + p2 + p34 +
  plot_layout(design = layout)


ggsave("output/imgs/figure_2_variance_decomposition_and_maps.pdf", fig2, width = 10, height = 10, units = "in", bg = "white")
