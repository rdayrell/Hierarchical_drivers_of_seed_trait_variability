# Crataegus germination by stratification ---------------------------------
# Run from the project root (open seed_itv.Rproj first).

dir.create("output/imgs", recursive = TRUE, showWarnings = FALSE)
dir.create("output/csv", recursive = TRUE, showWarnings = FALSE)


# Packages --------------------------------------------------------------------
library(brglm2)
library(emmeans)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# LOAD AND RESTRUCTURE DATA
# ------------------------------------------------------------

data <- read.csv("data/crataegus_germination.csv")

# ------------------------------------------------------------

# Check differences among incubation temperatures ----------
# Each embryo-containing seed contributes one record per round.
data_temp <- data %>%
  group_by(SeedID, dish, maxtemp, mintemp) %>%
  summarise(
    germinated = as.integer(any(germination == 1)),
    .groups = "drop"
  ) %>%
  group_by(dish, maxtemp, mintemp) %>%
  summarise(
    total = n(),
    sum_germ = sum(germinated),
    .groups = "drop"
  ) %>%
  mutate(
    germ_prop = sum_germ / total,
    temp_reg = factor(maxtemp)
  )

# FIT MODEL ----------
mod1 <- glm(
  cbind(sum_germ, total - sum_germ) ~ temp_reg,
  family = binomial(link = "logit"),
  data   = data_temp
)

summary(mod1)
anova(mod1)

mod0 <- glm(
  cbind(sum_germ, total - sum_germ) ~ 1,
  family = binomial(link = "logit"),
  data = data_temp
)

anova(mod0, mod1, test = "Chisq")

# ESTIMATES AND CONTRASTS ----------
(emm_pop <- emmeans(mod1, ~ temp_reg, type = "response"))
(pairs_pop <- pairs(emm_pop, adjust = "holm"))


# PLOT FIGURE - Supporting Information ----------

ps <- ggplot(data_temp, aes(x = temp_reg, y = germ_prop)) +
  geom_boxplot(fill = "grey90") +
  scale_y_continuous(
    limits = c(0, 0.7),
    breaks = seq(0, 0.7, 0.1)
  ) +
  scale_x_discrete(
    labels = c(
      "5",
      "10/5",
      "15/5",
      "20/5"
    )
  ) +
  labs(y = "Germination proportion",
       x = "Incubation temperature regime (12/12h, °C)") +
  theme_bw(base_size = 12) +
  theme(panel.border = element_blank(),
        legend.position = "none"
  )

ps

ggsave("output/imgs/figure_s2_crataegus_incubation_temperatures.png", 
       ps, width = 4, height = 3, units = "in", bg = "white")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# RESTRUCTURE DATA ------------------------------------------
# One row per seed with true outcome
data_seed <- data %>%
  group_by(SeedID, dish, tree) %>%
  summarise(
    outcome = case_when(
      any(germination == 1 & strat == "strat1") ~ "strat1",
      any(germination == 1 & strat == "strat2") ~ "strat2",
      TRUE ~ "none"
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# POPULATION-LEVEL
# ------------------------------------------------------------

# COUNTS PER DISH X STRAT ----------
# total is same for both strat rows
data_counts <- data_seed %>%
  group_by(dish) %>%
  summarise(
    total  = n(),
    strat1 = sum(outcome == "strat1"),
    strat2 = sum(outcome == "strat2"),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols      = c(strat1, strat2),
    names_to  = "strat",
    values_to = "germinated"
  ) %>%
  mutate(strat = factor(strat, levels = c("strat1", "strat2")))

# FIT MODEL ----------
mod <- glm(
  cbind(germinated, total - germinated) ~ strat,
  data   = data_counts,
  family = binomial(link = "logit"),
  method = "brglmFit"
)

summary(mod)

# ESTIMATES, CI, AND CONTRASTS ----------
emm_pop <- emmeans(mod, ~ strat, type = "response")
ci_pop  <- confint(emm_pop)
pairs_pop <- pairs(emm_pop, adjust = "holm")

pairs_pop <- as.data.frame(pairs_pop)
ci_pop<- as.data.frame(ci_pop)

pairs_pop$tree <- "population"
ci_pop$tree <- "population"

# PREPARE DATA FOR FIGURE ----------

(pairs_pop$stat_dif <- ifelse(pairs_pop$p.value < 0.05, "*", NA))
(sig_1vs2 <- subset(pairs_pop, pairs_pop$contrast == "strat1 - strat2"))
# No difference

pop_n <- data_counts %>%
  dplyr::filter(strat == "strat1") %>%
  ungroup() %>%
  summarise(n = sum(total))

plot_ci_pop <- ci_pop %>%
  mutate(y_text = max(asymp.UCL) + 0.05,
         n = paste0(pop_n),
         tree = as.factor(1))

# PLOT FIGURE ----------
pd <- position_dodge(width = 0.7)

p_pop <- ggplot(plot_ci_pop, aes(x = tree, y = prob, color = strat, fill = strat, group = strat)) +

  geom_segment(
    x = 1, xend = 1, y = 0, yend = 1,
    colour = "grey92",
    linewidth = 0.5
  ) +

  geom_text(aes(x = tree, y = y_text, label = n), size = 3.5, colour = "black") +
  geom_errorbar(aes(ymin = prob, ymax = asymp.UCL),
                position = pd, width = 0.2, colour = "grey30") +
  geom_bar(stat = "identity", position = pd, colour = "grey30", linewidth = 0.3, width = 0.5) +

  scale_y_continuous(limits = c(0,1.15), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), minor_breaks = seq(0.1, 0.9, by = 0.2)) +
  scale_fill_viridis_d(name = "Stratification", begin = 0.4, end = 0.8, option = "D",
                       labels = c("strat1" = "1",
                                  "strat2" = "2")) +
  labs(y = "Germination probability",
       x = "Population") +
  theme_bw(base_size = 12) +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(colour = "white"),
        legend.position = "none"
  )

p_pop

# ------------------------------------------------------------
# TREE-LEVEL
# ------------------------------------------------------------

# COUNTS PER DISH X STRAT ----------
# total is same for both strat rows

data_counts <- data_seed %>%
  group_by(tree, dish) %>%
  summarise(
    total  = n(),
    strat1 = sum(outcome == "strat1"),
    strat2 = sum(outcome == "strat2"),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols      = c(strat1, strat2),
    names_to  = "strat",
    values_to = "germinated"
  ) %>%
  mutate(tree  = factor(tree),
         strat = factor(strat, levels = c("strat1", "strat2")))

# FIT MODEL ----------

mod <- glm(
  cbind(germinated, total - germinated) ~ strat * tree,
  data   = data_counts,
  family = binomial(link = "logit"),
  method = "brglmFit"
)

# ESTIMATES, CI, AND CONTRASTS ----------

emm_tree <- emmeans(mod, ~ strat | tree, type = "response")
ci_tree  <- confint(emm_tree)
pairs_tree <- pairs(emm_tree, adjust = "holm")

pairs_tree <- as.data.frame(pairs_tree)
ci_tree <- as.data.frame(ci_tree)

# PREPARE DATA FOR FIGURE ----------

(pairs_tree$stat_dif <- ifelse(pairs_tree$p.value < 0.05, "*", NA))
(sig_1vs2 <- subset(pairs_tree, pairs_tree$contrast == "strat1 / strat2"))

sig_1vs2 <- sig_1vs2[, c("tree", "stat_dif")]
sig_1vs2$stat_dif <- ifelse(is.na(sig_1vs2$stat_dif),"", sig_1vs2$stat_dif)

tree_n <- data_counts %>%
  dplyr::filter(strat == "strat1") %>%
  group_by(tree) %>%
  summarise(n = sum(total))
tree_n

ci_tree$tree <- factor(ci_tree$tree, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

plot_ci_tree <- ci_tree %>%
  group_by(tree) %>%
  mutate(
    y_text = max(asymp.UCL) + 0.05) %>%
  right_join(tree_n) %>%
  right_join(sig_1vs2)

# Order tree levels by n
levs <- unique(plot_ci_tree$tree[order(plot_ci_tree$n)])
plot_ci_tree$tree_ord_n <- factor(plot_ci_tree$tree, levels = levs)
plot_ci_tree$n <- ifelse(plot_ci_tree$strat == "strat2", "", plot_ci_tree$n)
plot_ci_tree$stat_dif <- ifelse(plot_ci_tree$strat == "strat2", "", plot_ci_tree$stat_dif)

as.data.frame(plot_ci_tree)
pd <- position_dodge(width = 0.7)

# PLOT FIGURE ----------

x_grid <- data.frame(tree_ord_n = unique(plot_ci_tree$tree_ord_n))

p_tree <- ggplot(plot_ci_tree, aes(x = tree, y = prob, color = strat, fill = strat, group = strat)) +

  geom_segment(
    data = x_grid,
    aes(x = tree_ord_n, xend = tree_ord_n, y = 0, yend = 1),
    inherit.aes = FALSE,
    colour = "grey92",
    linewidth = 0.5
  ) +

  geom_text(aes(x = tree_ord_n, y = y_text, label = n), size = 3.5, colour = "black") +
  geom_text(
    aes(x = tree_ord_n, y = y_text, label = stat_dif),
    fontface = "bold",   # bold entire label
    size = 5,
    colour = "black",
    position = position_nudge(x = +0.3)
  ) +

  geom_errorbar(aes(ymin = prob, ymax = asymp.UCL),
                position = pd, width = 0.2, colour = "grey30") +
  geom_bar(stat = "identity", position = pd, colour = "grey30", linewidth = 0.3, width = 0.7) +

  scale_y_continuous(limits = c(0,1.15), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), minor_breaks = seq(0.1, 0.9, by = 0.2)) +
  scale_fill_viridis_d(name = "Stratification", begin = 0.4, end = 0.8, option = "D",
                       labels = c("strat1" = "1",
                                  "strat2" = "2")) +
  labs(y = "Germination probability",
       x = "Tree") +
  theme_bw(base_size = 12) +
  theme(#panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    legend.box.spacing =  unit(0.1, 'cm'),
    legend.margin = margin(0, 0, 0, 0),
    legend.position = "top",
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.direction = "horizontal"
  )

p_tree


# ------------------------------------------------------------
# POPULATION + TREE-LEVEL FIGURES AND TABLE
# ------------------------------------------------------------

# Extract and combine
tbl_pop <- as.data.frame(pairs_pop) %>%
  rename(group = tree) %>%
  dplyr::select(-null, -df)

tbl_tree <- as.data.frame(pairs_tree) %>%
  rename(group = tree) %>%
  dplyr::select(-null, -df) %>%
  mutate(
    group = paste0("tree ", group)
  )

tbl_combined <- bind_rows(tbl_pop, tbl_tree) %>%
  select(contrast, group, odds.ratio, SE, z.ratio, p.value) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

write.csv(tbl_combined, "output/csv/crataegus_stratification_contrast.csv", row.names = F)

# ------------------------------------------------------------
# SAMPLE SIZE PER TREE vs CI WIDTH
# ------------------------------------------------------------

# Obtain number of seeds per tree
tree_n <- data_counts %>%
  dplyr::filter(strat == "strat1") %>%
  group_by(tree) %>%
  summarise(n_total_tree = sum(total))

# Obtain CI width for each tree x strat
df_ci <- ci_tree %>%
  left_join(tree_n) %>%
  mutate(
    CI_width = asymp.UCL - asymp.LCL
  )

# Plot -----------

plot_ci <- ggplot(df_ci, aes(x = n_total_tree, y = CI_width)) +
  geom_smooth(method = "lm", se = FALSE, aes(colour = strat, linetype = strat), linewidth = 0.7) +
  geom_point(size = 2, aes(colour = strat, fill = strat), stroke = 0.5, alpha = .8) + #
  scale_fill_viridis_d(name = "  ", begin = 0.4, end = 0.8, option = "D",
                       labels = c("strat1" = "1",
                                  "strat2" = "2")) +
  scale_colour_viridis_d(name = "  ", begin = 0.4, end = 0.8, option = "D",
                         labels = c("strat1" = "1",
                                    "strat2" = "2")) +
  scale_linetype_discrete(name = "  ",
                          labels = c("strat1" = "1",
                                     "strat2" = "2")) +
  labs( x = "Seeds per tree (n)",
        y = "95% CI width") +
  theme_bw(base_size = 12) +
  theme(panel.border = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.box.spacing =  unit(0.1, 'cm'),
        legend.margin = margin(0, 0, 0, 0),
        legend.position = "top",
        legend.justification = c(1, 0),
        # legend.title.position = "top",
        legend.background = element_blank(),
        legend.direction = "horizontal"
  )

plot_ci

p_pop_tag   <- p_pop   + labs(tag = "a") +  theme(plot.tag.position = c(0.02, 0.94))
p_tree_tag  <- p_tree  + labs(tag = "b") +  theme(plot.tag.position = c(0.02, 0.94))
plot_ci_tag <- plot_ci + labs(tag = "c") +  theme(plot.tag.position = c(0.02, 0.94))

layout <- c(
  patchwork::area(t = 0, l = 0, b = 5, r = 2),
  patchwork::area(t = 0, l = 3, b = 5, r = 15),
  patchwork::area(t = 0, l = 16, b = 5, r = 21)
)
fig4 = wrap_plots(p_pop_tag, p_tree_tag, plot_ci_tag) + plot_layout(axes = "collect_y") +
    plot_layout(design = layout)

fig4
ggsave("output/imgs/figure_4_crataegus_germination.pdf", fig4, width = 11, height = 3.8, units = "in", bg = "white")
