# Seed morphology and hierarchical variance partitioning ---------------------------------
# Run from the project root (open seed_itv.Rproj first).

dir.create("output/rds", recursive = TRUE, showWarnings = FALSE)
dir.create("output/csv", recursive = TRUE, showWarnings = FALSE)

# Packages --------------------------------------------------------------------
library(dplyr)
library(lme4)
library(ggplot2)
library(DHARMa)
library(corrr)
library(viridis)
library(tidyverse)
library(lmerTest)

# Import and prepare data -----------------------------------------------------
embryoless_data <- read.csv("data/embryoless_seeds.csv")
sum_emb <- embryoless_data %>%
  select(-test_tree_number, -Empty, -Sown) %>%
  distinct()
  
names(embryoless_data)

morph_data = read.csv("data/morphology_tpif.csv")
names(morph_data)

# merge datasets
dat = inner_join(sum_emb, morph_data, by = c("Species", "SBD"))

dat = dat %>%
  dplyr::filter(Species == "Betula pubescens" |
                Species == "Pinus sylvestris") %>%
  mutate(uniTree = paste0(SBD,"_", Tree))

dat$AR = 1/dat$AR

(n_seeds = dat %>%
    ungroup() %>%
    group_by(Species, SBD, uniTree) %>%
    summarise(n_seeds = n())
  )

n_seeds %>%
  ungroup() %>%
  summarise(tot_seeds = sum(n_seeds))


# Sample sizes for Supporting Information Notes S1
(n_trees = n_seeds %>%
    ungroup() %>%
    group_by(Species,SBD) %>%
    summarise(n_trees = length(unique(uniTree)),
    ) %>%
    group_by(Species) %>%
    summarise(n_populations = length(unique(SBD)),
              n_all_trees = sum(n_trees),
              min_trees = min(n_trees),
              max_trees = max(n_trees)
    ))

(n_seeds_sum = n_seeds %>%
  ungroup() %>%
  group_by(Species) %>%
  summarise(median_seeds = median(n_seeds),
            min_seeds = min(n_seeds),
            max_seeds = max(n_seeds)
  ))

# colour -----------------------------------------------------------------------

# # convert sRGB to Lab

Lab <- convertColor(dat[,c("Red", "Green", "Blue")]/255, from="sRGB", to="Lab")
dat = cbind(dat, Lab)

# Correlation -------------------------------------------------------------------

names(dat)

var_morph <- c("Feret", "MinFeret", "AR", "Circ.", "L", "a", "b", "Area", "Perim.", "Round", "Solidity" )

# dat %>%
#   dplyr::filter(Species == "Betula pubescens") %>%
#   dplyr::select_at(var_morph) %>%
#   # removed (highly correlation > 0.8): Perim., Area, Round, Solidity
#   correlate() %>%
#   # Re-arrange a correlation data frame to group highly correlated variables closer together
#   rearrange(method = "MDS", absolute = FALSE) %>%
#   shave() %>%
#   rplot(shape = 19, colors = turbo(3), print_cor = T)

# dat %>%
#   dplyr::filter(Species == "Pinus sylvestris") %>%
#   dplyr::select_at(var_morph) %>%
#   correlate() %>%
#   # Re-arrange a correlation data frame to group highly correlated variables closer together
#   rearrange(method = "MDS", absolute = FALSE) %>%
#   shave() %>%
#   rplot(shape = 19, colors = turbo(3), print_cor = T)

# Remove highly correlated (highly correlation > 0.8):
# AR and Circ. correlated in pine but not in Betula
# "Solidity only correlated in Betula, but has a very narrow distribution
highly_correlated <- c("Perim.", "Area", "Round", "Solidity")

# Exclude highly correlated predictors
selected_predictors <- setdiff(var_morph, highly_correlated)

summary(dat[selected_predictors])

# # Transform variables ----------------------------------------------------------------

# log transform: size measurements
log_cols <- c("Feret", "MinFeret", "L", "AR")
dat[paste0(log_cols, "_log")] <- lapply(dat[log_cols], function(x) ifelse(x > 0, log(x), NA))

# logit transform: bounded ratio variables (Aspect ratio, Circularity) + L (rescaled to 0-1)
logit <- function(x) log(x / (1 - x))
logit_cols <- c("AR", "Circ.")
dat[paste0(logit_cols, "_logit")] <- lapply(dat[logit_cols], function(x) {
  ifelse(x > 0 & x < 1, logit(x), NA)
})

dat$L_logit <- with(dat, ifelse(L > 0 & L < 100, logit(L / 100), NA))

dat$L_sqrt <- with(dat, sqrt(L))

# ----------------------------------------------------------------------------
# Create histograms for each variable separately
dat_melted <- dat %>%
  dplyr::select_at(vars(Species, selected_predictors, ends_with("_log"), ends_with("_logit"), "L_sqrt")) %>%
  reshape2::melt()

# Create a histogram with ggplot
ggplot(dat_melted, aes(x = value, fill = variable)) +
  geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
  facet_wrap(variable ~ Species, scales = "free") +
  labs(title = "Histograms of Seed Measurements", x = "Value", y = "Frequency") +
  theme_minimal()

# Approach by Messier et al. 2010 ----------------------------------------------
# https://doi.org/10.1111/j.1461-0248.2010.01476.x

# evaluate models --------------------------------------------------------------
# Variables to test

all_variables <- c(selected_predictors, paste0(log_cols, "_log"), paste0(logit_cols, "_logit"), "L_logit", "L_sqrt")
all_variables <- sort(all_variables)

# Function to fit model and run DHARMa diagnostics
check_LMM_DHARMa <- function(var, sp, data) {
  # Fit the model
  model <- lmer(as.formula(paste(var, "~ 1 + (1 | SBD / uniTree)")),
                data = data, na.action = na.omit)

  # Simulate residuals using DHARMa
  simulationOutput <- simulateResiduals(fittedModel = model) #, n = 1000

  # Plot residual diagnostics
  plot(simulationOutput)
  title(sub = paste(var, "|", sp))
}

# # Apply DHARMa diagnostics to all variables
# for (species in unique(dat$Species)) {
#   species_data <- dat[dat$Species == species, ]
#   lapply(all_variables, function(v) check_LMM_DHARMa(v, species, species_data))
# }

## get variance partition -----------------------------------------------------

final_variables <- c("Feret", "MinFeret", "AR_log", "Circ._logit","L_log","a","b")

var_partition <- function(var, data, species) {
  model <- lmer(as.formula(paste(var, "~ 1 + (1 | SBD / uniTree)")),
                data = data, na.action = na.omit)

  if (isSingular(model)) {
    message("Warning: Model for ", species, ", variable '", var, "' has a boundary (singular) fit.")
  }

  var_comp <- as.data.frame(VarCorr(model))
  var_comp$proportion <- var_comp$vcov / sum(var_comp$vcov)
  var_comp$variable <- var
  var_comp$Species <- species

  return(var_comp)
}

all_results <- list()

for (species in unique(dat$Species)) {
  species_data <- dat[dat$Species == species, ]

  species_results <- do.call(rbind, lapply(final_variables, function(v) var_partition(v, species_data, species)))

  all_results[[species]] <- species_results
}

# Plot ---------------------------------------------------
combined_results <- do.call(rbind, all_results)

combined_results <- combined_results %>%
  mutate(
    grp = case_match(grp,
                     "SBD"         ~ "population",
                     "uniTree:SBD" ~ "tree",
                     "Residual"    ~ "within-tree and error",
                     .default      = grp
    ),
         grp = factor(grp, levels = c("population", "tree", "within-tree and error")))

unique(combined_results$variable)
combined_results <- combined_results %>%
  mutate(variable = case_match(variable,
                      "Feret" ~ "Length",
                      "MinFeret" ~ "Width",
                      "AR_log" ~ "Aspect Ratio",
                      "Circ._logit" ~ "Circularity",
                      "L_log" ~ "Lightness (L)",
                      "a" ~ "Red-Green Value (a)",
                      "b" ~ "Yellow-Blue Value (b)",
                      .default      = variable),
         variable = factor(variable, levels = c("Yellow-Blue Value (b)", "Red-Green Value (a)", "Lightness (L)", "Circularity","Aspect Ratio", "Width", "Length")))

pal_1 <- c("population" = "#1b5f8a", "tree" = "#7fcdbb", "within-tree and error" = "#e2ecf4")
(p2 = ggplot(combined_results, aes(x = variable, y = proportion, fill = fct_rev(grp))) +
  facet_grid(.~Species) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.7) +
  scale_fill_manual(name = "", values = pal_1, guide = guide_legend(reverse = TRUE))+
  labs(
    x = "",
    y = "Fraction of total variance explained by each scale"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey30", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 12),
    strip.text = element_text(face = "italic"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12)
  ) +
  coord_flip())

saveRDS(p2, "output/rds/plot_variance_partition_seed_traits.rds")

write.csv(combined_results, "output/csv/morphology_variance_partition.csv")
