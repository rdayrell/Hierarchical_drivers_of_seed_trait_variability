# Seed quality, population maps and bioclimatic effects -------------------
# Run from the project root (open seed_itv.Rproj first).

dir.create("output/rds", recursive = TRUE, showWarnings = FALSE)
dir.create("output/csv", recursive = TRUE, showWarnings = FALSE)
dir.create("output/imgs", recursive = TRUE, showWarnings = FALSE)


# Packages --------------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(patchwork)
library(caret)
library(ggpattern)

# Import data ------------------------------------------------------------------

original_data <- read.csv("data/embryoless_seeds.csv")

# ------------------------------------------------------------------------------
# Maps -------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Fig 2a

library(terra)
library(tidyterra)
library(sf)
library(geodata)

## add info on populations with morphological analysis
morph = data.frame("SBD" = as.integer(c("892915", "895732", "950323", "1010271", "912563", "912585", "912596", "925037", "925509")), "morph" = rep(1, 9)) 
# "893048" removed due to genetic analysis: contains 1 pentaploid, 1 diploid

df_location = original_data %>%
    full_join(., morph, by = c("SBD")) %>%
  select(Species, SBD, n_trees, lat_dd, long_dd, Geographical.Location, morph) %>%
  distinct() %>%
  mutate(fil = ifelse(is.na(morph), "x", "y"),
         Species = str_replace_all(Species, " ", "\n"))

crdref <- "+proj=longlat +datum=WGS84"

w <- world(path="data", resolution = 3, version = "3.6")
uk <- w[w$GID_0=="GBR", ]
crs(uk) <- crdref

lonlat <- cbind(df_location$long_dd, df_location$lat_dd)

filt_pts <- vect(lonlat, crs=crdref, atts = df_location)

p0 <- filt_pts %>%
ggplot() +
  facet_wrap(vars(Species), nrow = 1, strip.position = "bottom", dir = "h") + #dir = "v"
  geom_spatvector(data = uk, color = "grey82", fill ="grey82") +
  geom_spatvector(aes(colour = fil, fill = fil), size=2, show.legend = NA, shape = 21) +
  scale_colour_manual(values = c("#1f5d8b", "darkorange3")) +
  scale_fill_manual(values = c(adjustcolor("#1f5d8b", alpha.f = 0.75), adjustcolor("darkorange3", alpha.f = 0.75))) +
  coord_sf() +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "italic", size = 11),
        panel.grid = element_blank(),
        strip.background = element_blank(), #element_rect(fill="white"),
        panel.border = element_blank(),
        panel.background = element_blank(), #element_rect(fill = "white"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0.01, "lines"))

p0

saveRDS(p0, "output/rds/plot_map_seed_empty.rds")

#-------------------------------------------------------------------------------
# Bioclim names ----------------------------------------------------------------

library(stringr)
bioclim_names = read.csv("data/bioclim_names.csv")

bioclim_names <- bioclim_names %>%
  mutate(
    term_name = str_replace_all(term_name, c(
      "m-2" = "m^{-2}",
      "month-1" = "month^{-1}",
      "year-1" = "year^{-1}",
      "°C" = "degree*C",
      "/100" = " 10^{-2}",
      "mean annual air temperature" = 'mean annual air,temperature',
      "mean daily mean air temperature of the" = 'mean air temperature,of the',
      "mean daily mean air temperatures of the" = 'mean air temperatures,of the',
      "precipitation amount of the" = 'precipitation amount,of the',
      'seasonality ' = 'seasonality,'
    )),
    term_name = ifelse(
      grepl(",", term_name),
      paste0(
        'expression(atop(',
        gsub(",", ', ', term_name),
        '))'
      ),
      paste0("expression(", term_name, ")")
    ),
    term_name = str_replace_all(term_name, " ", "~")
  )


#-------------------------------------------------------------------------------
# Prepare and inspect data for analyses ----------------------------------------
#-------------------------------------------------------------------------------

# removes population that is included in morphological analyses but
# doesn't meet inclusion criteria for embryoless seed analyses
df <- subset(original_data, !is.na(n_populations)) 

glimpse(df)

# Table for Supplementary Material (Notes S1)
df %>%
  select("Species","SBD","n_populations","n_trees") %>%
  distinct() %>%
  group_by(Species, n_populations) %>%
  summarise(
    sum = sum(n_trees),
    min = min(n_trees),
    max = max(n_trees)
  )

df %>%
  select("Species","n_populations") %>%
  distinct() %>%
  ungroup() %>%
  summarise(
    n_pop = sum(n_populations)
  )

df %>%
  select("Species","SBD","n_trees") %>%
  distinct() %>%
  ungroup() %>%
  summarise(
    tot = sum(n_trees)
  )


# Apply logit transformation
# https://doi.org/10.1890/10-0340.1
df$per_empty <- df$Empty / df$Sown
df$logit_empty <- car::logit(df$per_empty, percents = FALSE, adjust = 0.005)

# Function to invert logit for plotting
inv.logit <- function(f,a) {
  a <- (1-2*a)
  (a*(1+exp(f))+(exp(f)-1))/(2*a*(1+exp(f)))
}


all_plot_list <- list() # List to store plots for all species
models_fixed <- data.frame()
results_list <- list()  # List to store variance results
correlation_data_all <- data.frame(Variable = character(),
                                   Component = character(),
                                   Correlation = numeric(),
                                   stringsAsFactors = FALSE)

for (species in unique(df$Species)) {

  bioclim_data <- df %>%
  filter(Species == species) %>%
  select(Species, SBD, test_tree_number, starts_with("CHELSA_bio"), logit_empty, per_empty, Empty, Sown) %>%
  distinct() %>%
  mutate(across(starts_with("CHELSA_bio"), ~ ., .names = "unscaled_{.col}")) %>%
  mutate(across(starts_with("CHELSA_bio"), scale))

  scaling_parameters <- bioclim_data %>%
    select(starts_with("CHELSA_bio")) %>%
    reframe(across(everything(), ~ list(attr(.x, "scaled:scale"), attr(.x, "scaled:center"))))

  cor_data = bioclim_data %>%
    select(-Species, -test_tree_number, -logit_empty, -per_empty, -Empty, -Sown, -starts_with("unscaled_")) %>%
    distinct() %>%
    select(-SBD)

  # Correlation analysis to identify collinear predictors
  cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

  highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.7, names = TRUE)

  selected_predictors <- setdiff(paste0("CHELSA_bio", 1:19), highly_correlated)

  formula <- as.formula(
    paste("logit_empty ~", paste(selected_predictors, collapse = " + "), "+ (1 | SBD)")
  )
  model <- lmer(formula, data = bioclim_data)
  plot(model)

  model_summary <- summary(model)

  significant_terms <- model_summary$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(`Pr(>|t|)` < 0.05 & term != "(Intercept)") %>%
    select(term, Estimate,`Pr(>|t|)`) %>%
    left_join(., bioclim_names)

  if (nrow(significant_terms) > 0) {
    species_plot_list <- list()  # Store plots for the current species


    for (i in 1:nrow(significant_terms)) {
      term <- significant_terms$term[i]
      term_name <- significant_terms$term_name[i]
      estimate <- round(significant_terms$Estimate[i], 3)

      term_min <- min(bioclim_data[term], na.rm = TRUE) - 0.05
      term_max <- max(bioclim_data[term], na.rm = TRUE) + 0.1
      term_pred <- ggpredict(
        model,
        terms = paste0(term, " [", term_min, ":", term_max, " by = .1]"),
        interval = "confidence"
      )

      scale_factor <- scaling_parameters[[term]][[1]]
      center_factor <- scaling_parameters[[term]][[2]]

      term_pred <- term_pred %>%
        group_by(x) %>%
        mutate(
          original = inv.logit(predicted, a = 0.005),
          conf.low = inv.logit(conf.low, a = 0.005),
          conf.high = inv.logit(conf.high, a = 0.005),
          x_unscaled =  x * scale_factor + center_factor
        )

      (p <- ggplot() +
        geom_line(data = term_pred, aes(x = x_unscaled, y = original), color = "black") +
        geom_ribbon(data = term_pred, aes(x = x_unscaled, ymin = conf.low, ymax = conf.high), fill = "black", alpha = 0.2) +
        geom_point(data = bioclim_data, aes(x = !!sym(paste0("unscaled_",term)), y = per_empty, color = factor(SBD)), alpha = 0.5, size = 2) +
        labs(
          x = eval(parse(text = term_name)),
          y = "Proportion of embryoless seeds",
          subtitle = bquote(atop(italic(.(species)), "Coefficient: " * .(estimate)))
        ) +
        scale_y_continuous(limits = c(-0.005,1)) +
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.subtitle = element_text(size = 11, face = "plain"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10)
        )
      )

      species_plot_list[[i]] <- p  # Add plot to the species list

    }

    all_plot_list[[species]] <- species_plot_list
  }

  model_fixed <- model_summary$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(sig = ifelse(`Pr(>|t|)` < 0.05, 1, 0),
           Species = species)

  models_fixed <- rbind(models_fixed, model_fixed)

  data <- subset(df, Species == species)
  bioclim_vars <- grep("CHELSA_bio", names(data), value = TRUE)

  for (selected_predictor in selected_predictors) {
    for (bioclim_var in bioclim_vars) {
      correlation_value <- cor(data[[selected_predictor]], data[[bioclim_var]], use = "pairwise.complete.obs")
      correlation_data_all <- rbind(correlation_data_all, data.frame(
        Species = species,
        Variable = bioclim_var,
        Component = selected_predictor,
        Correlation = correlation_value
      ))
    }
  }
}

# ------------------------

all_plots_flat <- unlist(all_plot_list, recursive = FALSE)


(p2 <- all_plots_combined <- wrap_plots(all_plots_flat, ncol = 3) +
    plot_layout(axes = "collect_y")
    )


saveRDS(p2, "output/rds/plot_fixed_effects_embryoless_seeds.rds")

# Table fixed effects models
write.csv(models_fixed, "output/csv/fixed_effects_embryoless_seeds.csv", row.names = FALSE)

# ------------------------

correlation_data_all$Correlation_lab <- round(correlation_data_all$Correlation, 2)
correlation_data_all <- correlation_data_all %>%
  mutate(Variable = str_replace(Variable, "CHELSA_bio", ""),
         Variable = str_pad(Variable, width = 2, pad = "0"),
         Variable = fct_rev(as.factor(Variable)),
         Component = str_replace(Component, "CHELSA_bio", ""),
         Component = str_pad(Component, width = 2, pad = "0"),
         Correlation_lab = round(Correlation, 2),
         Correlation_lab = ifelse(Variable == Component, NA, Correlation_lab)
         )

p3 <- ggplot(correlation_data_all, aes(x = Component, y = Variable, fill = Correlation_lab)) +
  facet_wrap(~ Species, scales = "free", ncol = 2) +
  geom_tile(color = "white") +
  geom_tile_pattern(aes(fill = Correlation_lab),
                    pattern = ifelse(is.na(correlation_data_all$Correlation_lab), "stripe", "none"),
                    pattern_density = 0.1, pattern_spacing = 0.01,
                    pattern_fill = "grey60") +
  scale_fill_gradient2(low = "red4", high = "green4", mid = "white", midpoint = 0, na.value = "grey90",
                       limit = c(-1, 1), space = "Lab", name = "Correlation") +
    geom_text(aes(label = Correlation_lab), size = 3) +
  labs(x = "Tested Bioclim Variable",
       y = "Bioclim Variable") +
  theme_minimal() +
  theme(strip.text = element_text(face = "italic"))


ggsave("output/imgs/figure_s1_bioclim_correlation.png", p3, width = 8, height = 14, units = "in", bg = "white")


#-------------------------------------------------------------------------------
# Variance decomposition
# Approach by Messier et al. 2010 ----------------------------------------------
# doi: https://doi.org/10.1111/j.1461-0248.2010.01476.x
#-------------------------------------------------------------------------------

# 
# # https://doi.org/10.1890/10-0340.1

results_list <- list()

for (species in unique(df$Species)) {

  bioclim_data <- df %>%
    filter(Species == species) %>%
    select(Species, SBD, test_tree_number, logit_empty) %>%
    distinct() %>%
    droplevels()

  formula <- as.formula(
    paste("logit_empty ~ 1 + (1 | SBD)")
  )
  model <- lmer(formula, data = bioclim_data)

  # Variance partition after correlation
  var_comp <- as.data.frame(VarCorr(model))
  random_variance <- var_comp$vcov[1]       # Random effect variance (SBD)
  residual_variance <- var_comp$vcov[2]     # Residual variance

  total_variance <- random_variance + residual_variance
  prop_random <- random_variance / total_variance
  prop_residual <- residual_variance / total_variance

  results_list[[species]] <- data.frame(
    Species = species,
    Component = c("Random Effects", "Residuals"),
    Variance = c(random_variance, residual_variance),
    Proportion = c(prop_random, prop_residual)
  )
}

# ------------------------

variance_results <- do.call(rbind, results_list)
variance_plot <- variance_results %>%
  mutate(Species = fct_rev(as.factor(Species)),
         Component = dplyr::recode(Component,
                            "Random Effects" = "population",
                            "Residuals" = "tree and error"),
         Component = factor(Component, levels = c("population", "tree and error"))
  )

pal_1 <- c("population" = "#1b5f8a", "tree and error" = "#7fcdbb") #"bioclim" = "#3b4a5d",
(p1 <- ggplot(variance_plot, aes(x = Species, y = Proportion, fill = fct_rev(Component))) +
    geom_bar(stat = "identity",position = "stack", color = "black",  width = 0.7) +
    labs(y = "Fraction of total variance of embryoless\nseeds explained by each scale", x = "") +
    scale_fill_manual(name = "", values = pal_1, guide = guide_legend(reverse = TRUE)) + #Variance Components
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey30", linewidth = 0.6),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.justification = "right",
      legend.title = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.text.y = element_text(face = "italic")
    )+
    coord_flip())

saveRDS(p1, "output/rds/plot_variance_partition_seed_quality_random_only.rds")

# Table random effects models
write.csv(variance_plot, "output/csv/variance_partition_embryoless_seed.csv", row.names = FALSE)
