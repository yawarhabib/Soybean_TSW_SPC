# ==================================================================================================
#                  Script 02 - Phenotypic analysis of the Khabarovsk soybean panel
#                Field trials (2017–2024): thousand-seed weight (TSW)
#                                 and seed protein content (SPC)
# --------------------------------------------------------------------------------------------------
#                           Author :-  Yawar Habib, habib.yawr.10@gmail.com
#                              Center for Bio and Medical Technologies
#                    Skolkovo Institute of Science and Technology, Moscow, Russia
# --------------------------------------------------------------------------------------------------
# Input : one observation per genotype and season
#         data/raw/phenotype/TSW_SPC_season_means.csv
# Steps : data quality checks and descriptive statistics
#         outlier detection (Bernal-Vasquez et al., 2016)
#         genotype BLUEs
#         entry-mean and Cullis heritability
#         genetic correlation between TSW and SPC (bivariate model)
#         model diagnostics
# Output: data/processed/BLUE_tsw_GWAS.csv, BLUE_protein_GWAS.csv and Origin_Accessions.csv
# ==================================================================================================


### Clean the working environment
rm(list = ls())
graphics.off()
gc()

### Libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(sommer)
library(patchwork)
library(lme4)
library(emmeans)
library(here)


### Load the data from working environment
season_means <- read.csv(here("Paper_TSW_SPC/data/raw/phenotype/TSW_SPC_season_means.csv"), header = TRUE)
View(season_means)

str(season_means)


# Protein
season_means %>%
  filter(Genotype %in% c("Marinata", "Batya")) %>%
  dplyr::select(Year, Genotype, Protein.content.in.seeds) %>%
  pivot_wider(names_from = Genotype, values_from = Protein.content.in.seeds) %>%
  arrange(Year)

# TSW
season_means %>%
  filter(Genotype %in% c("Marinata", "Batya")) %>%
  dplyr::select(Year, Genotype, Weight.of.1000.seeds) %>%
  pivot_wider(names_from = Genotype, values_from = Weight.of.1000.seeds) %>%
  arrange(Year)

# Rename traits
season_means <- season_means %>% 
  rename(TSW = Weight.of.1000.seeds, Protein_content = Protein.content.in.seeds)

# Convert types
season_means <- season_means %>% 
  mutate(Year = as.factor(Year), Genotype = as.factor(Genotype))

#### Origin and genotype ID of each accession, needed for supplementary stuff.
Origin <- season_means %>%
  distinct(Genotype, ID_in_GneoFile, Origin_Eng)

nrow(Origin)

## Save the Origin information
write.csv(Origin, here("Paper_TSW_SPC/data/processed/Origin_Accessions.csv"), row.names = FALSE)

# Basic structure check
str(season_means)
nlevels(season_means$Genotype)   # number of unique genotypes
nlevels(season_means$Year)       # number of years

# Observations with a value per year and trait
season_means %>%
  group_by(Year) %>%
  summarise(TSW = sum(!is.na(TSW)), Protein_content = sum(!is.na(Protein_content)))


# TSW
TSW_data <- season_means %>%
  filter(!is.na(TSW))
TSW_data %>%
  summarise(n_obs = n(), n_genotypes = n_distinct(Genotype), n_years = n_distinct(Year))

# Protein
Protein_data <- season_means %>%
  filter(!is.na(Protein_content))
Protein_data %>%
  summarise(n_obs = n(), n_genotypes = n_distinct(Genotype), n_years = n_distinct(Year))


#### Outlier screen (Bernal-Vasquez, Utz & Piepho 2016, Theor Appl Genet 129:787-804, doi:10.1007/s00122-016-2666-6)
# Residuals of the genotype + year model are standardized by the median absolute deviation (rescaled x 1.4826,
# robust to the outliers themselves), turned into two-sided normal p-values and tested with the Bonferroni-Holm
# correction across all observations of a trait. Values with adjusted p < 0.05 are removed. One pass, no iteration.
# TSW
m_TSW_res <- lm(TSW ~ Genotype + Year, data = TSW_data)
TSW_data$res_std <- (residuals(m_TSW_res) - median(residuals(m_TSW_res))) / mad(residuals(m_TSW_res))
TSW_data$p_holm  <- p.adjust(2 * pnorm(-abs(TSW_data$res_std)), method = "holm")

TSW_outliers <- TSW_data %>% 
  filter(p_holm < 0.05)
TSW_outliers

TSW_data <- TSW_data %>% 
  filter(p_holm >= 0.05) %>% 
  dplyr::select(-res_std, -p_holm)

# Protein
m_Prot_res <- lm(Protein_content ~ Genotype + Year, data = Protein_data)
Protein_data$res_std <- (residuals(m_Prot_res) - median(residuals(m_Prot_res))) / mad(residuals(m_Prot_res))
Protein_data$p_holm  <- p.adjust(2 * pnorm(-abs(Protein_data$res_std)), method = "holm")

Protein_outliers <- Protein_data %>% filter(p_holm < 0.05)
Protein_outliers

Protein_data <- Protein_data %>% filter(p_holm >= 0.05) %>% dplyr::select(-res_std, -p_holm)


# ---- 1. Descriptive stats per trait, overall and per year ----
TSW_data %>%
  summarise(n = n(), mean = mean(TSW), sd = sd(TSW), min = min(TSW), max = max(TSW))

TSW_data %>%
  group_by(Year) %>%
  summarise(n = n(), mean = mean(TSW), sd = sd(TSW))

Protein_data %>%
  summarise(n = n(), mean = mean(Protein_content), sd = sd(Protein_content), min = min(Protein_content), max = max(Protein_content))

Protein_data %>%
  group_by(Year) %>%
  summarise(n = n(), mean = mean(Protein_content), sd = sd(Protein_content))

# ---- 2. Distribution of traits over the seasons ----

trait_colors <- c(TSW = "#1F77B4", Protein_content = "#CC0000")   

TSW_long <- TSW_data %>%
  dplyr::select(Year, Value = TSW) %>%
  mutate(Trait = "TSW")

Protein_long <- Protein_data %>%
  dplyr::select(Year, Value = Protein_content) %>%
  mutate(Trait = "Protein_content")

set.seed(1)
trait_long <- bind_rows(TSW_long, Protein_long) %>%
  mutate(Trait = factor(Trait, levels = c("TSW", "Protein_content")),
         x_point = as.numeric(Year) - 0.32 + runif(n(), -0.07, 0.07))

n_season <- trait_long %>%
  group_by(Trait, Year) %>%
  summarise(n = n(), .groups = "drop")

season_plot <- ggplot(trait_long, aes(as.numeric(Year), Value, group = Year, colour = Trait, fill = Trait)) +
  geom_point(aes(x = x_point), size = 0.8, alpha = 0.6) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.7) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.6, linewidth = 0.7) +
  geom_text(data = n_season, aes(as.numeric(Year), -Inf, label = paste0("n = ", n)),
            vjust = -0.6, size = 2.6, colour = "grey35") +
  scale_fill_manual(values = trait_colors) +
  scale_colour_manual(values = trait_colors) +
  scale_x_continuous(breaks = 1:nlevels(trait_long$Year), labels = levels(trait_long$Year)) +
  scale_y_continuous(expand = expansion(mult = c(0.14, 0.05))) +
  facet_wrap(~ Trait, ncol = 1, scales = "free_y", strip.position = "left",
             labeller = as_labeller(c(TSW = "Thousand-seed weight (g)",
                                      Protein_content = "Seed protein content (%)"))) +
  labs(x = "Growing season", y = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        strip.placement  = "outside",
        strip.background = element_blank(),
        strip.text       = element_text(size = 12),
        axis.line        = element_line(colour = "black", linewidth = 0.7),
        axis.ticks       = element_line(colour = "black", linewidth = 0.7),
        axis.text        = element_text(colour = "black"))

ggsave(here("Paper_TSW_SPC/results/02_phenotypes/Fig_traits_by_season.tiff"), season_plot,
       width = 17, height = 14, units = "cm", dpi = 600, compression = "lzw")


################ BLUEs per genotype

# ---- TSW: genotype and year fixed ----
m_TSW_blue <- lm(TSW ~ Genotype + Year, data = TSW_data)

summary(m_TSW_blue)

BLUE_TSW <- emmeans(m_TSW_blue, ~ Genotype) %>%
  as.data.frame() %>%
  rename(BLUE_TSW = emmean) %>%
  dplyr::select(Genotype, BLUE_TSW, SE)

# Genotype -> ID lookup, then ID next to Genotype
id_lookup <- TSW_data %>%
  dplyr::select(Genotype, ID_in_GneoFile) %>%
  distinct()

BLUE_TSW <- BLUE_TSW %>%
  left_join(id_lookup, by = "Genotype") %>%
  rename(ID = ID_in_GneoFile) %>%
  relocate(ID, .after = Genotype)

head(BLUE_TSW)
nrow(BLUE_TSW)
summary(BLUE_TSW$BLUE_TSW)

# genotypes without an ID, and genotypes with more than one ID across years
sum(is.na(BLUE_TSW$ID))

TSW_data %>%
  distinct(Genotype, ID_in_GneoFile) %>%
  count(Genotype) %>%
  filter(n > 1)


# ---- Protein: genotype and year fixed ----
m_Prot_blue <- lm(Protein_content ~ Genotype + Year, data = Protein_data)

summary(m_Prot_blue)

BLUE_Prot <- emmeans(m_Prot_blue, ~ Genotype) %>%
  as.data.frame() %>%
  rename(BLUE_Protein = emmean) %>%
  dplyr::select(Genotype, BLUE_Protein, SE)

# Genotype -> ID lookup, then ID next to Genotype
id_lookup_prot <- Protein_data %>%
  dplyr::select(Genotype, ID_in_GneoFile) %>%
  distinct()

BLUE_Prot <- BLUE_Prot %>%
  left_join(id_lookup_prot, by = "Genotype") %>%
  rename(ID = ID_in_GneoFile) %>%
  relocate(ID, .after = Genotype)

head(BLUE_Prot)
nrow(BLUE_Prot)
summary(BLUE_Prot$BLUE_Protein)

# genotypes without an ID, and genotypes with more than one ID across years
sum(is.na(BLUE_Prot$ID))

Protein_data %>%
  distinct(Genotype, ID_in_GneoFile) %>%
  count(Genotype) %>%
  filter(n > 1)

write.csv(BLUE_Prot, here("Paper_TSW_SPC/data/processed/BLUE_protein_GWAS.csv"), row.names = FALSE)
write.csv(BLUE_TSW,  here("Paper_TSW_SPC/data/processed/BLUE_tsw_GWAS.csv"),     row.names = FALSE)

################ Phenotypic correlation of the BLUEs

BLUE_combined <- BLUE_TSW %>%
  dplyr::select(Genotype, BLUE_TSW) %>%
  inner_join(BLUE_Prot %>% dplyr::select(Genotype, BLUE_Protein), by = "Genotype")

nrow(BLUE_combined)

cor.test(BLUE_combined$BLUE_TSW, BLUE_combined$BLUE_Protein)

#### Variance-components model: genotype RANDOM
m_TSW_vc  <- lmer(TSW ~ 1 + (1 | Genotype) + (1 | Year),
                  data = TSW_data, REML = TRUE)

m_Prot_vc <- lmer(Protein_content ~ 1 + (1 | Genotype) + (1 | Year),
                  data = Protein_data, REML = TRUE)

summary(m_TSW_vc)
summary(m_Prot_vc)

# Null models: drop the genotype random effect
m_TSW_null  <- lmer(TSW ~ 1 + (1 | Year),
                    data = TSW_data, REML = TRUE)

m_Prot_null <- lmer(Protein_content ~ 1 + (1 | Year),
                    data = Protein_data, REML = TRUE)

# REML log-likelihoods
logLik_TSW_full  <- as.numeric(logLik(m_TSW_vc,  REML = TRUE))
logLik_TSW_null  <- as.numeric(logLik(m_TSW_null, REML = TRUE))

logLik_Prot_full <- as.numeric(logLik(m_Prot_vc,  REML = TRUE))
logLik_Prot_null <- as.numeric(logLik(m_Prot_null, REML = TRUE))

# LRT statistic
LRT_TSW  <- -2 * (logLik_TSW_null  - logLik_TSW_full)
LRT_Prot <- -2 * (logLik_Prot_null - logLik_Prot_full)

# Standard chi-square (1 df) p-value
p_chi_TSW  <- pchisq(LRT_TSW,  df = 1, lower.tail = FALSE)
p_chi_Prot <- pchisq(LRT_Prot, df = 1, lower.tail = FALSE)

# Boundary-corrected p-value (50:50 mixture of chi^2_0 and chi^2_1)
p_mix_TSW  <- 0.5 * pchisq(LRT_TSW,  df = 1, lower.tail = FALSE)
p_mix_Prot <- 0.5 * pchisq(LRT_Prot, df = 1, lower.tail = FALSE)

LRT_results <- tibble::tibble(
  Trait        = c("TSW", "Protein"),
  LRT_stat     = c(LRT_TSW, LRT_Prot),
  p_chi1df     = c(p_chi_TSW, p_chi_Prot),
  p_boundary   = c(p_mix_TSW, p_mix_Prot)
)

LRT_results


# Harmonic mean of years per genotype, per trait
harmonic_mean <- function(x) length(x) / sum(1 / x)

r_TSW <- TSW_data %>%
  count(Genotype, name = "n_years") %>%
  pull(n_years) %>%
  harmonic_mean()

r_Prot <- Protein_data %>%
  count(Genotype, name = "n_years") %>%
  pull(n_years) %>%
  harmonic_mean()

# Variance components
vc_TSW  <- as.data.frame(VarCorr(m_TSW_vc))
vc_Prot <- as.data.frame(VarCorr(m_Prot_vc))

sigma2_G_TSW  <- vc_TSW [vc_TSW $grp == "Genotype", "vcov"]
sigma2_e_TSW  <- vc_TSW [vc_TSW $grp == "Residual", "vcov"]
sigma2_G_Prot <- vc_Prot[vc_Prot$grp == "Genotype", "vcov"]
sigma2_e_Prot <- vc_Prot[vc_Prot$grp == "Residual", "vcov"]

mean_TSW  <- mean(TSW_data$TSW,             na.rm = TRUE)
mean_Prot <- mean(Protein_data$Protein_content, na.rm = TRUE)

# Heritability, CVg, CVr
results <- tibble::tibble(
  Trait = c("TSW", "Protein"),
  Mean  = c(mean_TSW, mean_Prot),
  sigma2_G = c(sigma2_G_TSW, sigma2_G_Prot),
  sigma2_e = c(sigma2_e_TSW, sigma2_e_Prot),
  r_bar    = c(r_TSW, r_Prot),
  H2       = c(sigma2_G_TSW  / (sigma2_G_TSW  + sigma2_e_TSW  / r_TSW),
               sigma2_G_Prot / (sigma2_G_Prot + sigma2_e_Prot / r_Prot)),
  CVg = c(sqrt(sigma2_G_TSW)  / mean_TSW  * 100,
          sqrt(sigma2_G_Prot) / mean_Prot * 100),
  CVr = c(sqrt(sigma2_e_TSW)  / mean_TSW  * 100,
          sqrt(sigma2_e_Prot) / mean_Prot * 100)
)

results

##
#### Cullis Heritability
# ---- Cullis heritability from an lme4 fit (Cullis et al. 2006; lme4 computation as in Schmidt et al. 2019) ----
# Builds the mixed-model coefficient matrix, inverts it, and takes the prediction-error
# variance of the genotype BLUPs. Random terms with zero variance are dropped.
cullisH2 <- function(m, gterm = "Genotype", tol = 1e-8) {
  if (!inherits(m, "merMod")) return(NA_real_)
  vc <- as.data.frame(VarCorr(m))
  gv <- function(g) vc$vcov[vc$grp == g & is.na(vc$var2)][1]
  vg <- gv(gterm)
  if (is.na(vg) || vg <= tol) return(NA_real_)
  s2e  <- sigma(m)^2
  X    <- getME(m, "X")
  Ztl  <- getME(m, "Ztlist")
  grp  <- sub(".(Intercept)", "", names(Ztl), fixed = TRUE)
  keep <- which(vapply(grp, function(g) { v <- gv(g); !is.na(v) && v > tol }, logical(1)))
  if (!(which(grp == gterm)[1] %in% keep)) return(NA_real_)
  Ztl  <- Ztl[keep]; grp <- grp[keep]
  Z    <- t(do.call(rbind, Ztl))
  vpar <- unlist(lapply(seq_along(Ztl), function(i) rep(gv(grp[i]), nrow(Ztl[[i]]))))
  C <- rbind(cbind(crossprod(X),    crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + s2e * Diagonal(x = 1 / vpar)))
  Cinv <- try(solve(C), silent = TRUE)
  if (inherits(Cinv, "try-error")) return(NA_real_)
  off <- c(0, cumsum(vapply(Ztl, nrow, numeric(1))))
  gi  <- which(grp == gterm)[1]
  idx <- ncol(X) + ((off[gi] + 1):off[gi + 1])
  PEV <- s2e * as.matrix(Cinv[idx, idx, drop = FALSE])
  n   <- nrow(PEV)
  as.numeric(1 - (2 * (n * sum(diag(PEV)) - sum(PEV)) / (n * (n - 1))) / (2 * vg))
}

# ---- Apply to both traits ----
cullis_TSW  <- cullisH2(m_TSW_vc,  gterm = "Genotype")
cullis_Prot <- cullisH2(m_Prot_vc, gterm = "Genotype")

cullis_TSW
cullis_Prot

# ---- Combine with standard H2 results table ----
results_H2 <- results %>%
  mutate(H2_Cullis = c(cullis_TSW, cullis_Prot))

results_H2

################ Genetic correlation between TSW and protein, bivariate model

# Screened data of both traits, rows where at least one trait is present
biv_data <- full_join(TSW_data %>% dplyr::select(Genotype, Year, TSW),
                      Protein_data %>% dplyr::select(Genotype, Year, Protein_content),
                      by = c("Genotype", "Year"))

# Long format for mmes; record pairs the two traits of the same genotype x season, empty cells are dropped by mmes
DTL <- stackTraits(biv_data, traits = c("TSW", "Protein_content"), keep = c("Genotype", "Year"))

# Trait means fixed; genotype and residual unstructured between traits, year with one variance per trait
ansm <- mmes(value ~ trait,
             random = ~ vsm(usm(trait), ism(Genotype)) + vsm(dsm(trait), ism(Year)),
             rcov   = ~ vsm(usm(trait), ism(record)),
             data   = DTL,
             control = mmes_control(diagnostics = list(verbose = FALSE)))

# Variance components: V1-V3 genotype (TSW, covariance, protein), V4-V5 year, V6-V8 residual
vc_table <- summary(ansm)$varcomp
vc_table

# Genetic and residual correlation between the traits, delta-method SE
rg_TSW_Prot <- vpredict(ansm, rg ~ V2 / sqrt(V1 * V3))
rg_TSW_Prot

re_TSW_Prot <- vpredict(ansm, re ~ V7 / sqrt(V6 * V8))
re_TSW_Prot


############### Supplementary figure: residual diagnostics of the BLUE models

# Standardized residuals of both models; genotypes seen in one season are fitted exactly (leverage 1) and give NaN
diag_df <- bind_rows(
  data.frame(Trait = "TSW",             fitted = fitted(m_TSW_blue),  resid_std = rstandard(m_TSW_blue)),
  data.frame(Trait = "Protein_content", fitted = fitted(m_Prot_blue), resid_std = rstandard(m_Prot_blue))) %>%
  mutate(Trait = factor(Trait, levels = c("TSW", "Protein_content"))) %>%
  filter(is.finite(resid_std))

trait_labels <- c(TSW = "Thousand-seed weight (g)", Protein_content = "Seed protein content (%)")

p_fit <- ggplot(diag_df, aes(fitted, resid_std, colour = Trait)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Trait, scales = "free_x", labeller = as_labeller(trait_labels)) +
  labs(x = "Fitted values", y = "Standardized residuals")

p_qq <- ggplot(diag_df, aes(sample = resid_std, colour = Trait)) +
  stat_qq(size = 0.8, alpha = 0.6) +
  stat_qq_line(colour = "black") +
  facet_wrap(~ Trait, labeller = as_labeller(trait_labels)) +
  labs(x = "Theoretical quantiles", y = "Standardized residuals")

# & applies the colour scale and the theme to every panel of the stacked figure
diag_plot <- (p_fit / p_qq) &
  scale_colour_manual(values = trait_colors) &
  theme_classic(base_size = 12) &
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.line  = element_line(colour = "black", linewidth = 0.7),
        axis.ticks = element_line(colour = "black", linewidth = 0.7),
        axis.text  = element_text(colour = "black"))

ggsave(here("Paper_TSW_SPC/results/02_phenotypes/FigS_residual_diagnostics.tiff"), diag_plot,
       width = 17, height = 14, units = "cm", dpi = 600, compression = "lzw")
