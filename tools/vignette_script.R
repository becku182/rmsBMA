################################################################################
##
##  rmsBMA -- runnable companion script to the package vignette
##
##  Every R chunk from vignettes/rmsBMA.Rmd, in order, with the section
##  headings kept as comments so you can navigate and run it piece by piece.
##  Do not edit by hand: regenerate with tools/make_vignette_script.R.
##
##  HOW TO USE
##    - Run top to bottom, or step through a section at a time in RStudio.
##    - Fold sections with Alt+O (Windows) / Cmd+Alt+O (macOS).
##    - Takes roughly fifteen seconds with RUN_SLOW = FALSE.
##
##  RUN_SLOW = TRUE additionally builds the two large K = 17 model spaces by
##  exhaustive enumeration, and rebuilds the packaged `modelSpace` object from
##  scratch to check it still matches (about a minute in total). Those chunks
##  are eval = FALSE in the vignette, which is why its own build is quick.
##
################################################################################

RUN_SLOW <- FALSE

library(rmsBMA)

section <- function(x) message("\n=== ", x, " ===")


###############################################################################
## INTRODUCTION
###############################################################################
section("Introduction")


###############################################################################
## BAYESIAN MODEL AVERAGING: A THEORETICAL OVERVIEW
###############################################################################
section("Bayesian model averaging: A theoretical overview")


## -- Model setup and Bayesian estimation ---------------------------------
section("Model setup and Bayesian estimation")


## -- Bayesian model averaging --------------------------------------------
section("Bayesian model averaging")


## -- g-priors ------------------------------------------------------------
section("g-priors")


## -- Model priors --------------------------------------------------------
section("Model priors")


## -- Dilution priors -----------------------------------------------------
section("Dilution priors")


## -- Jointness measures --------------------------------------------------
section("Jointness measures")


## -- Extreme Bounds Analysis ---------------------------------------------
section("Extreme Bounds Analysis")


###############################################################################
## DATA PREPARATION
###############################################################################
section("Data preparation")

# (interactive / installation step, left commented)
# install.packages("rmsBMA")


library(rmsBMA)

# (interactive / installation step, left commented)
# ?migration_panel


migration_panel[1:10,1:7]


data <- data_preparation(migration_panel,
                         time = "Year_0",
                         id = "Pair_ID",
                         fixed_effects = TRUE,
                         effect = "twoway",
                         standardize = TRUE)
data[1:10,1:6]


data <- data_preparation(Trade_data,
                         standardize = TRUE)
data[1:10,1:6]

# (interactive / installation step, left commented)
# ?Trade_data


###############################################################################
## ESTIMATION OF THE MODEL SPACE
###############################################################################
section("Estimation of the model space")

if (RUN_SLOW) {
  modelSpace10 <- model_space(Trade_data, M = 10, g = "Benchmark")
}

if (RUN_SLOW) {
  modelSpace10_none <- model_space(Trade_data, M = 10, g = "None", HC = TRUE)
}

if (RUN_SLOW) {
  modelSpace <- model_space(Trade_data_small, M = 7, g = "UIP")
  # Same strings, but the packaged object carries a stale dim attribute on
  # its dimnames, so compare values and names rather than attributes.
  message("  rebuilt modelSpace matches the packaged object: ",
          isTRUE(all.equal(unname(modelSpace[[2]]),
                           unname(rmsBMA::modelSpace[[2]]))) &&
          identical(as.character(colnames(modelSpace[[2]])),
                    as.character(colnames(rmsBMA::modelSpace[[2]]))))
}


## -- Sampling the model space: MC^3 --------------------------------------
section("Sampling the model space: MC^3")

set.seed(1)
mc3Space <- model_space(Trade_data, mc3 = TRUE, draws = 50000, burn = 25000, g = "UIP")

mc3_results <- bma(mc3Space, EMS = 5, round = 4)
mc3_results[[1]]


## ---- Chain diagnostics ---------------------------------------------------
section("Chain diagnostics")

str(mc3Space[[6]][c("draws", "burn", "acceptance", "cor_pmp", "mean_size", "space_size")])


## ---- Reduced model spaces ------------------------------------------------
section("Reduced model spaces")


## ---- Checking the sampler against enumeration ----------------------------
section("Checking the sampler against enumeration")

exact_red   <- model_space(Trade_data_small, M = 5, g = "UIP")
exact_bma   <- bma(exact_red, EMS = 5, round = 6)

set.seed(4)
sampled_red <- model_space(Trade_data_small, M = 5, mc3 = TRUE,
                           draws = 50000, burn = 25000, g = "UIP")
sampled_bma <- bma(sampled_red, EMS = 5, round = 6)

round(cbind(Enumerated = exact_bma[[1]][, "PIP"],
            MC3        = sampled_bma[[1]][, "PIP"]), 4)


## ---- Extreme Bounds Analysis ---------------------------------------------
section("Extreme Bounds Analysis")

is.null(mc3_results[[3]])


###############################################################################
## PERFORMING BAYESIAN MODEL AVERAGING
###############################################################################
section("Performing Bayesian model averaging")


## -- Bayesian model averaging: The bma function --------------------------
section("Bayesian model averaging: The bma function")


bma_results <- bma(modelSpace, round = 3)


bma_results[[1]]


bma_results[[2]]


bma_results[[3]]


bma_results[[4]]


## -- Prior and posterior model probabilities -----------------------------
section("Prior and posterior model probabilities")


for_models <- model_pmp(bma_results)


for_models <- model_pmp(bma_results, top = 10)


size_graphs <- model_sizes(bma_results)


## -- Selecting the best models -------------------------------------------
section("Selecting the best models")


best_8_models <- best_models(bma_results, criterion = 1, best = 8)
best_8_models[[1]]


best_3_models <- best_models(bma_results, criterion = 2, best = 3)
best_3_models[[4]]


best_3_models <- best_models(bma_results, criterion = 2, best = 3)
grid::grid.draw(best_3_models[[6]])


## -- Calculating jointness measures --------------------------------------
section("Calculating jointness measures")


jointness(bma_results)[1:9,1:9]

jointness(bma_results, measure = "LS")[1:9,1:9]

jointness(bma_results, measure = "DW")[1:9,1:9]


## -- Visualizing model coefficients and posterior distributions ----------
section("Visualizing model coefficients and posterior distributions")


bin_sizes <- matrix(80, nrow = 11, ncol = 1)
coef_plots <- coef_hist(bma_results, BN = 1, num = bin_sizes)
coef_plots[[3]]


coef_plots2 <- coef_hist(bma_results, kernel = 1)
coef_plots2[[5]]


library(gridExtra)
grid.arrange(coef_plots[[3]], coef_plots[[5]], coef_plots2[[3]],
             coef_plots2[[5]], nrow = 2, ncol = 2)


coef_plots3 <- coef_hist(bma_results, weight = "beta", BN = 1, num = bin_sizes)
coef_plots3[[5]]


distPlots <- posterior_dens(bma_results, prior = "binomial")
grid.arrange(distPlots[[3]], distPlots[[5]], nrow = 2, ncol = 1)


###############################################################################
## CHANGES IN MODEL PRIORS
###############################################################################
section("Changes in model priors")


## -- Changing expected model size ----------------------------------------
section("Changing expected model size")


bma_results2 <- bma(modelSpace, round = 3, EMS = 2)


bma_results2[[4]]


size_graphs2 <- model_sizes(bma_results2)


model_graphs2 <- model_pmp(bma_results2)


bma_results2[[1]]


bma_results2[[2]]


jointness(bma_results2, measure = "HCGHM", rho = 0.5, round = 3)[1:9,1:9]


bma_results8 <- bma(modelSpace, round = 3, EMS = 8)
bma_results8[[4]]


size_graphs8 <- model_sizes(bma_results8)


model_graphs8 <- model_pmp(bma_results8)


bma_results8[[1]]


bma_results8[[2]]


jointness(bma_results8, measure = "HCGHM", rho = 0.5, round = 3)[1:9,1:9]


## -- Dilution prior ------------------------------------------------------
section("Dilution prior")


bma_results_dil <- bma(
  modelSpace = modelSpace,
  round       = 3,
  dilution    = 1
  )


size_graphs_dil <- model_sizes(bma_results_dil)


bma_results_dil01 <- bma(
  modelSpace = modelSpace,
  round       = 3,
  dilution    = 1,
  dil.Par     = 0.1
)
size_graphs_dil01 <- model_sizes(bma_results_dil01)


bma_results_dil2 <- bma(
  modelSpace = modelSpace,
  round       = 3,
  dilution    = 1,
  dil.Par     = 2
)
size_graphs_dil2 <- model_sizes(bma_results_dil2)


bma_results_dil2[[2]]


group_vec <- c(1,0,1,0,0,0,2,2,3,3)


cbind(modelSpace[[1]],group_vec)


par_vec <- c(0.8,0.6,0.4)


bma_results_dil3 <- bma(
  modelSpace = modelSpace,
  Narrative  = 1,
  Nar_vec    = group_vec,
  p          = par_vec,
  round      = 3
)
bma_results_dil3[[1]]


bma_results_dil3[[2]]


###############################################################################
## CONCLUDING REMARKS
###############################################################################
section("Concluding remarks")


###############################################################################
## REFERENCES
###############################################################################
section("References")


message("\n=== script complete ===")

