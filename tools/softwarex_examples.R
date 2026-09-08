# ---------------------------------------------------------------------------
# Reproduces every example and figure in
#
#   Beck, K. "rmsBMA: Bayesian model averaging over reduced model spaces",
#   SoftwareX.
#
# Run from any directory; figures are written to the working directory as
# fig_sizes.pdf, fig_mc3.pdf and fig_sim_pmp.pdf. Every number quoted in the
# paper is printed to the console.
#
# Takes two to three minutes. Requires rmsBMA >= 0.2.0 and ggplot2.
# ---------------------------------------------------------------------------

library(rmsBMA)
library(ggplot2)

data(Trade_data)

# ===========================================================================
# 1. Enumeration, and the graphical output of model_sizes  (Figure 2)
# ===========================================================================

cat("== enumerated model space, K = 17 ==\n")
t_enum <- system.time(ms_full <- model_space(Trade_data))[["elapsed"]]
b_full <- bma(ms_full)

cat("models         :", format(ms_full$MS, big.mark = ","), "\n")
cat("elapsed        :", round(t_enum, 1), "s\n\n")

# model_sizes draws as a side effect AND returns the assembled plot as its
# third element. Saving that element avoids the blank first page a bare
# pdf() device would produce.
ggsave("fig_sizes.pdf", model_sizes(b_full)[[3]], width = 7, height = 7.2)

# A reduced space over the same data
t_red <- system.time(ms_red <- model_space(Trade_data, M = 5))[["elapsed"]]
cat("with M = 5     :", format(ms_red$MS, big.mark = ","), "models in",
    round(t_red, 1), "s\n\n")

# ===========================================================================
# 2. MC^3 against exact enumeration  (Figure 3)
# ===========================================================================

cat("== MC3 on the same space ==\n")
set.seed(20260906)
t_mc3  <- system.time(
  ms_mc3 <- model_space(Trade_data, mc3 = TRUE, draws = 100000, burn = 20000)
)[["elapsed"]]
b_mc3 <- bma(ms_mc3)

cat("elapsed        :", round(t_mc3, 1), "s\n")
cat("distinct models:", format(ms_mc3$MS, big.mark = ","), "\n")
cat("acceptance     :", round(ms_mc3$info$acceptance, 3), "\n")
cat("cor_pmp        :", round(ms_mc3$info$cor_pmp, 4), "\n")

t1 <- summary(b_full)$table
t2 <- summary(b_mc3)$table
k  <- setdiff(intersect(rownames(t1), rownames(t2)), "CONST")
d  <- data.frame(enum = t1[k, "PIP"], mc3 = t2[k, "PIP"])

cat("cor(PIP)       :", round(cor(d$enum, d$mc3), 5), "\n")
cat("max |diff|     :", round(max(abs(d$enum - d$mc3)), 4), "\n\n")

p_mc3 <- ggplot(d, aes(enum, mc3)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 2.4, colour = "darkred") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "PIP, full enumeration (131,072 models)",
       y = expression(paste("PIP, ", MC^3, " (100,000 draws)"))) +
  theme_bw(base_size = 12)
ggsave("fig_mc3.pdf", p_mc3, width = 5.2, height = 5.0)

# ===========================================================================
# 3. More regressors than observations  (Figure 4)
# ===========================================================================

cat("== K > n: n = 20, K = 30, M = 5 ==\n")
set.seed(20260907)

n <- 20; K <- 30; M <- 5
true_names <- sprintf("x%02d", 1:4)
true_beta  <- c(1.5, -1.2, 1.0, -0.8)

X <- matrix(rnorm(n * K), nrow = n, ncol = K)
colnames(X) <- sprintf("x%02d", seq_len(K))
y   <- 2 + X[, 1:4] %*% true_beta + rnorm(n, sd = 0.5)
sim <- cbind(y = y, X)

cat("full space     :", format(2^K, big.mark = ",", scientific = FALSE), "models\n")
cat("largest estimable model:", n - 1, "regressors\n")
cat("admissible     :", format(sum(choose(K, 0:M)), big.mark = ","), "\n")
cat("uniform prior PMP:", signif(1 / sum(choose(K, 0:M)), 2), "\n")

t_sim <- system.time(ms_sim <- model_space(sim, M = M, g = "UIP"))[["elapsed"]]
res   <- bma(ms_sim, EMS = 3, round = 4)
cat("elapsed        :", round(t_sim, 1), "s\n\n")

tab   <- res[[1]]
pip   <- tab[rownames(tab) != "CONST", "PIP"]
noise <- pip[setdiff(names(pip), true_names)]

print(round(data.frame(beta = true_beta,
                       PIP  = pip[true_names],
                       PM   = tab[true_names, "PM"]), 3))
cat("\nnoise (26): highest PIP", round(max(noise), 4),
    " median", round(median(noise), 4), "\n")

cat("\n-- the five highest-probability models --\n")
bm <- best_models(res, best = 5, round = 4)
comp <- bm[[1]]
keep <- rownames(comp) %in% c("Const", sprintf("x%02d", 1:5), "PMP", "R^2") |
        rowSums(abs(as.matrix(comp))) > 0
print(comp[keep, , drop = FALSE])

cat("\n-- prior and posterior mean model size (bound M =", M, ") --\n")
print(round(res[[4]], 3))

ggsave("fig_sim_pmp.pdf",
       model_pmp(res, top = 5, type = "histogram")[[3]],
       width = 7, height = 4.4)

cat("\nfigures written: fig_sizes.pdf, fig_mc3.pdf, fig_sim_pmp.pdf\n")
