# ---------------- Libraries ----------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)

# ---------------- Parameters ----------------
I <- 120        # sample size (no. of matched pairs)
K <- 30         # total number of outcomes
no.sims <- 5000 # number of iterations in the simulation
mu.grid <- seq(0.3, 0.6, by = 0.05)  # effect sizes
K1.grid <- c(1, 3, 5, 10)            # number of affected outcomes

t.test.pval.func <- function(x) t.test(x)$p.value

# ---------------- Simulation ----------------
results <- data.frame()

for (K1 in K1.grid) {
  for (mu in mu.grid) {
    
    reject.vec.two.team.global <- numeric(no.sims)
    reject.vec.three.team.global <- numeric(no.sims)
    reject.vec.four.team.global <- numeric(no.sims)
    
    reject.vec.two.team.rep <- numeric(no.sims)
    reject.vec.three.team.rep <- numeric(no.sims)
    reject.vec.four.team.rep <- numeric(no.sims)
    
    for (s in 1:no.sims) {
      # Generate data
      y.mat <- cbind(matrix(rnorm(I * K1, mean = mu, sd = 1), ncol = K1),
                     matrix(rnorm(I * (K - K1), mean = 0, sd = 1), ncol = (K - K1)))
      
      ## ----- Two team cross-screening---------------
      y.mat.fh <- y.mat[1:(I/2), ]
      y.mat.sh <- y.mat[(I/2 + 1):I, ]
      
      pvals.fh <- apply(y.mat.fh, 2, t.test.pval.func)
      pvals.sh <- apply(y.mat.sh, 2, t.test.pval.func)
      
      true.rej.num.global2 <- 0
      true.rej.num.rep2 <- 0
      
      for (f in 1:K1) {
        cond_two <- ((pvals.fh[f] < .025 & pvals.sh[f] < .025 / max(1, sum(pvals.fh < .025))) |
                       (pvals.sh[f] < .025 & pvals.fh[f] < .025 / max(1, sum(pvals.sh < .025))))
        
        if (cond_two) {
          true.rej.num.global2 <- true.rej.num.global2 + 1
        }
        
        if ((pvals.fh[f] < .025 & pvals.sh[f] < .025 / max(1, sum(pvals.fh < .025))) &
            (pvals.sh[f] < .025 & pvals.fh[f] < .025 / max(1, sum(pvals.sh < .025)))) {
          true.rej.num.rep2 <- true.rej.num.rep2 + 1
        }
      }
      
      reject.vec.two.team.global[s] <- true.rej.num.global2
      reject.vec.two.team.rep[s] <- true.rej.num.rep2
      
      ## ----- Three team cross-screening (1/3 for planning) -----
      y.mat.ft <- y.mat[1:(I/3), ]
      y.mat.st <- y.mat[(I/3+1):(2*I/3), ]
      y.mat.tt <- y.mat[(2*I/3+1):I, ]
      
      pvals.ft <- apply(y.mat.ft, 2, t.test.pval.func)
      pvals.st <- apply(y.mat.st, 2, t.test.pval.func)
      pvals.tt <- apply(y.mat.tt, 2, t.test.pval.func)
      
      pvals.ftst <- apply(rbind(y.mat.ft, y.mat.st), 2, t.test.pval.func)
      pvals.fttt <- apply(rbind(y.mat.ft, y.mat.tt), 2, t.test.pval.func)
      pvals.sttt <- apply(rbind(y.mat.st, y.mat.tt), 2, t.test.pval.func)
      
      true.rej.num.global3 <- 0
      true.rej.num.rep3 <- 0
      
      for (f in 1:K1) {
        cond1 <- (pvals.ft[f] < .0167 & pvals.sttt[f] < .0167 / max(1, sum(pvals.ft < .0167)))
        cond2 <- (pvals.st[f] < .0167 & pvals.fttt[f] < .0167 / max(1, sum(pvals.st < .0167)))
        cond3 <- (pvals.tt[f] < .0167 & pvals.ftst[f] < .0167 / max(1, sum(pvals.tt < .0167)))
        
        if (cond1 | cond2 | cond3) {
          true.rej.num.global3 <- true.rej.num.global3 + 1
        }
        if (cond1 & cond2 & cond3) {
          true.rej.num.rep3 <- true.rej.num.rep3 + 1
        }
      }
      
      reject.vec.three.team.global[s] <- true.rej.num.global3
      reject.vec.three.team.rep[s] <- true.rej.num.rep3
      
      ## ----- Four team cross-screening (1/4 for planning) -----
      y.mat.1f <- y.mat[1:(I/4), ]
      y.mat.2f <- y.mat[(I/4+1):(2*I/4), ]
      y.mat.3f <- y.mat[(2*I/4+1):(3*I/4), ]
      y.mat.4f <- y.mat[(3*I/4+1):I, ]
      
      pvals.1f <- apply(y.mat.1f, 2, t.test.pval.func)
      pvals.2f <- apply(y.mat.2f, 2, t.test.pval.func)
      pvals.3f <- apply(y.mat.3f, 2, t.test.pval.func)
      pvals.4f <- apply(y.mat.4f, 2, t.test.pval.func)
      
      pvals.234 <- apply(rbind(y.mat.2f, y.mat.3f, y.mat.4f), 2, t.test.pval.func)
      pvals.134 <- apply(rbind(y.mat.1f, y.mat.3f, y.mat.4f), 2, t.test.pval.func)
      pvals.124 <- apply(rbind(y.mat.1f, y.mat.2f, y.mat.4f), 2, t.test.pval.func)
      pvals.123 <- apply(rbind(y.mat.1f, y.mat.2f, y.mat.3f), 2, t.test.pval.func)
      
      alpha_t <- 0.05/4
      
      true.rej.num.global4 <- 0
      true.rej.num.rep4 <- 0
      
      for (f in 1:K1) {
        cond1 <- (pvals.1f[f] < alpha_t & pvals.234[f] < alpha_t / max(1, sum(pvals.1f < alpha_t)))
        cond2 <- (pvals.2f[f] < alpha_t & pvals.134[f] < alpha_t / max(1, sum(pvals.2f < alpha_t)))
        cond3 <- (pvals.3f[f] < alpha_t & pvals.124[f] < alpha_t / max(1, sum(pvals.3f < alpha_t)))
        cond4 <- (pvals.4f[f] < alpha_t & pvals.123[f] < alpha_t / max(1, sum(pvals.4f < alpha_t)))
        
        if (cond1 | cond2 | cond3 | cond4) {
          true.rej.num.global4 <- true.rej.num.global4 + 1
        }
        if (cond1 & cond2 & cond3 & cond4) {
          true.rej.num.rep4 <- true.rej.num.rep4 + 1
        }
      }
      
      reject.vec.four.team.global[s] <- true.rej.num.global4
      reject.vec.four.team.rep[s] <- true.rej.num.rep4
    }
    
    # ---- Summarize results ----
    methods <- c("Two Teams Global", "Three Teams Global", "Four Teams Global",
                 "Two Teams Replicability", "Three Teams Replicability", "Four Teams Replicability")
    
    power <- c(mean(reject.vec.two.team.global / K1),
               mean(reject.vec.three.team.global / K1),
               mean(reject.vec.four.team.global / K1),
               mean(reject.vec.two.team.rep / K1),
               mean(reject.vec.three.team.rep / K1),
               mean(reject.vec.four.team.rep / K1))
    
    se <- c(sd(reject.vec.two.team.global / K1) / sqrt(no.sims),
            sd(reject.vec.three.team.global / K1) / sqrt(no.sims),
            sd(reject.vec.four.team.global / K1) / sqrt(no.sims),
            sd(reject.vec.two.team.rep / K1) / sqrt(no.sims),
            sd(reject.vec.three.team.rep / K1) / sqrt(no.sims),
            sd(reject.vec.four.team.rep / K1) / sqrt(no.sims))
    
    df_temp <- data.frame(K1 = K1, mu = mu, Method = methods,
                          Power = power, SE = se)
    results <- rbind(results, df_temp)
  }
}


# ---------------- Plots --------------------------------

# ---------------- Prepare plotting data ----------------
results_plot <- results %>%
  separate(Method, into = c("Teams", "Type"), sep = " (?=[^ ]+$)") %>%
  mutate(
    Teams = factor(Teams, levels = c("Two Teams", "Three Teams", "Four Teams")),
    Type  = factor(Type, levels = c("Global", "Replicability")),
    # Set K1 factor for desired row order
    K1 = factor(K1, levels = c(1, 3, 5, 10)),
    K1_label = paste0("K[1] == ", K1)  # plotmath expression for subscript
  )
results_plot <- results_plot %>%
  mutate(
    # Define the order explicitly
    K1_label = factor(K1_label, levels = c(
      "K[1] == 1",
      "K[1] == 3",
      "K[1] == 5",
      "K[1] == 10"
    ))
  )
# ---------------- Define colors and legend labels ----------------
team_colors <- c(
  "Two Teams"   = "red",
  "Three Teams" = "green",
  "Four Teams"  = "blue"
)

team_labels <- c(
  "Two Teams"   = "Two team cross-screening",
  "Three Teams" = "Three team cross-screening",
  "Four Teams"  = "Four team cross-screening"
)

type_labels <- c(
  "Global"       = "Global null findings",
  "Replicability" = "Replicable findings"
)

# ---------------- Compute y-axis limits ----------------
y_min <- max(0, min(results_plot$Power - 2*results_plot$SE))
y_max <- min(1, max(results_plot$Power + 2*results_plot$SE))

# ---------------- Plot ----------------
p <- ggplot(results_plot, aes(x = mu, y = Power, color = Teams, group = Teams)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(0, Power - 2*SE),
                    ymax = pmin(1, Power + 2*SE)), width = 0.01) +
  facet_grid(
    K1_label ~ Type,
    scales = "fixed",
    labeller = labeller(Type = type_labels, K1_label = label_parsed)
  ) +
  scale_color_manual(values = team_colors, labels = team_labels) +
  scale_y_continuous(limits = c(y_min - 0.02, y_max + 0.02),
                     expand = expansion(mult = c(0,0))) +
  labs(x = expression(mu ~ "(Effect size)"),
       y = "Power",
       color = "Method") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

# ---------------- Display / Save ----------------
print(p)
#ggsave("matrix_plot.pdf", p, width = 14, height = 18)

