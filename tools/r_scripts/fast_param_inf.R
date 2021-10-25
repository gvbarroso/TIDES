#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(reshape2)
library(tidyverse)
library(abc)
library(rethinking)
library(scales)

setwd(getwd())

setSessionTimeLimit(cpu = Inf, elapsed = Inf) 

tolerance <- as.numeric(args[1])
true_s <- as.numeric(args[2])
true_h <- as.numeric(args[3])

obs_sum_stats <- read.table("obs_sum_stats.txt", header = T)
sim_sum_stats <- read.table(gzfile("sum_stats_dist.gz"), header = T)
sim_params <- read.table(gzfile("sim_params_dist.gz"), header = T)

nnet <- abc(target=obs_sum_stats, param=sim_params, sumstat=sim_sum_stats, tol=tolerance, method="neuralnet", numnet=6, transf="log")

num_accepted <- nrow(sim_sum_stats) * tolerance
prior_h <- runif(num_accepted, -0.1, 0.6)
prior_s <- 10^(-runif(num_accepted, 1, 6)) 

dat <- cbind.data.frame(nnet$unadj.values, prior_h, prior_s)
names(dat) <- c("nnet_h", "nnet_s", "prior_h", "prior_s")

write.table(dat$nnet_h, "posterior_dist_h.txt", quote = F, col.names = F, row.names = F, sep = "\n")
write.table(dat$nnet_s, "posterior_dist_s.txt", quote = F, col.names = F, row.names = F, sep = "\n")

s_hpdi <- HPDI(dat$nnet_s, prob = 0.95)
h_hpdi <- HPDI(dat$nnet_h, prob = 0.95)

dat_s <- dplyr::select(dat, ends_with("_s")) 
dat_s$sim <- 1:nrow(dat_s)
molten_vals <- melt(dat_s, id.vars = "sim")

dens_plot <- ggplot(molten_vals, aes(x = value, color = variable))
dens_plot <- dens_plot + stat_density(geom = "line", size = 1.125) + theme_bw() 
dens_plot <- dens_plot + scale_color_manual(values = c("red", "grey"))
dens_plot <- dens_plot + geom_vline(aes(xintercept = true_s), linetype = "dashed", color = "black", size = 0.85)
dens_plot <- dens_plot + labs(title = NULL, x = "s", y = "Posterior Density") 
dens_plot <- dens_plot + theme(text = element_text(size = 14), 
                               legend.position = "right",
                               legend.title = element_blank(),
                               axis.title.x = element_text(size = 18),
                               axis.title.y = element_text(size = 16))
dens_plot <- dens_plot + scale_x_log10(limits = c(1e-6, 0.1), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1))
ggsave("posterior_density_s.pdf", dens_plot, device = "pdf", dpi = 500)

dat_h <- dplyr::select(dat, ends_with("_h")) 
dat_h$sim <- 1:nrow(dat_h)
molten_vals <- melt(dat_h, id.vars = "sim")

dens_plot <- ggplot(molten_vals, aes(x = value, color = variable))
dens_plot <- dens_plot + stat_density(geom = "line", size = 1.125) + theme_bw() 
dens_plot <- dens_plot + scale_color_manual(values = c("red", "grey"))
dens_plot <- dens_plot + geom_vline(aes(xintercept = true_h), linetype = "dashed", color = "black", size = 1)
dens_plot <- dens_plot + labs(title = NULL, x = "h", y = "Posterior Density") 
dens_plot <- dens_plot + theme(text = element_text(size = 14), 
                               legend.position = "right",
                               legend.title = element_blank(),
                               axis.title.x = element_text(size = 18), 
                               axis.title.y = element_text(size = 16))
dens_plot <- dens_plot + scale_x_continuous(limits = c(-0.1, 0.6), breaks = pretty_breaks())
ggsave("posterior_density_h.pdf", dens_plot, device = "pdf", dpi = 500)

post_scatter <- ggplot(data = dat, aes(y = nnet_s, x = nnet_h))
post_scatter <- post_scatter + geom_point(shape = 16, size = 2,  alpha = 0.25) + theme_bw()
post_scatter <- post_scatter + labs(title = NULL, x = "h", y = "s")
post_scatter <- post_scatter + theme(axis.text.x = element_text(size = 14),
                                     axis.text.y = element_text(size = 14),
                                     axis.title.x = element_text(size = 24),
                                     axis.title.y = element_text(size = 24))
post_scatter <- post_scatter + geom_vline(aes(xintercept = true_h), color = "blue")
post_scatter <- post_scatter + geom_vline(aes(xintercept = h_hpdi[1]), color = "blue", linetype = "dashed", size = 0.5)
post_scatter <- post_scatter + geom_vline(aes(xintercept = h_hpdi[2]), color = "blue", linetype = "dashed", size = 0.5)
post_scatter <- post_scatter + geom_hline(aes(yintercept = true_s), color = "red")
post_scatter <- post_scatter + geom_hline(aes(yintercept = s_hpdi[1]), color = "red", linetype = "dashed", size = 0.5)
post_scatter <- post_scatter + geom_hline(aes(yintercept = s_hpdi[2]), color = "red", linetype = "dashed", size = 0.5)
post_scatter <- post_scatter + geom_density2d(colour = "green")
post_scatter <- post_scatter + scale_y_log10(breaks = pretty_breaks())
post_scatter <- post_scatter + scale_x_continuous(breaks = pretty_breaks())
ggsave("posterior_scatterplot.pdf", post_scatter, device = "pdf")

