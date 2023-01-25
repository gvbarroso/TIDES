# Date created: 23/05/2020
# Author: Gustavo V. Barroso
# This script combines posterior distributions of parameters inferred by TIDES
# from different simulated scenarios and 10 replicates each

library(reshape2)
library(cowplot)
library(tidyverse)
library(scales)
library(rethinking)

setSessionTimeLimit(cpu = Inf, elapsed = Inf) # some of the plots in the script can take a few seconds to generate

#setwd(getwd())

##################################
#
# posterior distributions of s for all reps
#
##################################

nreps <- 10
npost <- 18e+3 # number of posterior samples after rejection algorithm (nsim * tol)

# recyclable
post_dists_s <- as.data.frame(matrix(nrow = npost, ncol = 10))
names(post_dists_s) <- c(paste("rep_", seq(from = 1, to = 10), sep = ""))

# s = 0.01, h = 0.0
start_index <- 1
for(i in start_index:(start_index + nreps - 1)){ post_dists_s[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_s.txt", sep = ""))$V1 }

molten_s <- melt(post_dists_s)
dens_plot_s_1_10 <- ggplot(molten_s, aes(x = value, color = variable)) 
dens_plot_s_1_10 <- dens_plot_s_1_10 + geom_density(size = 1) + theme_bw() 
dens_plot_s_1_10 <- dens_plot_s_1_10 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                     "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_s_1_10 <- dens_plot_s_1_10 + geom_vline(aes(xintercept = 1e-2), linetype = "dashed", color = "black", size = 1)
dens_plot_s_1_10 <- dens_plot_s_1_10 + labs(title = NULL, x = NULL, y = "Density") 
dens_plot_s_1_10 <- dens_plot_s_1_10 + scale_x_log10(limits = c(1e-6, 0.1), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) # limits ~ range of prior
dens_plot_s_1_10 <- dens_plot_s_1_10 + theme(text = element_text(size = 14), 
                                             legend.position = "none",
                                             legend.title = element_blank(),
                                             axis.title.x = element_text(size = 18),
                                             axis.title.y = element_text(size = 16),
                                             axis.text.y = element_blank(),
                                             axis.ticks.y = element_blank())

s_median_1_10 <- apply(post_dists_s, 2, median)
s_median_1_10$model <- "strong/recessive"
s_hpdi_1_10 <- apply(post_dists_s, 2, HPDI, prob = 0.95)

# s = 0.01, h = 0.5
start_index <- 11
for(i in start_index:(start_index + nreps - 1)){ post_dists_s[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_s.txt", sep = ""))$V1 }

molten_s <- melt(post_dists_s)
dens_plot_s_11_20 <- ggplot(molten_s, aes(x = value, color = variable)) 
dens_plot_s_11_20 <- dens_plot_s_11_20 + geom_density(size = 1) + theme_bw() 
dens_plot_s_11_20 <- dens_plot_s_11_20 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                       "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_s_11_20 <- dens_plot_s_11_20 + geom_vline(aes(xintercept = 1e-2), linetype = "dashed", color = "black", size = 1)
dens_plot_s_11_20 <- dens_plot_s_11_20 + labs(title = NULL, x = " ", y = " ") 
dens_plot_s_11_20 <- dens_plot_s_11_20 + scale_x_log10(limits = c(1e-6, 0.1), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) # limits ~ range of prior
dens_plot_s_11_20 <- dens_plot_s_11_20 + theme(text = element_text(size = 14), 
                                               legend.position = "none",
                                               legend.title = element_blank(),
                                               axis.title.x = element_text(size = 18),
                                               axis.title.y = element_text(size = 16),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())

s_median_11_20 <- apply(post_dists_s, 2, median)
s_median_11_20$model <- "strong/additive"
s_hpdi_11_20 <- apply(post_dists_s, 2, HPDI, prob = 0.95)

# s = 0.001, h = 0.0
start_index <- 21
for(i in start_index:(start_index + nreps - 1)){ post_dists_s[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_s.txt", sep = ""))$V1 }

molten_s <- melt(post_dists_s)
dens_plot_s_21_30 <- ggplot(molten_s, aes(x = value, color = variable)) 
dens_plot_s_21_30 <- dens_plot_s_21_30 + geom_density(size = 1) + theme_bw() 
dens_plot_s_21_30 <- dens_plot_s_21_30 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                       "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_s_21_30 <- dens_plot_s_21_30 + geom_vline(aes(xintercept = 1e-3), linetype = "dashed", color = "black", size = 1)
dens_plot_s_21_30 <- dens_plot_s_21_30 + labs(title = NULL, x = " ", y = "Density") 
dens_plot_s_21_30 <- dens_plot_s_21_30 + scale_x_log10(limits = c(1e-6, 0.1), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) # limits ~ range of prior
dens_plot_s_21_30 <- dens_plot_s_21_30 + theme(text = element_text(size = 14), 
                                               legend.position = "none",
                                               legend.title = element_blank(),
                                               axis.title.x = element_text(size = 18),
                                               axis.title.y = element_text(size = 16),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())

s_median_21_30 <- apply(post_dists_s, 2, median)
s_median_21_30$model <- "moderate/recessive"
s_hpdi_21_30 <- apply(post_dists_s, 2, HPDI, prob = 0.95)

# s = 0.001, h = 0.5
start_index <- 31
for(i in start_index:(start_index + nreps - 1)){ post_dists_s[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_s.txt", sep = ""))$V1 }

molten_s <- melt(post_dists_s)
dens_plot_s_31_40 <- ggplot(molten_s, aes(x = value, color = variable)) 
dens_plot_s_31_40 <- dens_plot_s_31_40 + geom_density(size = 1) + theme_bw() 
dens_plot_s_31_40 <- dens_plot_s_31_40 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                       "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_s_31_40 <- dens_plot_s_31_40 + geom_vline(aes(xintercept = 1e-3), linetype = "dashed", color = "black", size = 1)
dens_plot_s_31_40 <- dens_plot_s_31_40 + labs(title = NULL, x = " ", y = " ") 
dens_plot_s_31_40 <- dens_plot_s_31_40 + scale_x_log10(limits = c(1e-6, 0.1), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) # limits ~ range of prior
dens_plot_s_31_40 <- dens_plot_s_31_40 + theme(text = element_text(size = 14), 
                                               legend.position = "none",
                                               legend.title = element_blank(),
                                               axis.title.x = element_text(size = 18),
                                               axis.title.y = element_text(size = 16),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())

s_median_31_40 <- apply(post_dists_s, 2, median)
s_median_31_40$model <- "moderate/additive"
s_hpdi_31_40 <- apply(post_dists_s, 2, HPDI, prob = 0.95)

# s = 0.0001, h = 0.0
start_index <- 41
for(i in start_index:(start_index + nreps - 1)){ post_dists_s[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_s.txt", sep = ""))$V1 }

molten_s <- melt(post_dists_s)
dens_plot_s_41_50 <- ggplot(molten_s, aes(x = value, color = variable)) 
dens_plot_s_41_50 <- dens_plot_s_41_50 + geom_density(size = 1) + theme_bw() 
dens_plot_s_41_50 <- dens_plot_s_41_50 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                       "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_s_41_50 <- dens_plot_s_41_50 + geom_vline(aes(xintercept = 1e-4), linetype = "dashed", color = "black", size = 1)
dens_plot_s_41_50 <- dens_plot_s_41_50 + labs(title = NULL, x = "|s|", y = "Density") 
dens_plot_s_41_50 <- dens_plot_s_41_50 + scale_x_log10(limits = c(1e-6, 0.1), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) # limits ~ range of prior
dens_plot_s_41_50 <- dens_plot_s_41_50 + theme(text = element_text(size = 14), 
                                               legend.position = "none",
                                               legend.title = element_blank(),
                                               axis.title.x = element_text(size = 18),
                                               axis.title.y = element_text(size = 16),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())

s_median_41_50 <- apply(post_dists_s, 2, median)
s_median_41_50$model <- "weak/recessive"
s_hpdi_41_50 <- apply(post_dists_s, 2, HPDI, prob = 0.95)

# s = 0.0001, h = 0.5
start_index <- 51
for(i in start_index:(start_index + nreps - 1)){ post_dists_s[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_s.txt", sep = ""))$V1 }

molten_s <- melt(post_dists_s)
dens_plot_s_51_60 <- ggplot(molten_s, aes(x = value, color = variable)) 
dens_plot_s_51_60 <- dens_plot_s_51_60 + geom_density(size = 1) + theme_bw() 
dens_plot_s_51_60 <- dens_plot_s_51_60 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                       "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_s_51_60 <- dens_plot_s_51_60 + geom_vline(aes(xintercept = 1e-4), linetype = "dashed", color = "black", size = 1)
dens_plot_s_51_60 <- dens_plot_s_51_60 + labs(title = NULL, x = "|s|", y = " ") 
dens_plot_s_51_60 <- dens_plot_s_51_60 + scale_x_log10(limits = c(1e-6, 0.1), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) # limits ~ range of prior
dens_plot_s_51_60 <- dens_plot_s_51_60 + theme(text = element_text(size = 14), 
                                               legend.position = "none",
                                               legend.title = element_blank(),
                                               axis.title.x = element_text(size = 18),
                                               axis.title.y = element_text(size = 16),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())

s_median_51_60 <- apply(post_dists_s, 2, median)
s_median_51_60$model <- "weak/additive"
s_hpdi_51_60 <- apply(post_dists_s, 2, HPDI, prob = 0.95)

# put dens plots together
post_dists <- plot_grid(dens_plot_s_1_10, dens_plot_s_11_20, 
                        dens_plot_s_21_30, dens_plot_s_31_40,
                        dens_plot_s_41_50, dens_plot_s_51_60,
                        ncol = 2, labels = "AUTO", label_size = 20)

cowplot::save_plot("post_dists_s.pdf", plot = post_dists, device = "pdf", dpi = 500, base_width = 12, base_height = 8)

# dotplot
s_hpdi <- rbind.data.frame(s_hpdi_1_10, s_hpdi_11_20, s_hpdi_21_30, s_hpdi_31_40, s_hpdi_41_50, s_hpdi_51_60)
s_hpdi$bound_95 <- rep(c("lower", "upper"), 6)
s_hpdi$model <- c(rep("strong/recessive", 2), rep("strong/additive", 2), rep("moderate/recessive", 2), rep("moderate/additive", 2), rep("weak/recessive", 2), rep("weak/additive", 2))
s_hpdi_lower <- s_hpdi[s_hpdi$bound_95 == "lower",]
s_hpdi_upper <- s_hpdi[s_hpdi$bound_95 == "upper",]
molten_ci_lower <- melt(s_hpdi_lower)
molten_ci_upper <- melt(s_hpdi_upper)

s_medians <- rbind.data.frame(s_median_1_10, s_median_11_20, s_median_21_30, s_median_31_40, s_median_41_50, s_median_51_60)
s_medians$model <- as.factor(s_medians$model)

molten_medians <- melt(s_medians)
molten_medians$lower <- molten_ci_lower$value
molten_medians$upper <- molten_ci_upper$value
molten_medians$model <- with(molten_medians, reorder(model, value))

dotplot_s <- ggplot(molten_medians, aes(x = model, y = value, color = model, fill = model))
dotplot_s <- dotplot_s + scale_y_log10(limits = c(1e-6, 0.1), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) 
dotplot_s <- dotplot_s + geom_pointrange(position = position_jitter(0.3),
                                         aes(y = value, group = variable, ymin = lower, ymax = upper, color = model),
                                         alpha = 0.5, fatten = 5, size = 0.5) 
dotplot_s <- dotplot_s + theme_bw()
dotplot_s <- dotplot_s + labs(title = NULL, x = "Model", y = "Inferred |s|") 
dotplot_s <- dotplot_s + scale_color_manual(values = c("#000000", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#e7298a"))
dotplot_s <- dotplot_s + geom_segment(aes(x = 0.7, xend = 1.3, y = 1e-4, yend = 1e-4, color = "black"))
dotplot_s <- dotplot_s + geom_segment(aes(x = 1.7, xend = 2.3, y = 1e-4, yend = 1e-4, color = "black"))
dotplot_s <- dotplot_s + geom_segment(aes(x = 2.7, xend = 3.3, y = 1e-3, yend = 1e-3, color = "black"))
dotplot_s <- dotplot_s + geom_segment(aes(x = 3.7, xend = 4.3, y = 1e-3, yend = 1e-3, color = "black"))
dotplot_s <- dotplot_s + geom_segment(aes(x = 4.7, xend = 5.3, y = 1e-2, yend = 1e-2, color = "black"))
dotplot_s <- dotplot_s + geom_segment(aes(x = 5.7, xend = 6.3, y = 1e-2, yend = 1e-2, color = "black"))
dotplot_s <- dotplot_s + theme(legend.position = "none",
                               legend.title = element_blank(),
                               axis.title.x = element_text(size = 18),
                               axis.title.y = element_text(size = 16),
                               axis.text.x = element_text(angle = 45, vjust = 0.6, size = 10),
                               axis.text.y = element_text(size = 10)) 
dotplot_s
ggsave("dotplot_s.pdf", dotplot_s, device = "pdf", dpi = 500, width = 5, height = 5)


# change s from 1e-3 to 1e-4 in last generation
start_index <- 61
for(i in start_index:(start_index + nreps - 1)){ post_dists_s[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_s.txt", sep = ""))$V1 }

molten_s <- melt(post_dists_s)
dens_plot_s_61_70 <- ggplot(molten_s, aes(x = value, color = variable)) 
dens_plot_s_61_70 <- dens_plot_s_61_70 + geom_density(size = 1) + theme_bw() 
dens_plot_s_61_70 <- dens_plot_s_61_70 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                       "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_s_61_70 <- dens_plot_s_61_70 + geom_vline(aes(xintercept = 1e-4), linetype = "dashed", color = "black", size = 1)
dens_plot_s_61_70 <- dens_plot_s_61_70 + geom_vline(aes(xintercept = 1e-3), linetype = "dotted", color = "black", size = 1)
dens_plot_s_61_70 <- dens_plot_s_61_70 + labs(title = NULL, x = "|s|", y = "Density") 
dens_plot_s_61_70 <- dens_plot_s_61_70 + scale_x_log10(limits = c(1e-6, 0.1), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) # limits ~ range of prior
dens_plot_s_61_70 <- dens_plot_s_61_70 + theme(text = element_text(size = 14), 
                                                 legend.position = "none",
                                                 legend.title = element_blank(),
                                                 axis.title.x = element_text(size = 18),
                                                 axis.title.y = element_text(size = 16),
                                                 axis.text.y = element_blank(),
                                                 axis.ticks.y = element_blank())

ggsave("post_dist_s_change.pdf", dens_plot_s_61_70, device = "pdf", dpi = 500)


# neutral sims (s = 0), prior [-0.1, 0.1], must be plotted with x-axis in linear scale
start_index <- 91
for(i in start_index:(start_index + nreps - 1)){ post_dists_s[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_s.txt", sep = ""))$V1 }

molten_s <- melt(post_dists_s)
dens_plot_s_91_100 <- ggplot(molten_s, aes(x = value, color = variable)) 
dens_plot_s_91_100 <- dens_plot_s_91_100 + stat_density(geom = "line", size = 1) + theme_bw() 
dens_plot_s_91_100 <- dens_plot_s_91_100 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                         "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_s_91_100 <- dens_plot_s_91_100 + geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black", size = 1)
dens_plot_s_91_100 <- dens_plot_s_91_100 + labs(title = NULL, x = "s", y = "Density") 
dens_plot_s_91_100 <- dens_plot_s_91_100 + xlim(-0.1, 0.1) 
dens_plot_s_91_100 <- dens_plot_s_91_100 + theme(text = element_text(size = 14), 
                                               legend.position = "none",
                                               legend.title = element_blank(),
                                               axis.title.x = element_text(size = 18),
                                               axis.title.y = element_text(size = 16),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())

ggsave("post_dist_neutral.pdf", dens_plot_s_91_100, device = "pdf", dpi = 500)

##################################
#
# posterior distributions of h for all reps
#
##################################

nreps <- 10
npost <- 18e+3

# recyclable
post_dists_h <- as.data.frame(matrix(nrow = npost, ncol = 10))
names(post_dists_h) <- c(paste("rep_", seq(from = 1, to = 10), sep = ""))

# s = 0.01, h = 0.0
start_index <- 1
for(i in start_index:(start_index + nreps - 1)){ post_dists_h[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_h.txt", sep = ""))$V1 }

molten_h <- melt(post_dists_h)
dens_plot_h_1_10 <- ggplot(molten_h, aes(x = value, color = variable)) 
dens_plot_h_1_10 <- dens_plot_h_1_10 + geom_density(size = 1) + theme_bw() 
dens_plot_h_1_10 <- dens_plot_h_1_10 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                     "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_h_1_10 <- dens_plot_h_1_10 + geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black", size = 1)
dens_plot_h_1_10 <- dens_plot_h_1_10 + labs(title = NULL, x = NULL, y = "Density") 
dens_plot_h_1_10 <- dens_plot_h_1_10 + scale_x_continuous(limits = c(0, 0.6), breaks = pretty_breaks())
dens_plot_h_1_10 <- dens_plot_h_1_10 + theme(text = element_text(size = 14), 
                                             legend.position = "none",
                                             legend.title = element_blank(),
                                             axis.title.x = element_text(size = 18),
                                             axis.title.y = element_text(size = 16),
                                             axis.text.y = element_blank(),
                                             axis.ticks.y = element_blank())


# s = 0.01, h = 0.5
start_index <- 11
for(i in start_index:(start_index + nreps - 1)){ post_dists_h[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_h.txt", sep = ""))$V1 }

molten_h <- melt(post_dists_h)
molten_h[which(molten_h$value > 0.6), 2] <- 0.6 # because of posterior shrinkage using abc package, some values are slightly about prior limits

dens_plot_h_11_20 <- ggplot(molten_h, aes(x = value, color = variable)) 
dens_plot_h_11_20 <- dens_plot_h_11_20 + geom_density(size = 1) + theme_bw() 
dens_plot_h_11_20 <- dens_plot_h_11_20 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                       "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_h_11_20 <- dens_plot_h_11_20 + geom_vline(aes(xintercept = 0.5), linetype = "dashed", color = "black", size = 1)
dens_plot_h_11_20 <- dens_plot_h_11_20 + labs(title = NULL, x = NULL, y = " ") 
dens_plot_h_11_20 <- dens_plot_h_11_20 + scale_x_continuous(limits = c(0, 0.6), breaks = pretty_breaks())
dens_plot_h_11_20 <- dens_plot_h_11_20 + theme(text = element_text(size = 14), 
                                               legend.position = "none",
                                               legend.title = element_blank(),
                                               axis.title.x = element_text(size = 18),
                                               axis.title.y = element_text(size = 16),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())


# s = 0.001, h = 0.0
start_index <- 21
for(i in start_index:(start_index + nreps - 1)){ post_dists_h[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_h.txt", sep = ""))$V1 }

molten_h <- melt(post_dists_h)
dens_plot_h_21_30 <- ggplot(molten_h, aes(x = value, color = variable)) 
dens_plot_h_21_30 <- dens_plot_h_21_30 + geom_density(size = 1) + theme_bw() 
dens_plot_h_21_30 <- dens_plot_h_21_30 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                       "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_h_21_30 <- dens_plot_h_21_30 + geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black", size = 1)
dens_plot_h_21_30 <- dens_plot_h_21_30 + labs(title = NULL, x = " ", y = "Density") 
dens_plot_h_21_30 <- dens_plot_h_21_30 + scale_x_continuous(limits = c(0, 0.6), breaks = pretty_breaks())
dens_plot_h_21_30 <- dens_plot_h_21_30 + theme(text = element_text(size = 14), 
                                               legend.position = "none",
                                               legend.title = element_blank(),
                                               axis.title.x = element_text(size = 18),
                                               axis.title.y = element_text(size = 16),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())

# s = 0.001, h = 0.5
start_index <- 31
for(i in start_index:(start_index + nreps - 1)){ post_dists_h[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_h.txt", sep = ""))$V1 }

molten_h <- melt(post_dists_h)
molten_h[which(molten_h$value > 0.6), 2] <- 0.6 # because of posterior shrinkage using abc package, some values are slightly about prior limits

dens_plot_h_31_40 <- ggplot(molten_h, aes(x = value, color = variable)) 
dens_plot_h_31_40 <- dens_plot_h_31_40 + geom_density(size = 1) + theme_bw() 
dens_plot_h_31_40 <- dens_plot_h_31_40 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                       "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_h_31_40 <- dens_plot_h_31_40 + geom_vline(aes(xintercept = 0.5), linetype = "dashed", color = "black", size = 1)
dens_plot_h_31_40 <- dens_plot_h_31_40 + labs(title = NULL, x = " ", y = " ") 
dens_plot_h_31_40 <- dens_plot_h_31_40 + scale_x_continuous(limits = c(0, 0.6), breaks = pretty_breaks())
dens_plot_h_31_40 <- dens_plot_h_31_40 + theme(text = element_text(size = 14), 
                                               legend.position = "none",
                                               legend.title = element_blank(),
                                               axis.title.x = element_text(size = 18),
                                               axis.title.y = element_text(size = 16),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())

# s = 0.0001, h = 0
start_index <- 41
for(i in start_index:(start_index + nreps - 1)){ post_dists_h[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_h.txt", sep = ""))$V1 }

molten_h <- melt(post_dists_h)
molten_h[which(molten_h$value > 0.6), 2] <- 0.6 # because of posterior shrinkage using abc package, some values are slightly about prior limits

dens_plot_h_41_50 <- ggplot(molten_h, aes(x = value, color = variable)) 
dens_plot_h_41_50 <- dens_plot_h_41_50 + geom_density(size = 1) + theme_bw() 
dens_plot_h_41_50 <- dens_plot_h_41_50 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                       "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_h_41_50 <- dens_plot_h_41_50 + geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black", size = 1)
dens_plot_h_41_50 <- dens_plot_h_41_50 + labs(title = NULL, x = "h", y = "Density") 
dens_plot_h_41_50 <- dens_plot_h_41_50 + scale_x_continuous(limits = c(0, 0.6), breaks = pretty_breaks())
dens_plot_h_41_50 <- dens_plot_h_41_50 + theme(text = element_text(size = 14), 
                                               legend.position = "none",
                                               legend.title = element_blank(),
                                               axis.title.x = element_text(size = 18),
                                               axis.title.y = element_text(size = 16),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())

# s = 0.0001, h = 0.5
start_index <- 51
for(i in start_index:(start_index + nreps - 1)){ post_dists_h[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_h.txt", sep = ""))$V1 }

molten_h <- melt(post_dists_h)
molten_h[which(molten_h$value > 0.6), 2] <- 0.6 # because of posterior shrinkage using abc package, some values are slightly about prior limits

dens_plot_h_51_60 <- ggplot(molten_h, aes(x = value, color = variable)) 
dens_plot_h_51_60 <- dens_plot_h_51_60 + geom_density(size = 1) + theme_bw() 
dens_plot_h_51_60 <- dens_plot_h_51_60 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                       "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_h_51_60 <- dens_plot_h_51_60 + geom_vline(aes(xintercept = 0.5), linetype = "dashed", color = "black", size = 1)
dens_plot_h_51_60 <- dens_plot_h_51_60 + labs(title = NULL, x = "h", y = " ") 
dens_plot_h_51_60 <- dens_plot_h_51_60 + scale_x_continuous(limits = c(0, 0.6), breaks = pretty_breaks())
dens_plot_h_51_60 <- dens_plot_h_51_60 + theme(text = element_text(size = 14), 
                                               legend.position = "none",
                                               legend.title = element_blank(),
                                               axis.title.x = element_text(size = 18),
                                               axis.title.y = element_text(size = 16),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())

# put them together
post_dists <- plot_grid(dens_plot_h_1_10, dens_plot_h_11_20, 
                        dens_plot_h_21_30, dens_plot_h_31_40,
                        dens_plot_h_41_50, dens_plot_h_51_60,
                        ncol = 2, labels = "AUTO", label_size = 20)

cowplot::save_plot("post_dists_h.pdf", plot = post_dists, device = "pdf", dpi = 500, base_width = 12, base_height = 8)


# change s from 1e-3 to 1e-4 in last generation
start_index <- 61
for(i in start_index:(start_index + nreps - 1)){ post_dists_h[,(i - start_index + 1)] <- read.table(paste("script_", i, "/posterior_dist_h.txt", sep = ""))$V1 }

molten_h <- melt(post_dists_h)
dens_plot_h_61_70 <- ggplot(molten_h, aes(x = value, color = variable)) 
dens_plot_h_61_70 <- dens_plot_h_61_70 + geom_density(size = 1) + theme_bw() 
dens_plot_h_61_70 <- dens_plot_h_61_70 + scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                                                       "#e7298a", "#a65628", "#f781bf", "#fb8072", "#762a83"))
dens_plot_h_61_70 <- dens_plot_h_61_70 + geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black", size = 1)
dens_plot_h_61_70 <- dens_plot_h_61_70 + labs(title = NULL, x = "h", y = "Density") 
dens_plot_h_61_70 <- dens_plot_h_61_70 + theme(text = element_text(size = 14), 
                                               legend.position = "none",
                                               legend.title = element_blank(),
                                               axis.title.x = element_text(size = 18),
                                               axis.title.y = element_text(size = 16),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())

ggsave("post_dist_h_change.pdf", dens_plot_h_61_70, device = "pdf", dpi = 500)


##################################
#
# model posterior probs 
#
##################################

#library(ggtern)

#mpp_simplex <- ggtern(data = mpp, aes(x = additive, y = neutral, z = recessive, shape = model, color = s)) + geom_mask() #+ theme_bw()
#mpp_simplex <- mpp_simplex + scale_color_manual(values = c("#7fc97f", "#beaed4")) + theme_showarrows()
#mpp_simplex <- mpp_simplex + stat_density_tern(aes(alpha = ..level.., fill = s), 
#                                               geom = 'polygon', 
#                                               bins = 10,
#                                               color = "grey",
#                                               base = "identity")
#mpp_simplex <- mpp_simplex + geom_point(size = 3, alpha = 1)
#mpp_simplex
#ggsave("mpp_simplex.png", mpp_simplex, device = "png", dpi = 500)


##################################
#
# sum stats
#
##################################

# number of simulations in second step
nsims <- 1.8e+7

# model indicators
models <- c(rep("additive", nsims), rep("recessive", nsims), rep("neutral", nsims))

start_indices <- seq(from = 1, to = 71, by = 10)
true_model <- c(rep(c("recessive", "additive"), 3), "recessive")

pb <- txtProgressBar(min = 0, max = length(start_indices), style = 3)
for(i in 1:length(start_indices)) {
  
  setTxtProgressBar(pb, i)
  
  plot_list <- list(length = nreps)
  start_index <- start_indices[i]
  
  for(j in start_index:(start_index + nreps - 1)) {
  
    obs_sum_stats <- read.table(paste("script_", j, "/full_model/obs_sum_stats.txt", sep = ""), header = T)
  
    script <- paste("script_", j, sep = "")
  
    sim_sum_stats_additive <- read.table(gzfile(paste(script, "/additive_model/sum_stats_dist.gz", sep = "")), header = T)
    sim_sum_stats_recessive <- read.table(gzfile(paste(script, "/recessive_model/sum_stats_dist.gz", sep = "")), header = T)
    sim_sum_stats_neutral <- read.table(gzfile(paste(script, "/neutral_model/sum_stats_dist.gz", sep = "")), header = T)
  
    sim_sum_stats_models <- rbind.data.frame(sim_sum_stats_additive, sim_sum_stats_recessive, sim_sum_stats_neutral)
  
    sum_stats_2d <- cbind.data.frame(sim_sum_stats_models, models)
    sum_stats_2d_ds <- sample_n(sum_stats_2d, 6e+4) # downsample so the plot isn't too big
  
    obs <- data.frame(obs_sum_stats$diff_het, obs_sum_stats$diff_homo, true_model[i])
    names(obs) <- names(sum_stats_2d_ds)
    sum_stats_2d_ds <- rbind.data.frame(sum_stats_2d_ds, obs)
  
    ss2d <- ggplot(data = sum_stats_2d_ds, aes(x = diff_het, y = diff_homo, color = models))
    ss2d <- ss2d + geom_point(shape = 16, size = 2,  alpha = 0.15) + theme_bw()
    ss2d <- ss2d + labs(title = NULL, x = paste(expression(delta),"_het", sep = ""), y = paste(expression(delta),"_homo", sep = ""))
    ss2d <- ss2d + scale_color_manual(values = c("#7fc97f", "#beaed4", "#fdc086"))
    ss2d <- ss2d + geom_point(x = obs_sum_stats$diff_het, y = obs_sum_stats$diff_homo, col = 'red', shape = 8, size = 6)
    ss2d <- ss2d + guides(colour = guide_legend(override.aes = list(alpha = 1))) 
    ss2d <- ss2d + theme(#legend.position = "none",
                         legend.title = element_blank(),
                         axis.title.x = element_text(size = 18),
                         axis.title.y = element_text(size = 18),
                         axis.ticks.y = element_blank(),
                         axis.ticks.x = element_blank())
    
    plot_list[[j - start_index + 1]] <- ss2d
  }
  
  ss2d_reps <- plot_grid(plotlist = plot_list, ncol = 2, labels = "AUTO", label_size = 20)
  plot_name <- paste("ss2d_", as.character(start_index), "-", as.character(start_index + nreps - 1), ".pdf", sep = "")
  cowplot::save_plot(plot_name, ss2d_reps, device = "pdf", dpi = 500, base_height = 32, base_width = 16, limitsize = F)
}
close(pb)

######################
#
# tight linkage
#
######################

setwd("~/Data/TIDES/paper_1/tight_linkage/")

start_index <- 1
nreps <- 10
npost <- 2e+3

# s
post_dists_s_10X <- as.data.frame(matrix(nrow = npost, ncol = 10))
names(post_dists_s_10X) <- c(paste("rep_", seq(from = 1, to = 10), sep = ""))

# focal SNP with 10X stronger selection as linked SNPs
for(i in 1:nreps){
  post_dists_s_10X[,i] <- read.table(paste("10X_bgs/rep_", i, "/params_post_dist.tsv", sep = ""), header = T)$loclinear_s
}

molten_s_10X <- melt(post_dists_s_10X)

s_median_10X <- apply(post_dists_s_10X, 2, median)
s_hpdi_10X <- apply(post_dists_s_10X, 2, HPDI, prob = 0.95)

post_dists_s_1X <- as.data.frame(matrix(nrow = npost, ncol = 10))
names(post_dists_s_1X) <- c(paste("rep_", seq(from = 1, to = 10), sep = ""))

for(i in 1:nreps){
  post_dists_s_1X[,i] <- read.table(paste("1X_bgs/rep_", i, "/params_post_dist.tsv", sep = ""), header = T)$neuralnet_s
}

molten_s_1X <- melt(post_dists_s_1X)

s_median_1X <- apply(post_dists_s_1X, 2, median)
s_hpdi_1X <- apply(post_dists_s_1X, 2, HPDI, prob = 0.95)

s_hpdi <- rbind.data.frame(s_hpdi_1X, s_hpdi_10X)
s_hpdi$bound_95 <- c("upper", "lower", "upper", "lower")
s_hpdi$model <- as.factor(c("1X", "1X", "10X", "10X"))
s_hpdi_lower <- s_hpdi[s_hpdi$bound_95 == "lower",]
s_hpdi_upper <- s_hpdi[s_hpdi$bound_95 == "upper",]
molten_ci_lower <- melt(s_hpdi_lower)
molten_ci_upper <- melt(s_hpdi_upper)

s_medians <- rbind.data.frame(s_median_1X, s_median_10X)
names(s_medians) <- paste("rep_", 1:nreps, sep = "")
s_medians$model <- c(1, 10)

molten_medians <- melt(s_medians, id.vars = "model")
molten_medians$lower <- molten_ci_lower$value
molten_medians$upper <- molten_ci_upper$value

dotplot_s <- ggplot(molten_medians, aes(x = as.factor(model), y = abs(value)))
dotplot_s <- dotplot_s + scale_y_log10(limits = c(1e-3, 1), breaks = c(1e-3, 1e-2, 1e-1, 1)) 
dotplot_s <- dotplot_s + geom_pointrange(position = position_jitter(0.3),
                                         aes(y = abs(value), group = variable, ymin = abs(lower), ymax = abs(upper)),
                                         alpha = 0.5, fatten = 5, size = 1) 
dotplot_s <- dotplot_s + theme_bw()
dotplot_s <- dotplot_s + labs(title = NULL, x = "Model", y = "Inferred |s|") 
dotplot_s <- dotplot_s + geom_segment(aes(x = 0.7, xend = 1.3, y = 1e-1, yend = 1e-1))
dotplot_s <- dotplot_s + geom_segment(aes(x = 1.7, xend = 2.3, y = 1e-1, yend = 1e-1))
dotplot_s <- dotplot_s + theme(legend.position = "none",
                               legend.title = element_blank(),
                               axis.title.x = element_text(size = 18),
                               axis.title.y = element_text(size = 16),
                               axis.text.x = element_text(vjust = 0.6, size = 10),
                               axis.text.y = element_text(size = 10)) 
dotplot_s
ggsave("dotplot_s.pdf", dotplot_s, device = "pdf", dpi = 500, width = 5, height = 5)

# h
post_dists_h_10X <- as.data.frame(matrix(nrow = npost, ncol = 10))
names(post_dists_h_10X) <- c(paste("rep_", seq(from = 1, to = 10), sep = ""))

# focal SNP with 10X stronger selection as linked SNPs
for(i in 1:nreps){
  post_dists_h_10X[,i] <- read.table(paste("10X_bgs/rep_", i, "/params_post_dist.tsv", sep = ""), header = T)$loclinear_h
}

molten_h_10X <- melt(post_dists_h_10X)

h_median_10X <- apply(post_dists_h_10X, 2, median)
h_hpdi_10X <- apply(post_dists_h_10X, 2, HPDI, prob = 0.95)

post_dists_h_1X <- as.data.frame(matrix(nrow = npost, ncol = 10))
names(post_dists_h_1X) <- c(paste("rep_", seq(from = 1, to = 10), sep = ""))

for(i in 1:nreps){
  post_dists_h_1X[,i] <- read.table(paste("1X_bgs/rep_", i, "/params_post_dist.tsv", sep = ""), header = T)$neuralnet_h
}

molten_h_1X <- melt(post_dists_h_1X)

h_median_1X <- apply(post_dists_h_1X, 2, median)
h_hpdi_1X <- apply(post_dists_h_1X, 2, HPDI, prob = 0.95)

h_hpdi <- rbind.data.frame(h_hpdi_1X, h_hpdi_10X)
h_hpdi$bound_95 <- c("upper", "lower", "upper", "lower")
h_hpdi$model <- as.factor(c("1X", "1X", "10X", "10X"))
h_hpdi_lower <- h_hpdi[h_hpdi$bound_95 == "lower",]
h_hpdi_upper <- h_hpdi[h_hpdi$bound_95 == "upper",]
molten_ci_lower <- melt(h_hpdi_lower)
molten_ci_upper <- melt(h_hpdi_upper)

h_medians <- rbind.data.frame(h_median_1X, h_median_10X)
names(h_medians) <- paste("rep_", 1:nreps, sep = "")
h_medians$model <- c(1, 10)

molten_medians <- melt(h_medians, id.vars = "model")
molten_medians$lower <- molten_ci_lower$value
molten_medians$upper <- molten_ci_upper$value

dotplot_h <- ggplot(molten_medians, aes(x = as.factor(model), y = value))
dotplot_h <- dotplot_h + scale_y_continuous(limits = c(-0.6, 0.6), breaks = pretty_breaks()) 
dotplot_h <- dotplot_h + geom_pointrange(position = position_jitter(0.3),
                                         aes(y = value, group = variable, ymin = lower, ymax = upper,),
                                         alpha = 0.5, fatten = 5, size = 1) 
dotplot_h <- dotplot_h + theme_bw()
dotplot_h <- dotplot_h + labs(title = NULL, x = "Model", y = "Inferred h") 
dotplot_h <- dotplot_h + geom_segment(aes(x = 0.7, xend = 1.3, y = 0, yend = 0))
dotplot_h <- dotplot_h + geom_segment(aes(x = 1.7, xend = 2.3, y = 0, yend = 0))
dotplot_h <- dotplot_h + theme(legend.position = "none",
                               legend.title = element_blank(),
                               axis.title.x = element_text(size = 18),
                               axis.title.y = element_text(size = 16),
                               axis.text.x = element_text(vjust = 0.6, size = 10),
                               axis.text.y = element_text(size = 10)) 
dotplot_h
ggsave("dotplot_h.pdf", dotplot_h, device = "pdf", dpi = 500, width = 5, height = 5)

post_dists <- plot_grid(dotplot_s, dotplot_h, ncol = 1, labels = "AUTO", label_size = 20)

cowplot::save_plot("post_dists_tight.pdf", plot = post_dists, device = "pdf", dpi = 500, base_width = 8, base_height = 12)


