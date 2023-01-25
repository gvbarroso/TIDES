# Created on 24/01/2023
# Author: Gustavo V. Barroso

# https://datawookie.dev/blog/2013/08/fitting-a-model-by-maximum-likelihood/

library(stats4)
library(iterpc)

counts <- read.table("genotype_counts.txt")

obs <- c(a, b, c, d, e, f, g)
obs <- c(100, 82, 43, 71, 34, 36, 24)
obs <- obs

LogLik <- function(s, h, a, b, c, d, e, f, g) {
  
  w1 <- (2 + h*s ) / 2
  w2 <- (4 + 2*h*s + s) / 4
  w3 <- (2 + h*s + s) / 2
  
  #mult <- multichoose(c(a, b, c, d, e, f, g), bigz=TRUE) 
  unifProb <- ((1/(2*w1))^a) * (((1 + h*s)/(2*w1))^b) * ((1/(4*w2))^c) * (((1 + h*s)/(2*w2))^d) * (((1 + s)/(4*w2))^e) * (((1 + h*s)/(2*w3))^f) * (((1 + s)/(2*w3))^g)

  cat(paste(unifProb, "\n"))
  
  #return(-log(mult * unifProb))
  return(-log(unifProb))
}

fit <- mle(LogLik, method="L-BFGS-B", start=list(s=-0.005, h=0.25), lower=list(s=-0.1, h=0), upper=(list(s=0, h=1)), 
           fixed=list(a=obs[1], b=obs[2], c=obs[3], d=obs[4], e=obs[5], f=obs[6], g=obs[7]))

summary(fit)
