# Created on 24/01/2023
# Author: Gustavo V. Barroso

# https://datawookie.dev/blog/2013/08/fitting-a-model-by-maximum-likelihood/

library(stats4)
library(iterpc)

file <- readLines("genotype_counts.txt")
obs <- cbind.data.frame(t(rbind.data.frame(strsplit(counts[9:15], " "))))[11]
row.names(obs) <- paste("c", 1:7, sep="")
colnames(obs) <- "count"
obs[,1] <- as.numeric(obs[,1])

LogLik <- function(s, h, c1, c2, c3, c4, c5, c6, c7) {
  
  w1 <- (2 + h*s ) / 2
  w2 <- (4 + 2*h*s + s) / 4
  w3 <- (2 + h*s + s) / 2
  
  mult <- multichoose(c(c1, c2, c3, c4, c5, c6, c7), bigz=TRUE) 
  
  # the following runs into underflow issues:
  #orderedProb <- ((1/(2*w1))^c1) * (((1 + h*s)/(2*w1))^c2) * ((1/(4*w2))^c3) * (((1 + h*s)/(2*w2))^c4) * (((1 + s)/(4*w2))^c5) * (((1 + h*s)/(2*w3))^c6) * (((1 + s)/(2*w3))^c7)
  #return(-log(mult * orderedProb)) # R cannot multiply a ridiculously large number by a ridiculously small one
  
  # instead we work directly in log space:
  logOrderedProb <- c1*log(1/(2*w1)) + c2*log((1 + h*s)/(2*w1)) + c3*log(1/(4*w2)) + c4*log((1 + h*s)/(2*w2)) + c5*log((1 + s)/(4*w2)) + c6*log((1 + h*s)/(2*w3)) + c7*log((1 + s)/(2*w3))
  
  
  return(-(log(mult) + logOrderedProb)) 
}

fit <- mle(LogLik, method="L-BFGS-B", start=list(s=-0.005, h=0.25), lower=list(s=-0.1, h=0), upper=(list(s=0, h=1)), 
           fixed=list(c1=obs[1,], c2=obs[2,], c3=obs[3,], c4=obs[4,], c5=obs[5,], c6=obs[6,], c7=obs[7,]))

summary(fit)
