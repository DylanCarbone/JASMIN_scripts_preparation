# Trend utility functions
#'@export
expit <- function(x){ exp(x)/(1+exp(x)) }

#'@export
pcfunc <- function(pred){(tail(pred,1)-head(pred,1))/head(pred,1)*100}

#'@export
pcfunc2 <- function(pred1,pred2){(pred2-pred1)/pred1*100}

#'@export
sigfunc <- function(x){y <- ""; if(x <= 0.05) y <- "*";
  if(x <= 0.01) y <- "**";
  if(x <= 0.001) y <- "***";
  return(y)}
