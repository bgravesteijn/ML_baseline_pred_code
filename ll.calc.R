ll.calc <- function(data, lev=NULL, model=NULL){
  obs <- as.numeric(data[,"obs"])
  pred <- as.numeric(data[,"pred"])
  pred <- plogis(pred)
  ll <- sum(obs * log(pred) + (1-obs)*log(1-pred))
  names(ll) <- "LogLik"
  return(ll)
}

