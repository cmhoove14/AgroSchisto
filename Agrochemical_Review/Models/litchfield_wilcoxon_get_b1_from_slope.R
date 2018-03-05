library(rootSolve)

get_b1 <- function(slp){
  uniroot.all(f = function(b1){2*slp / 10^(-qnorm(.16)/b1) - (10^(qnorm(.84)/b1) / 10^(-qnorm(.16)/b1)) - 1}, interval = c(0, slp*100), n = 1e6)
}
