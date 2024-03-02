#####################################################################
#Capture-recapture history
#####################################################################

#Note: The script was modified for our case study.

#Capture-recapture history based on Kéry & Schaub (2011)

CRH_simulations <- function(n.occasions, N, p.value, simulations, phi.value, ...) {
  b.value <- runif(1, min = 0.05, max = 0.15)
  b <- c(b.value, rep((1 - b.value)/(n.occasions - 1), n.occasions - 1))
  p <- rep(p.value, n.occasions)
  P <- matrix(rep(p, n.occasions * N), ncol = n.occasions, nrow = N, byrow = T)
  data.sim <- list()
  for (i in 1:simulations) {
    phi <- rep(phi.value, n.occasions - 1)
    PHI <- matrix(rep(phi, (n.occasions - 1) * N), ncol = n.occasions - 1, nrow = N, byrow = T)
    simul.js <- function(PHI, P, b, N) {
      B <- rmultinom(1, N, b)
      n.occasions <- dim(PHI)[2] + 1
      CH.sur <- CH.p <- matrix(0, ncol = n.occasions, nrow = N)
      ent.occ <- numeric()
      for (t in 1:n.occasions) {
        ent.occ <- c(ent.occ, rep(t, B[t]))
      }
      for (i in 1:N) {
        CH.sur[i, ent.occ[i]] <- 1
        if (ent.occ[i] == n.occasions)
          next
        for (t in (ent.occ[i] + 1):n.occasions) {
          sur <- rbinom(1, 1, PHI[i, t - 1])
          ifelse(sur == 1, CH.sur[i, t] <- 1, break)
        }
      }
      for (i in 1:N) {
        CH.p[i, ] <- rbinom(n.occasions, 1, P[i, ])
      }
      CH <- CH.sur * CH.p
      cap.sum <- rowSums(CH)
      never <- which(cap.sum == 0)
      if (length(never > 0)) {
        CH <- CH[-never, ]
      }
      Nt <- colSums(CH.sur)
      return(list(CH = CH, B = B, N = Nt))
    }
    sim <- simul.js(PHI, P, b, N)
    CH <- sim$CH
    data.sim[[i]] <- data.frame(sim$CH)
  }
  return(data.sim)
}

## Bibliography
#
# Kéry, M., & Schaub, M. (2012). Bayesian population analysis using WinBUGS: A hierarchical
# perspective (1st ed). Academic Press.
