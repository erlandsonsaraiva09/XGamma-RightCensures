###############################################################################
###############################################################################
##### FUNCTIONS FOR ESTIMATION OF THE PARAMETER OF THE XGAMMA #################
###############################################################################
###############################################################################

######################################
##### Survival Function ##############
######################################

       SXG <- function(tt, theta){((1 + theta + theta*tt + theta^2*tt^2/2) * exp(-theta*tt))/(1 + theta)}


###############################################################################
##### MLE FOR THE XGAMMA DISTRIBUTION WITH RIGHT-CENSORED DATA ################
###############################################################################

MLE.XGamma <- function(t, d)
           {
  
###################################
##### Negative log-likelihood #####
###################################
  
    x.gama <- function(theta)
           {
     theta <- theta[1]
         if(theta <= 0)
     return(1e10)
    
    loglik <- sum(d*(2*log(theta) - log(1 + theta) + log(1 + theta*t^2/2) - theta*t)) + sum((1-d)*(log(1 + theta + theta*t + theta^2*t^2/2) - log(1 + theta) - theta*t))
    
     return(-loglik)
           }
  
###################################
##### Score function ##############
###################################
  
     score <- function(theta)
           {
         r <- sum(d)
         n <- length(t)
    
        A1 <- 2*r/theta
        A2 <- n/(1 + theta)
    
        A3 <- sum(d*t^2/(2 + theta*t^2))
        A4 <- sum((1-d)*(1 + t + theta*t^2)/(1 + theta + theta*t + theta^2*t^2/2))
        A5 <- sum(t)
           (A1 - A2 + A3 - A5 + A4)
           }
  
###################################
##### MLE #########################
###################################
  
    est.ml <- optim(par = 1/mean(t), fn = x.gama, method = "L-BFGS-B", lower = 1e-8, hessian = TRUE)
  
###################################
##### Results #####################
###################################
  
  THETA.ML <- est.ml$par
     SE.ML <- sqrt(1/as.numeric(est.ml$hessian))
   INFO.ML <- as.numeric(est.ml$hessian)
  
     IC.ML <- c(THETA.ML - 1.96*SE.ML, THETA.ML + 1.96*SE.ML)
  
###################################
##### Output ######################
###################################
  
        cat("\n")
        cat("MLE =", round(THETA.ML, 4), "\n")
        cat("SE  =", round(SE.ML, 4), "\n")
        cat("Observed Information =", round(INFO.ML, 4), "\n")
        cat("95% Wald CI = (", round(IC.ML[1], 4), ", ", round(IC.ML[2], 4), ")\n")
  
        cat("\n")
        cat("Score at MLE =", round(score(THETA.ML), 8), "\n")
  
###################################
##### Return ######################
###################################
  
     return(
       list(theta.ml = THETA.ML,
            se.ml = SE.ML,
            info.ml = INFO.ML,
            ci.ml = IC.ML,
            score.ml = score(THETA.ML),
            optim.output = est.ml
           )
           )
           }

###############################################################################
##### IMPORTANCE SAMPLING #####################################################
###############################################################################

  IS.XGamma <- function(t, d, gamma = 0.01, beta = 0.01, M = 5000)
            {
     library(coda)
  
          n <- length(t)
          r <- sum(d)
  
   pp.gamma <- gamma + 2*r
    pp.beta <- beta + sum(t)
  
###########################################################################
##### log A(theta) ########################################################
###########################################################################
  
 logA.theta <- function(theta)
            {
            -n*log(1 + theta) + sum(d*log(1 + theta*t^2/2)) + sum((1-d)*log(1 + theta + theta*t + theta^2*t^2/2))
            }
  
###########################################################################
##### Generate proposal sample ###########################################
###########################################################################
  
      theta <- rgamma(M, shape = pp.gamma, rate = pp.beta)
  
###########################################################################
##### Importance weights ##################################################
###########################################################################
  
         lw <- sapply(theta, logA.theta)
         lw <- lw - max(lw)
          w <- exp(lw)
          w <- w/sum(w)
     ESS.IS <- 1/sum(w^2)
    RESS.IS <- ESS.IS/length(w)
  
###########################################################################
##### Bayes estimator under SE loss ######################################
###########################################################################
  
   theta.SE <- sum(theta*w)
  
###########################################################################
##### Equal-tail credible interval #######################################
###########################################################################
  
        ord <- order(theta)
 theta.sort <- theta[ord]
     w.sort <- w[ord]
         Fw <- cumsum(w.sort)
         Li <- theta.sort[which(Fw >= 0.025)[1]]
         Ls <- theta.sort[which(Fw >= 0.975)[1]]
  
      IC.SE <- c(Li, Ls)  

###########################################################################
##### HPD credible interval ###############################################
###########################################################################
  
 theta.post <- sample(x = theta, size = 10000, replace = TRUE, prob = w)
  
        HPD <- HPDinterval(mcmc(theta.post), prob = 0.95)
  
     IC.HPD <- c(HPD[1,"lower"], HPD[1,"upper"])
  
###########################################################################
##### Return ##############################################################
###########################################################################
  
      return(
        list(sample = theta,
             weights = w,
             theta.se = theta.SE,
             ci.se = IC.SE,
             ci.hpd = IC.HPD,
             ESS = ESS.IS,
             RelativeESS = RESS.IS
            )
            )
            }

###############################################################################
##### INDEPENDENT METROPOLIS-HASTINGS #########################################
###############################################################################

 IMH.XGamma <- function(t, d, gamma = 0.01, beta = 0.01, L = 30000, B = 5000, J = 5)
            {
  
     library(coda)
  
###########################################################################
##### POSTERIOR HYPERPARAMETERS ###########################################
###########################################################################
  
          n <- length(t)
          r <- sum(d)
  
   pp.gamma <- gamma + 2*r
   pp.beta  <- beta + sum(t)
  
###########################################################################
##### log(A(theta)) #######################################################
###########################################################################
  
 logA.theta <- function(theta)
            {
          if(theta <= 0)
      return(-Inf)
            -n*log(1 + theta) + sum(d*log(1 + theta*t^2/2)) + sum((1-d)*log(1 + theta + theta*t + theta^2*t^2/2))
            }
  
###############################################################################
##### CHAIN 1 #################################################################
###############################################################################
  
       repeat{
      theta1 <- numeric(L)
   theta1[1] <- rgamma(1, shape = pp.gamma, rate  = pp.beta)
   n.accept1 <- 0
    
for(l in 2:L){
  theta.curr <- theta1[l-1]
  theta.cand <- rgamma(1, shape = pp.gamma, rate  = pp.beta)
      
   log.ratio <- logA.theta(theta.cand) - logA.theta(theta.curr)
      
       alpha <- min(1, exp(log.ratio))
           if(runif(1) < alpha)
             {
   theta1[l] <- theta.cand
   n.accept1 <- n.accept1 + 1
             }
         else
             {
   theta1[l] <- theta.curr
             }
             }
    
 theta1.keep <- theta1[seq(B+1, L, J)]
           if(sd(theta1.keep) > 1e-8)
      break
          cat("Chain 1 is constant. Regenerating...\n")
             }
  
   ACC.RATE1 <- n.accept1/(L-1)
  
###############################################################################
##### CHAIN 2 #################################################################
###############################################################################
  
       repeat{
      theta2 <- numeric(L)
   theta2[1] <- rgamma(1, shape = pp.gamma, rate  = pp.beta)
   n.accept2 <- 0
    
for(l in 2:L){
  theta.curr <- theta2[l-1]
  theta.cand <- rgamma(1, shape = pp.gamma, rate  = pp.beta)
      
   log.ratio <- logA.theta(theta.cand) - logA.theta(theta.curr)
       alpha <- min(1, exp(log.ratio))
           if(runif(1) < alpha)
             {
   theta2[l] <- theta.cand
   n.accept2 <- n.accept2 + 1
             }
         else
             {
   theta2[l] <- theta.curr
             }
             }
    
 theta2.keep <- theta2[seq(B+1, L, J)]
    
           if(sd(theta2.keep) > 1e-8)
        break
          cat("Chain 2 is constant. Regenerating...\n")
             }
  
   ACC.RATE2 <- n.accept2/(L-1)
  
###########################################################################
##### BURN-IN AND THINNING ################################################
###########################################################################
  
 theta1.keep <- theta1[seq(B+1, L, J)]
 theta2.keep <- theta2[seq(B+1, L, J)]
  
           N <- length(theta1.keep)
  
###########################################################################
##### CODA OBJECTS ########################################################
###########################################################################
  
      chain1 <- mcmc(theta1.keep)
      chain2 <- mcmc(theta2.keep)
  
      CHAINS <- mcmc.list(chain1, chain2)
  
###########################################################################
##### GELMAN-RUBIN ########################################################
###########################################################################
  
          GR <- gelman.diag(CHAINS, autoburnin = FALSE)
  
###########################################################################
##### ESS #################################################################
###########################################################################
  
        ESS1 <- as.numeric(effectiveSize(chain1))
        ESS2 <- as.numeric(effectiveSize(chain2))
  
    ESS.REL1 <- ESS1/N
    ESS.REL2 <- ESS2/N
  
###########################################################################
##### BAYES ESTIMATOR (SE LOSS) ###########################################
###########################################################################
  
    theta.SE <- mean(theta1.keep)
       IC.SE <- quantile(theta1.keep, probs = c(0.025,0.975))
  
###########################################################################
##### OUTPUT ##############################################################
###########################################################################
  
          cat("\n")
          cat("--------------------------------------------------\n")
          cat("IMH RESULTS\n")
          cat("--------------------------------------------------\n")
          cat("Acceptance Rate Chain 1 =", round(ACC.RATE1,4), "\n")
          cat("Acceptance Rate Chain 2 =", round(ACC.RATE2,4), "\n\n")
          cat("ESS Chain 1 =", round(ESS1,2), "\n")
          cat("ESS Chain 2 =", round(ESS2,2), "\n\n")
          cat("Relative ESS Chain 1 =", round(ESS.REL1,4), "\n")
          cat("Relative ESS Chain 2 =", round(ESS.REL2,4), "\n\n")
          cat("Posterior Mean =", round(theta.SE,4), "\n")
          cat("95% Credible Interval = (", round(IC.SE[1],4), ", ", round(IC.SE[2],4), ")\n")
          cat("Retained Sample Size =", N, "\n")
  
###########################################################################
##### RETURN ##############################################################
###########################################################################
  
       return(
         list(sample = theta1.keep,
              theta.se = theta.SE,
              ci.se = IC.SE,
              acceptance.rate = c(Chain1 = ACC.RATE1, Chain2 = ACC.RATE2),
              ESS = c(Chain1 = ESS1, Chain2 = ESS2),
              RelativeESS = c(Chain1 = ESS.REL1, Chain2 = ESS.REL2),
              GelmanRubin = GR,
              chains = CHAINS,
              chain1 = chain1,
              chain2 = chain2
             )
             )
             }

###############################################################################
##### ADAPTIVE RANDOM-WALK METROPOLIS #########################################
###############################################################################

    RWM.XGamma <- function(t, d, gamma = 0.01, beta = 0.01, L = 30000, B = 5000, J = 5, sigma0 = 0.10, adapt.block = 200)
               {
  
        library(coda)
  
###########################################################################
##### SAMPLE CHARACTERISTICS ##############################################
###########################################################################
  
             n <- length(t)
             r <- sum(d)
  
###########################################################################
##### LOG-POSTERIOR #######################################################
###########################################################################
  
      log.post <- function(theta)
               {
             if(theta <= 0)
         return(-Inf)
          LOGA <- -n*log(1 + theta) + sum(d*log(1 + theta*t^2/2)) + sum((1-d)*log(1 + theta + theta*t + theta^2*t^2/2))
          LOGP <- (gamma + 2*r - 1)*log(theta) - (beta + sum(t))*theta
          LOGA + LOGP
               }
  
###########################################################################
##### CHAIN 1 #############################################################
###########################################################################
  
        theta1 <- numeric(L)
     theta1[1] <- 1/mean(t)
         sigma <- sigma0
     n.accept1 <- 0
n.accept.block <- 0
  
  for(l in 2:L){
    theta.curr <- theta1[l-1]
    theta.cand <- abs(theta.curr + rnorm(1,0,sigma))
     log.ratio <- log.post(theta.cand) - log.post(theta.curr)
   
         alpha <- min(1,exp(log.ratio))
             if(runif(1) < alpha)
               {
     theta1[l] <- theta.cand
     n.accept1 <- n.accept1 + 1
n.accept.block <- n.accept.block + 1
               }
           else
               {
     theta1[l] <- theta.curr
               }
    
             if(l <= B && l %% adapt.block == 0)
               {
     acc.local <- n.accept.block/adapt.block
             if(acc.local < 0.20)
         sigma <- 0.9*sigma
             if(acc.local > 0.50)
         sigma <- 1.1*sigma
n.accept.block <- 0
               }
               }
  
     ACC.RATE1 <- n.accept1/(L-1)
  
###########################################################################
##### CHAIN 2 #############################################################
###########################################################################
  
        theta2 <- numeric(L)
     theta2[1] <- 2/mean(t)
         sigma <- sigma0
     n.accept2 <- 0
n.accept.block <- 0
  
  for(l in 2:L){
    theta.curr <- theta2[l-1]
    theta.cand <- abs(theta.curr + rnorm(1,0,sigma))
     log.ratio <- log.post(theta.cand) - log.post(theta.curr)
  
         alpha <- min(1,exp(log.ratio))
             if(runif(1) < alpha)
               {
     theta2[l] <- theta.cand
     n.accept2 <- n.accept2 + 1
n.accept.block <- n.accept.block + 1
               }
           else
               {
     theta2[l] <- theta.curr
               }
    
             if(l <= B && l %% adapt.block == 0)
               {
     acc.local <- n.accept.block/adapt.block
             if(acc.local < 0.20)
         sigma <- 0.9*sigma
             if(acc.local > 0.50)
         sigma <- 1.1*sigma
n.accept.block <- 0
               }
               }
  
     ACC.RATE2 <- n.accept2/(L-1)
  
###########################################################################
##### BURN-IN AND THINNING ################################################
###########################################################################
  
   theta1.keep <- theta1[seq(B+1,L,J)]
   theta2.keep <- theta2[seq(B+1,L,J)]
             N <- length(theta1.keep)
  
###########################################################################
##### CODA OBJECTS ########################################################
###########################################################################
  
        chain1 <- mcmc(theta1.keep)
        chain2 <- mcmc(theta2.keep)
        CHAINS <- mcmc.list(chain1,chain2)
  
###########################################################################
##### GELMAN-RUBIN ########################################################
###########################################################################
  
            GR <- gelman.diag(CHAINS, autoburnin = FALSE)
  
###########################################################################
##### EFFECTIVE SAMPLE SIZE ###############################################
###########################################################################
  
          ESS1 <- as.numeric(effectiveSize(chain1))
          ESS2 <- as.numeric(effectiveSize(chain2))
  
      ESS.REL1 <- ESS1/N
      ESS.REL2 <- ESS2/N
  
###########################################################################
##### BAYES ESTIMATOR (SE LOSS) ###########################################
###########################################################################
  
      theta.SE <- mean(theta1.keep)
         IC.SE <- quantile(theta1.keep, probs = c(0.025,0.975))
  
###########################################################################
##### RETURN ##############################################################
###########################################################################
  
         return(
           list(sample = theta1.keep,
                theta.se = theta.SE,
                ci.se = IC.SE,
                acceptance.rate = c(
                Chain1 = ACC.RATE1,
                Chain2 = ACC.RATE2),
                ESS = c(Chain1 = ESS1, Chain2 = ESS2),
                RelativeESS = c(Chain1 = ESS.REL1, Chain2 = ESS.REL2),
                GelmanRubin = GR,
                chains = CHAINS,
                chain1 = chain1,
                chain2 = chain2
               )
               )
               }

###############################################################################
##### CALIBRATION OF LINEX AND GE ESTIMATORS ##################################
###############################################################################

     AB.XGamma <- function(t, d, theta.IS, w.IS, theta.MH, grid.min = -10, grid.max = 10, by = 0.1)
               {
  
        library(survival)
  
###########################################################################
##### KAPLAN-MEIER ########################################################
###########################################################################
  
            KM <- survfit(Surv(t,d) ~ 1)
  
          t.KM <- KM$time
          S.KM <- KM$surv
  
###########################################################################
##### SURVIVAL FUNCTION ###################################################
###########################################################################
  
            Sy <- function(tt, theta)
               {
               ((1 + theta + theta*tt + theta^2*tt^2/2)*exp(-theta*tt))/(1 + theta)
               }
  
###########################################################################
##### RMSE ################################################################
###########################################################################
  
        RMSE.S <- function(theta.est)
               {
         S.est <- Sy(t.KM, theta.est)
           sqrt(mean((S.est - S.KM)^2))
               }
  
###########################################################################
##### GRID SEARCH #########################################################
###########################################################################
  
     AB.Search <- function(theta.sample, weights = NULL, grid.min, grid.max, by)
               {
        LAST.A <- NA
        LAST.B <- NA
         repeat{
          GRID <- seq(grid.min, grid.max, by)
      
###############################################################
##### LINEX ###################################################
###############################################################
      
      RMSE.LI  <- rep(NA, length(GRID))
      THETA.LI <- rep(NA, length(GRID))
        STOP.A <- FALSE
      
            for(i in seq_along(GRID)){
             a <- GRID[i]
             if(abs(a) > 1e-10)
               {
             if(is.null(weights))
               {
      theta.li <- -(1/a)*log(mean(exp(-a*theta.sample)))
               }
           else
               {
      theta.li <- -(1/a)*log(sum(weights*exp(-a*theta.sample)))
               }
          
             if(!is.finite(theta.li))
               {
        STOP.A <- TRUE
        a.star <- LAST.A
         break
               }
          
        LAST.A <- a
   THETA.LI[i] <- theta.li
   RMSE.LI[i]  <- RMSE.S(theta.li)
               }
               }
      
###############################################################
##### GE ######################################################
###############################################################
      
      RMSE.GE  <- rep(NA, length(GRID))
      THETA.GE <- rep(NA, length(GRID))
      
        STOP.B <- FALSE
      
            for(i in seq_along(GRID))
               {
             b <- GRID[i]
             if(abs(b) > 1e-10)
               {
             if(is.null(weights))
               {
   theta.ge <- (mean(theta.sample^(-b)))^(-1/b)
               }
           else
               {
      theta.ge <- (sum(weights*theta.sample^(-b)))^(-1/b)
               }
             if(!is.finite(theta.ge))
               {
        STOP.B <- TRUE
        b.star <- LAST.B
         break
               }
          
        LAST.B <- b
   THETA.GE[i] <- theta.ge
   RMSE.GE[i]  <- RMSE.S(theta.ge)
               }
               }
      
###############################################################
##### OVERFLOW STOP ###########################################
###############################################################
      
             if(STOP.A || STOP.B)
               {
        warning(paste("Numerical overflow detected.", "Returning a =", round(a.star,4), "and b =", round(b.star,4)))
          break
               }
      
###############################################################
##### BEST VALUES #############################################
###############################################################
      
            IA <- which.min(RMSE.LI)
            IB <- which.min(RMSE.GE)
      
        a.star <- GRID[IA]
        b.star <- GRID[IB]
      
          OK.A <- !(IA %in% c(1, length(GRID)))
          OK.B <- !(IB %in% c(1, length(GRID)))
      
###############################################################
##### SUCCESS #################################################
###############################################################
      
             if(OK.A && OK.B)
          break
      
###############################################################
##### EXPAND GRID #############################################
###############################################################
      
             if(!OK.A)
               {
             if(a.star <= 0)
               {
      grid.min <- a.star - 10
      grid.max <- a.star + 0.1
               }
           else
               {
      grid.min <- a.star - 0.1
      grid.max <- a.star + 10
               }
               }
      
             if(!OK.B)
               {
             if(b.star <= 0)
               {
      grid.min <- min(grid.min, b.star - 10)
      grid.max <- max(grid.max, b.star + 0.1)
               }
          else
               {
      grid.min <- min(grid.min, b.star - 0.1)
      grid.max <- max(grid.max, b.star + 10)
               }
               }
      
###############################################################
##### GRID LIMIT STOP #########################################
###############################################################
      
             if(grid.min <= -500 || grid.max >= 500)
               {
        warning(
          paste("Grid limit reached.", "Returning a =", round(LAST.A,4), "and b =", round(LAST.B,4)))
        a.star <- LAST.A
        b.star <- LAST.B
        break
               }
               }
    
###############################################################
##### CENTERED GRID ###########################################
###############################################################
    
       WIDTH.A <- 10
    
        GRID.A <- seq(a.star - WIDTH.A, a.star + WIDTH.A, by = by)
        RMSE.A <- rep(NA, length(GRID.A))
    
            for(i in seq_along(GRID.A))
               {
             a <- GRID.A[i]
             if(abs(a) > 1e-10)
               {
             if(is.null(weights))
               {
      theta.li <- -(1/a)*log(
           mean(exp(-a*theta.sample)))
               }
           else
               {
      theta.li <- -(1/a)*log(sum(weights*exp(-a*theta.sample)))
               }
        
     RMSE.A[i] <- RMSE.S(theta.li)
               }
               }
  
       WIDTH.B <- 10
    
        GRID.B <- seq(b.star - WIDTH.B, b.star + WIDTH.B, by = by)
        RMSE.B <- rep(NA, length(GRID.B))
    
            for(i in seq_along(GRID.B))
               {
             b <- GRID.B[i]
             if(abs(b) > 1e-10)
               {
             if(is.null(weights))
               {
      theta.ge <- (mean(theta.sample^(-b)))^(-1/b)
               }
           else
              {
     theta.ge <- (sum(weights*theta.sample^(-b)))^(-1/b)
              }
        
    RMSE.B[i] <- RMSE.S(theta.ge)
              }
              }
        
        return(
          list(a.star = a.star,
               b.star = b.star,
               theta.LI = THETA.LI[IA],
               theta.GE = THETA.GE[IB],
               RMSE.LI = RMSE.LI[IA],
               RMSE.GE = RMSE.GE[IB],
               GRID = GRID,
               RMSE.LI.ALL = RMSE.LI,
               RMSE.GE.ALL = RMSE.GE,
               THETA.LI.ALL = THETA.LI,
               THETA.GE.ALL = THETA.GE,
               GRID.A = GRID.A,
               RMSE.A = RMSE.A,
               GRID.B = GRID.B,
               RMSE.B = RMSE.B
              )
              )
              }
  
###########################################################################
##### IS ##################################################################
###########################################################################
  
        AB.IS <- AB.Search(theta.sample = theta.IS, weights = w.IS, grid.min = grid.min, grid.max = grid.max, by = by)
  
###########################################################################
##### MH ##################################################################
###########################################################################
  
        AB.MH <- AB.Search(theta.sample = theta.MH, grid.min = grid.min, grid.max = grid.max, by = by)
  
###########################################################################
##### OUTPUT ##############################################################
###########################################################################
  
           cat("\n")
           cat("a*_IS =", AB.IS$a.star,"\n") 
           cat("b*_IS =", AB.IS$b.star,"\n")
           cat("a*_MH =", AB.MH$a.star,"\n")
           cat("b*_MH =", AB.MH$b.star,"\n")

###########################################################################
##### RETURN ##############################################################
###########################################################################
  
        return(
          list(IS = AB.IS,
               MH = AB.MH,
               KM = KM
              )
              )
              }

###############################################################################
##### INTERVAL COMPARISON #####################################################
###############################################################################

Intervals.XGamma <- function(MLE, IS = NULL, IMH = NULL, RWM = NULL, prob = 0.95)
                 {
  
          library(coda)
  
           alpha <- 1 - prob
  
###########################################################################
##### MLE #################################################################
###########################################################################
  
            WALD <- MLE$ci.ml
  
###########################################################################
##### INITIALIZE TABLE ####################################################
###########################################################################
  
       INTERVALS <- data.frame(Method = character(), Interval = character(), Lower = numeric(), Upper = numeric(), Length = numeric())
  
###########################################################################
##### MLE #################################################################
###########################################################################
  
       INTERVALS <- rbind(INTERVALS, data.frame(Method = "MLE", Interval = "Wald", Lower = WALD[1], Upper = WALD[2], Length = WALD[2] - WALD[1]))
  
###########################################################################
##### IS ##################################################################
###########################################################################
  
               if(!is.null(IS))
                 {
    IS.EqualTail <- IS$ci.se
          IS.HPD <- IS$ci.hpd
    
#####################################################################
##### IS EQUAL-TAIL #################################################
#####################################################################
    
       INTERVALS <- rbind(INTERVALS, data.frame(Method = "IS", Interval = "Equal-tail", Lower = IS.EqualTail[1], Upper = IS.EqualTail[2], Length = IS.EqualTail[2] - IS.EqualTail[1]))
    
#####################################################################
##### IS HPD ########################################################
#####################################################################
    
       INTERVALS <- rbind( INTERVALS, data.frame(Method = "IS", Interval = "HPD", Lower = IS.HPD[1], Upper = IS.HPD[2], Length = IS.HPD[2] - IS.HPD[1]))
                 }
  
###########################################################################
##### IMH #################################################################
###########################################################################
  
               if(!is.null(IMH))
                 {
   IMH.EqualTail <- quantile(IMH$sample, probs = c(alpha/2, 1 - alpha/2))
         IMH.HPD <- HPDinterval(mcmc(IMH$sample), prob = prob)
    
       INTERVALS <- rbind(INTERVALS,
       data.frame(Method = "IMH", Interval = "Equal-tail", Lower = IMH.EqualTail[1], Upper = IMH.EqualTail[2],Length = IMH.EqualTail[2] - IMH.EqualTail[1]))
    
       INTERVALS <- rbind(INTERVALS, data.frame(Method = "IMH", Interval = "HPD", Lower = IMH.HPD[1,"lower"], Upper = IMH.HPD[1,"upper"],Length = IMH.HPD[1,"upper"] - IMH.HPD[1,"lower"]))
                 }
  
###########################################################################
##### RWM #################################################################
###########################################################################
  
               if(!is.null(RWM))
                 {
   RWM.EqualTail <- quantile(RWM$sample, probs = c(alpha/2, 1 - alpha/2))
    
         RWM.HPD <- HPDinterval(mcmc(RWM$sample), prob = prob)
    
       INTERVALS <- rbind(INTERVALS, data.frame(Method = "RWM", Interval = "Equal-tail", Lower = RWM.EqualTail[1], Upper = RWM.EqualTail[2], Length = RWM.EqualTail[2] - RWM.EqualTail[1]))
    
       INTERVALS <- rbind(INTERVALS, data.frame(Method = "RWM", Interval = "HPD", Lower = RWM.HPD[1,"lower"], Upper = RWM.HPD[1,"upper"], Length = RWM.HPD[1,"upper"] - RWM.HPD[1,"lower"]))
                 }
  
###########################################################################
##### ROUND ###############################################################
###########################################################################
  
 INTERVALS[,3:5] <- round(INTERVALS[,3:5], 4)
  
###########################################################################
##### RETURN ##############################################################
###########################################################################
  
           return(INTERVALS)
                 }

###############################################################################
##### CONVERGENCE DIAGNOSTICS #################################################
###############################################################################

Convergence.XGamma <- function(chain1, chain2 = NULL, lag.max = 100)
                   {
            library(coda)
            library(ggplot2)
  
###########################################################################
##### MCMC OBJECTS #########################################################
###########################################################################
  
             mcmc1 <- mcmc(chain1)
  
              ESS1 <- effectiveSize(mcmc1)
          ESS.REL1 <- ESS1/length(chain1)
  
           RESULTS <- list()
  
###########################################################################
##### TRACEPLOT ############################################################
###########################################################################
  
                 if(!is.null(chain2))
                   {
             mcmc2 <- mcmc(chain2)
            CHAINS <- mcmc.list(mcmc1, mcmc2)
    
                GR <- gelman.diag(CHAINS, autoburnin = FALSE)
    
              ESS2 <- effectiveSize(mcmc2)
          ESS.REL2 <- ESS2/length(chain2)
    
             TRACE <- traceplot(CHAINS, col = c(1,2), ylab = expression(theta))
                   } 
               else 
                  {
               GR <- NULL
             ESS2 <- NA
         ESS.REL2 <- NA
  
            TRACE <- traceplot(mcmc1, ylab = expression(theta))
                  }
  
###########################################################################
##### ERGODIC MEAN #########################################################
###########################################################################
  
             ERG1 <- cumsum(chain1)/(1:length(chain1))
  
           D.ERG1 <- data.frame(Iteration = 1:length(chain1), Mean = ERG1)
  
           P.ERG1 <- ggplot(D.ERG1, aes(x = Iteration, y = Mean)) +
                     geom_line() +
                     geom_hline(yintercept = mean(chain1), linetype = 2) +
                     theme_bw() +
                     ylab("Ergodic mean")
  
###########################################################################
##### ACF ##################################################################
###########################################################################
  
             ACF1 <- acf(chain1, lag.max = lag.max,plot = FALSE)
  
           D.ACF1 <- data.frame(Lag = as.numeric(ACF1$lag[-1]), ACF = as.numeric(ACF1$acf[-1]))
  
             crit <- 1.96/sqrt(length(chain1))
  
           P.ACF1 <- ggplot(D.ACF1, aes(x = Lag, y = ACF)) +
                     geom_segment(aes(xend = Lag, yend = 0)) +
                     geom_hline(yintercept = 0) +
                     geom_hline(yintercept = c(-crit,crit), linetype = 2) +
                     theme_bw()
  
###########################################################################
##### GELMAN-RUBIN PLOT ###################################################
###########################################################################
  
                if(!is.null(chain2))
                  {
             P.GR <- function()
                  {
       gelman.plot(CHAINS, autoburnin = FALSE)
                  }
                  } 
              else 
                  {
             P.GR <- NULL
                  }
  
###########################################################################
##### OUTPUT ###############################################################
###########################################################################
  
        RESULTS$ESS <- c(Chain1 = ESS1, Chain2 = ESS2)
  
RESULTS$RelativeESS <- c(Chain1 = ESS.REL1, Chain2 = ESS.REL2)
  
RESULTS$GelmanRubin <- GR
  
  RESULTS$Traceplot <- TRACE
  
RESULTS$ErgodicPlot <- P.ERG1
  
    RESULTS$ACFPlot <- P.ACF1
  
 RESULTS$GelmanPlot <- P.GR
  
              return(RESULTS)
                    }

################################################################################
################################################################################
################################################################################

         fit.XGamma <- function(t, d, gamma = 0.01, beta = 0.01, M = 20000, L = 120000, B = 20000, J = 5)
                    {
  
             library(coda)
  
###########################################################################
##### MLE #################################################################
###########################################################################
  
           invisible(
      capture.output(
                MLE <- MLE.XGamma(t = t, d = d)
                    )
                    )
  
###########################################################################
##### IS ##################################################################
###########################################################################
  
           invisible(
      capture.output(
                 IS <- IS.XGamma(t = t, d = d, gamma = gamma, beta = beta, M = M)
                    )
                    )
  
###########################################################################
##### IMH #################################################################
###########################################################################
  
           invisible(
      capture.output(
                IMH <- IMH.XGamma(t = t,d = d, gamma = gamma, beta = beta, L = L, B = B, J = J)
                    )
                    )
  

###########################################################################
##### RWM #################################################################
###########################################################################
  
          invisible(
     capture.output(
               RWM <- RWM.XGamma(t = t, d = d, gamma = gamma, beta = beta, L = L, B = B, J = J)
                   )
                   )
  
###########################################################################
##### CALIBRATION #########################################################
###########################################################################
  
          invisible(
     capture.output(
            AB.IMH <- AB.XGamma(t = t, d = d, theta.IS = IS$sample, w.IS = IS$weights, theta.MH = IMH$sample)
                   )
                   )
  
          invisible(
     capture.output(
            AB.RWM <- AB.XGamma(t = t, d = d, theta.IS = IS$sample, w.IS = IS$weights, theta.MH = RWM$sample)
                   )
                   )
  

###########################################################################
##### INTERVALS ###########################################################
###########################################################################
  
          invisible(
     capture.output(
         INTERVALS <- Intervals.XGamma(MLE = MLE, IS = IS, IMH = IMH, RWM = RWM)      
                   )
                   )
  
###########################################################################
##### ESTIMATES ###########################################################
###########################################################################
  
         ESTIMATES <- data.frame(Method = c("MLE","IS.SE","IS.LI","IS.GE", "IMH.SE","IMH.LI","IMH.GE", "RWM.SE","RWM.LI","RWM.GE"),
                                 Estimate = c(MLE$theta.ml,
                                              IS$theta.se,
                                              AB.IMH$IS$theta.LI,
                                              AB.IMH$IS$theta.GE,
                                              IMH$theta.se,
                                              AB.IMH$MH$theta.LI,
                                              AB.IMH$MH$theta.GE,
                                              RWM$theta.se,
                                              AB.RWM$MH$theta.LI,
                                              AB.RWM$MH$theta.GE)
                                )
  
###########################################################################
##### CALIBRATION #########################################################
###########################################################################
  
       CALIBRATION <- data.frame(Method = c("IS", "IMH", "RWM"),
                                 a.star = c(AB.IMH$IS$a.star, AB.IMH$MH$a.star, AB.RWM$MH$a.star),
                                 b.star = c(AB.IMH$IS$b.star, AB.IMH$MH$b.star, AB.RWM$MH$b.star)
                                )
  
###########################################################################
##### OUTPUT ##############################################################
###########################################################################
  
        OUT <- list(Estimates = ESTIMATES,
                    Intervals = INTERVALS,
                    Calibration = CALIBRATION,
                    AB.IMH = AB.IMH,
                    AB.RWM = AB.RWM,
                    MLE = MLE,
                    IS = IS,
                    IMH = IMH,
                    RWM = RWM
                    )
  
         class(OUT) <- "xgamma"
              return(OUT)
                    }

############################
###### MCMC DIAGNOSTIC #####
############################         

       diag.mcmc.xg <- function(object)
                    {
  
###########################################################################
##### IMH #################################################################
###########################################################################
  
           DIAG.IMH <- Convergence.XGamma(chain1 = as.numeric(object$IMH$chain1),chain2 = as.numeric(object$IMH$chain2))
  
###########################################################################
##### RWM #################################################################
###########################################################################
  
           DIAG.RWM <- Convergence.XGamma(chain1 = as.numeric(object$RWM$chain1), chain2 = as.numeric(object$RWM$chain2))
  
###########################################################################
##### SUMMARY TABLE #######################################################
###########################################################################
  
            SUMMARY <- data.frame(Method = c("IMH","RWM"),
                                  AcceptanceRate = c(mean(object$IMH$acceptance.rate), mean(object$RWM$acceptance.rate)),
                                  ESS = c(mean(as.numeric(DIAG.IMH$ESS)), mean(as.numeric(DIAG.RWM$ESS))),
                                  RelativeESS = c(mean(as.numeric(DIAG.IMH$RelativeESS)), mean(as.numeric(DIAG.RWM$RelativeESS))),
                                  Rhat = c(DIAG.IMH$GelmanRubin$psrf[1,1],
                                  DIAG.RWM$GelmanRubin$psrf[1,1])
                                 )
  
###########################################################################
##### ROUND ###############################################################
###########################################################################
  
SUMMARY$AcceptanceRate <- round(SUMMARY$AcceptanceRate, 4)
           SUMMARY$ESS <- round(SUMMARY$ESS, 2)
   SUMMARY$RelativeESS <- round(SUMMARY$RelativeESS, 4)
          SUMMARY$Rhat <- round(SUMMARY$Rhat, 4)
  
###############################################################################
##### TRACE PLOTS #############################################################
###############################################################################
  
                library(ggplot2)
  
             TRACE.IMH <- ggplot(
                           rbind(
                      data.frame(Iteration = seq_along(object$IMH$chain1),
                                 Theta = as.numeric(object$IMH$chain1),
                                 Chain = "Chain 1"),
                      data.frame(Iteration = seq_along(object$IMH$chain2),
                                 Theta = as.numeric(object$IMH$chain2),
                                 Chain = "Chain 2")
                                ),
                             aes(x = Iteration, y = Theta, colour = Chain)) +
                       geom_line(linewidth = 0.3) +
             scale_colour_manual(values = c("black", "red")) +
                      theme_gray() +
                           theme(legend.position = "top") +
                            labs(x = "Iteration",y = expression(theta))
  
             TRACE.RWM <- ggplot(
                           rbind(
                      data.frame(Iteration = seq_along(object$RWM$chain1),
                                 Theta = as.numeric(object$RWM$chain1),
                                 Chain = "Chain 1"),
                      data.frame(Iteration = seq_along(object$RWM$chain2),
                                 Theta = as.numeric(object$RWM$chain2),
                                 Chain = "Chain 2")
                                ),
                             aes(x = Iteration, y = Theta, colour = Chain)) +
                       geom_line(linewidth = 0.3) +
             scale_colour_manual(values = c("black", "red")) +
                      theme_gray() +
                           theme(legend.position = "top") +
                            labs(x = "Iteration", y = expression(theta))
  
###########################################################################
##### OUTPUT ##############################################################
###########################################################################
  
                   OUT <- list(Summary = SUMMARY,
                               IMH = DIAG.IMH,
                               RWM = DIAG.RWM,
                               TracePlot.IMH = TRACE.IMH,
                               TracePlot.RWM = TRACE.RWM
                              )
  
            class(OUT) <- "diag.xgamma"
                 return(OUT)
                       }

#########################################
########## GENERATOR OF XGAMMA ##########
#########################################

               rXGamma <- function(n, theta){
                     U <- runif(n)
                     X <- numeric(n)
                    I1 <- U <= theta/(1 + theta)
                 X[I1] <- rgamma(sum(I1), shape = 1, rate = theta)
                X[!I1] <- rgamma(sum(!I1), shape = 3, rate = theta)
                X
                       }

#########################################
########## PREDICTIVE APPROACH ##########
#########################################

        Predict.XGamma <- function(fit, N = 1000, alpha = 0.05, T.new = NULL)
                       {
  
##########################################################
##### MLE ################################################
##########################################################
  
                 T.MLE <- rXGamma(n = N, theta = fit$MLE$theta.ml)
  
##########################################################
##### IS #################################################
##########################################################
  
              theta.IS <- fit$IS$sample
                  w.IS <- fit$IS$weights
                  T.IS <- rXGamma(n = length(theta.IS), theta = theta.IS)
  
##########################################################
##### IMH ################################################
##########################################################
  
             theta.IMH <- fit$IMH$sample
                 T.IMH <- rXGamma(n = length(theta.IMH), theta = theta.IMH)
  
##########################################################
##### RWM ################################################
##########################################################
  
             theta.RWM <- fit$RWM$sample
                 T.RWM <- rXGamma(n = length(theta.RWM), theta = theta.RWM)
  
##########################################################
##### PREDICTIVE MEANS ###################################
##########################################################
  
              Mean.MLE <- mean(T.MLE)
               Mean.IS <- sum(w.IS*T.IS)
              Mean.IMH <- mean(T.IMH)
              Mean.RWM <- mean(T.RWM)
  
##########################################################
##### PREDICTIVE INTERVALS ###############################
##########################################################
  
                PI.MLE <- quantile(T.MLE, probs = c(alpha/2,1-alpha/2))
  
##########################################################
##### IS #################################################
##########################################################
  
                   ord <- order(T.IS)
             T.IS.sort <- T.IS[ord]
                w.sort <- w.IS[ord]
                    Fw <- cumsum(w.sort)
                 PI.IS <- c(T.IS.sort[which(Fw >= alpha/2)[1]], T.IS.sort[which(Fw >= 1-alpha/2)[1]])
  
##########################################################
##### IMH ################################################
##########################################################
  
                PI.IMH <- quantile(T.IMH, probs = c(alpha/2,1-alpha/2))
  
##########################################################
##### RWM ################################################
##########################################################
  
                PI.RWM <- quantile(T.RWM, probs = c(alpha/2,1-alpha/2))
  
##########################################################
##### SUMMARY ############################################
##########################################################
  
               SUMMARY <- data.frame(Method = c("MLE", "IS", "IMH", "RWM"),
                                     PredictiveMean = c(Mean.MLE, Mean.IS, Mean.IMH, Mean.RWM),
                                     Lower = c(PI.MLE[1], PI.IS[1], PI.IMH[1], PI.RWM[1]),
                                     Upper = c(PI.MLE[2], PI.IS[2], PI.IMH[2], PI.RWM[2]))
  
##########################################################
##### PREDICTIVE PERFORMANCE #############################
##########################################################
  
           PERFORMANCE <- NULL
                     if(!is.null(T.new))
                       {
           PERFORMANCE <- data.frame(Method = c("MLE", "IS", "IMH", "RWM"),
                                     PredictiveMean = c(Mean.MLE, Mean.IS, Mean.IMH, Mean.RWM),
                                     APE = c(abs(T.new - Mean.MLE), abs(T.new - Mean.IS), abs(T.new - Mean.IMH), abs(T.new - Mean.RWM)),
                                     SPE = c((T.new - Mean.MLE)^2, (T.new - Mean.IS)^2, (T.new - Mean.IMH)^2, (T.new - Mean.RWM)^2),
                                     Coverage = c(PI.MLE[1] <= T.new & T.new <= PI.MLE[2], PI.IS[1]  <= T.new & T.new <= PI.IS[2], PI.IMH[1] <= T.new & T.new <= PI.IMH[2], PI.RWM[1] <= T.new & T.new <= PI.RWM[2]),
                                     Length = c(PI.MLE[2] - PI.MLE[1], PI.IS[2]  - PI.IS[1], PI.IMH[2] - PI.IMH[1], PI.RWM[2] - PI.RWM[1]))
                       }
  
##########################################################
##### OUTPUT #############################################
##########################################################
  
                   OUT <- list(Summary = SUMMARY,
                               Performance = PERFORMANCE,
                               MLE = T.MLE,
                               IS = T.IS,
                               IMH = T.IMH,
                               RWM = T.RWM)
  
            class(OUT) <- "pred.xgamma"
           return(OUT)
                       }

###############################################################################
##### LEAVE-ONE-OUT PREDICTIVE PERFORMANCE (PARALLEL) #########################
###############################################################################
        
      CVPredict.XGamma <- function(t, d, NCORES = parallel::detectCores() - 1)
                       {
          
                library(parallel)
                library(pbmcapply)
          
###########################################################################
##### UNCENSORED OBSERVATIONS #############################################
###########################################################################
          
                     D <- which(d == 1)
                     r <- length(D)
          
                    cat("\n")
                    cat("Number of uncensored observations =", r, "\n")
                    cat("Using", NCORES, "cores\n")
                    cat("Starting leave-one-out cross-validation...\n\n")
          
###########################################################################
##### OBSERVED FAILURES ###################################################
###########################################################################
          
                 T.OBS <- t[D]
          
###########################################################################
##### PARALLEL LEAVE-ONE-OUT ##############################################
###########################################################################
          
                   OUT <- pbmcapply::pbmclapply(
                          X = seq_along(D),
                          FUN = function(j)
                       {
                     i <- D[j]
              
#######################################################################
##### REMOVE ONE FAILURE ##############################################
#######################################################################
              
               t.minus <- t[-i]
               d.minus <- d[-i]
              
#######################################################################
##### FIT MODEL #######################################################
#######################################################################
              
                   FIT <- fit.XGamma(t = t.minus, d = d.minus)
              
#######################################################################
##### PREDICTION ######################################################
#######################################################################
              
                    PR <- Predict.XGamma(fit = FIT)
              
#######################################################################
##### RETURN ##########################################################
#######################################################################
              
                   list(Mean  = PR$Summary$PredictiveMean,
                        Lower = PR$Summary$Lower,
                        Upper = PR$Summary$Upper)
                       },
              mc.cores = NCORES
                       )
          
                    cat("\nCross-validation completed.\n\n")
          
###########################################################################
##### ORGANIZE RESULTS ####################################################
###########################################################################
          
                  PRED <- do.call(rbind, lapply(OUT, function(x) x$Mean))
                 LOWER <- do.call(rbind, lapply(OUT, function(x) x$Lower))
                 UPPER <- do.call(rbind, lapply(OUT, function(x) x$Upper))
        colnames(PRED) <- c("MLE", "IS", "IMH", "RWM")
       colnames(LOWER) <- colnames(PRED)
       colnames(UPPER) <- colnames(PRED)
          
###########################################################################
##### RMSPE ###############################################################
###########################################################################
          
                 RMSPE <- apply(PRED, 2, function(x){sqrt(mean((T.OBS - x)^2, na.rm = TRUE))})
          
###########################################################################
##### COVERAGE PROBABILITY ################################################
###########################################################################
          
                    CP <- sapply(1:ncol(PRED), function(k){mean(T.OBS >= LOWER[,k] & T.OBS <= UPPER[,k], na.rm = TRUE)})
          
###########################################################################
##### AVERAGE LENGTH ######################################################
###########################################################################
          
                    AL <- sapply(1:ncol(PRED), function(k){mean(UPPER[,k] - LOWER[,k], na.rm = TRUE)})
          
###########################################################################
##### TABLE ###############################################################
###########################################################################
          
                   TAB <- data.frame(Method = colnames(PRED), RMSPE = round(RMSPE,4), Coverage = round(CP,4), AvgLength = round(AL,4))
          
###########################################################################
##### OUTPUT ##############################################################
###########################################################################
          
                 return(
                   list(Performance = TAB,
                        Observed = T.OBS,
                        Predicted = as.data.frame(PRED),
                        Lower = as.data.frame(LOWER),
                        Upper = as.data.frame(UPPER)
                       )
                       )
                       }
        