###############################################################################
########## SIMULATED DATA #####################################################
###############################################################################

# Run Fit.functionsXGamma.R

     set.seed(19)

      library(survival)

###############################################################################
##### SIMULATION STUDY 1 ######################################################
###############################################################################

  theta.real <- 0.6
           n <- 100

           t <- numeric(n)
           d <- rep(1, n)

for(i in 1:n){
           u <- runif(1)
           if(u < theta.real/(1 + theta.real))
        t[i] <- rexp(1, theta.real)
    else
        t[i] <- rgamma(1, 3, theta.real)
           C <- runif(1,1,4.9)
           if(t[i] > C)
             {
        t[i] <- C
        d[i] <- 0
             }
             }

           r <- sum(d)

          cat("Number of failures =",r,"\n")
          cat("Number of censored observations =",n-r,"\n")

################################################################################
#################### FIT OF XGAMMA MODEL #######################################
################################################################################

         fit <- fit.XGamma(t, d)

###############################################################
##### Point and Interval estimates and calibration values #####
###############################################################

# Calibration values
          
        C.VA <-  fit$Calibration
        C.VA
         
# Point Estimates

       P.EST <- fit$Estimates

# Absolute Errors

          AE <- round(abs(fit$Estimates$Estimate - theta.real), 4)

##### Points estimates, Absolute errors and RMSE.KM

          res <- rbind(Estimate = fit$Estimates$Estimate, AE = AE)
colnames(res) <- fit$Estimates$Method
          res <- round(res, 4)
          res
 
 ##### RMSE_KM #####

           KM <- survfit(Surv(t,d) ~ 1)
         t.KM <- KM$time
         S.KM <- KM$surv
    theta.hat <- fit$Estimates$Estimate
          
      RMSE.KM <- sapply(theta.hat,
                        function(th){S.fit <- SXG(t.KM, th)
                                     sqrt(mean((S.fit - S.KM)^2))}
                       )
          
                  round(RMSE.KM, 4)          
          
# Interval estimates

     INT.EST <- fit$Intervals
     INT.EST

#############################
##### CONVERGENCE CHECK #####
#############################

DIAG <- diag.mcmc.xg(fit)

DIAG$Summary

########
#### IMH 
########

DIAG$TracePlot.IMH

DIAG$IMH$GelmanPlot()

DIAG$IMH$ErgodicPlot

DIAG$IMH$ACFPlot

########
#### RWM
########

DIAG$TracePlot.RWM

DIAG$RWM$ErgodicPlot

DIAG$RWM$GelmanPlot()

DIAG$RWM$ACFPlot
