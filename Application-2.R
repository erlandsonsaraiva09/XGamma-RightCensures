###############################################################################
########## APPLICATION 2 ######################################################
###############################################################################

# Run Fit.functionsXGamma.R

    set.seed(19)

     library(KMsurv)

        data(larynx)

###############################################################
##### SURVIVAL DATA ###########################################
###############################################################

          t <- larynx$time/12

          d <- larynx$delta

###############################################################
##### SAMPLE INFORMATION ######################################
###############################################################

           n <- length(t)

           r <- sum(d)

        cens <- n - r

       pcens <- 100*cens/n

          cat("Sample size =", n, "\n")
          cat("Failures =", r, "\n")
          cat("Censored =", cens, "\n")
          cat("Percentage censored =", round(pcens,2), "%\n")

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
       P.EST
       
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
                  
##### Points estimates and RMSE.KM
                  
          res <- rbind(Estimate = fit$Estimates$Estimate, RMSE.KM = RMSE.KM)
                       colnames(res) <- fit$Estimates$Method
                       res <- round(res, 4)
          res
                  
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

###########################
##### PREDICTIVE PERFORMACE
###########################

####################################
##### PREDICTIVE LEAVE-ONE-OUT #####
####################################

     library(parallel)

 system.time({
        CV <- CVPredict.XGamma(t = t, d = d, NCORES = 6)})

CV$Performance


