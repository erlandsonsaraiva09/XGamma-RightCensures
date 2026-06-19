###############################################################################
########## APPLICATION 1 ######################################################
###############################################################################

# Run Fit.functionsXGamma.R

     set.seed(19)

      library(ggfortify)
      library(survminer)
      library(survival)

#################################
##### DATA SET ##################
#################################

          ta <- c(1760, 260, 129, 3639, 4331, 4252, 687, 3, 2070, 709, 3466, 616, 3896, 3830, 454, 2650, 3915, 2333, 3754, 1270, 3704, 383, 3578, 2902, 3518, 3485, 2119, 2502, 3425, 3403, 715, 3198, 3110, 3209, 268, 25, 3014, 460, 2762, 1306, 2053, 3006, 2861, 1227, 2264, 841, 917, 2765, 2738, 2757, 2639, 736, 630, 2464, 2428, 1443, 654, 2355, 2278, 843, 2344, 2171, 2133, 2220, 1322, 594, 1960, 1927, 1832, 1941, 99, 1714, 151, 1697, 1692, 214, 1624, 1566, 1528, 1520, 487, 1481, 1410, 3, 1259, 1205, 1180, 572, 1120, 1103, 1065, 498, 991, 991, 994, 898, 969, 895, 893, 701, 810, 742, 758)
           t <- ta/1000
           d <- c(1,1,1,1,0,0,1,0,1,1,1,1,0,0,1,1,0,1,0,1,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,1,0,1,1,1,1,0,0,1,1,1,1,0,0,0,0,1,1,0,0,1,1,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,rep(0,6),1,0,0,0,1,rep(0,11))

           n <- length(t)
           r <- sum(d)

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


