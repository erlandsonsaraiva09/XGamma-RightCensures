#####################################################################
##### MONTE CARLOS FUNCTION #########################################
#####################################################################

MonteCarlo.XGamma <- function(theta = 0.6, n = 50, p.cens = 0.30, MC = 1000, K.pred = 30, gamma = 0.01, beta = 0.01)
                  {
           library(survival)
           library(progress)
           library(coda)
  
#####################################################################
##### TRUE SURVIVAL FUNCTION ########################################
#####################################################################
  
               Sy <- function(tt, theta){((1 + theta + theta*tt + theta^2*tt^2/2)*exp(-theta*tt))/(1 + theta)}
  
#####################################################################
##### ESTIMATORS ####################################################
#####################################################################
  
        EST.NAMES <- c("MLE", "IS.SE", "IS.LI", "IS.GE", "IMH.SE", "IMH.LI", "IMH.GE", "RWM.SE", "RWM.LI", "RWM.GE")
             NEST <- length(EST.NAMES)
  
#####################################################################
##### STORAGE #######################################################
#####################################################################
  
              EST <- matrix(NA, MC, NEST)
               AE <- matrix(NA, MC, NEST)
          RMSE.KM <- matrix(NA, MC, NEST)
  
#####################################################################
##### INTERVALS #####################################################
#####################################################################

        INT.NAMES <- c("MLE.Wald", "IS.EqualTail", "IS.HPD", "IMH.EqualTail", "IMH.HPD", "RWM.EqualTail", "RWM.HPD")
  
            COVER <- matrix(NA, MC,length(INT.NAMES))
           LENGTH <- matrix(NA, MC, length(INT.NAMES))
       PRED.NAMES <- c("MLE", "IS", "IMH", "RWM")
  
              APE <- matrix(NA, MC, 4)
              SPE <- matrix(NA, MC, 4)
          P.COVER <- matrix(NA, MC, 4)
         P.LENGTH <- matrix(NA, MC, 4)
  
#####################################################################
##### MCMC DIAGNOSTICS ##############################################
#####################################################################
  
          MCMC.SUM <- matrix(NA, MC, 8)
colnames(MCMC.SUM) <- c("IMH.ACC", "IMH.ESS", "IMH.RESS", "IMH.GR",
                        "RWM.ACC", "RWM.ESS", "RWM.RESS", "RWM.GR")
  
#####################################################################
##### PROGRESS BAR ##################################################
#####################################################################
  
                pb <- progress_bar$new(format = " Monte Carlo [:bar] :percent ETA: :eta",
                                       total = MC,
                                       clear = FALSE,
                                       width = 60)
  
###############################################################
##### CALIBRATE CENSORING #####################################
###############################################################
  
                 if(p.cens == 0){c.max <- Inf}
             else
                   {
             X.cal <- rXGamma(100000, theta)
             U.cal <- runif(100000)
            f.cens <- function(cmax)
                   {
             C.cal <- U.cal*cmax
               mean(C.cal < X.cal) - p.cens
                   }
    
             c.max <- tryCatch(uniroot(f.cens, interval = c(1e-6,100))$root,error = function(e) 100)
                   }
  
###############################################################
##### STORAGE #################################################
###############################################################
  
           SUCCESS <- 0
             NCENS <- rep(NA, MC)
  
     for(m in 1:MC){
    
###############################################################
##### GENERATE XGAMMA SAMPLE ##################################
###############################################################
    
                 T <- rXGamma(n, theta)
    
###############################################################
##### CENSORING ###############################################
###############################################################
    
                 if(is.infinite(c.max)){C <- rep(Inf,n)}
             else
                   {C <- runif(n,0,c.max)}
    
###############################################################
##### OBSERVED DATA ###########################################
###############################################################
    
             t.obs <- pmin(T,C)
                 d <- as.numeric(T <= C)
    
###############################################################
##### NUMBER OF CENSORED OBSERVATIONS #########################
###############################################################
    
          NCENS[m] <- sum(d == 0)
    
###############################################################
##### CHECK ALL CENSORED ######################################
###############################################################
    
                 if(sum(d) == 0)
                   {
                   pb$tick()
                   next
                   }
    
###############################################################
##### FIT MODEL ###############################################
###############################################################
    
               fit <- tryCatch(fit.XGamma(t = t.obs, d = d, gamma = gamma,beta = beta),
                      error = function(e){
                     cat("\nERROR IN fit.XGamma:\n")
                     print(e)
                     stop(e)})
    
                 if(is.null(fit))
                   {
                   pb$tick()
                   next
                   }
    
###############################################################
##### ESTIMATES ###############################################
###############################################################
    
           EST[m,] <- c(fit$MLE$theta.ml,
                        fit$IS$theta.se,
                        fit$AB.IMH$IS$theta.LI,
                        fit$AB.IMH$IS$theta.GE,
      
                        fit$IMH$theta.se,
                        fit$AB.IMH$MH$theta.LI,
                        fit$AB.IMH$MH$theta.GE,
       
                        fit$RWM$theta.se,
                        fit$AB.RWM$MH$theta.LI,
                        fit$AB.RWM$MH$theta.GE)
    
###############################################################
##### ABSOLUTE ERRORS #########################################
###############################################################
    
            AE[m,] <- abs(EST[m,] - theta)
    
###############################################################
##### KAPLAN-MEIER ############################################
###############################################################
    
                KM <- survfit(Surv(t.obs,d) ~ 1)
              t.KM <- KM$time
              S.KM <- KM$surv
                 if(length(t.KM)==0)
                   {
                   pb$tick()
                   next
                   }
    
###############################################################
##### RMSE ####################################################
###############################################################
    
   for(j in 1:NEST){
             S.est <- Sy(t.KM, EST[m,j])
      
      RMSE.KM[m,j] <- sqrt(mean((S.est - S.KM)^2))}
    
###############################################################
##### INTERVALS ###############################################
###############################################################
    
               INT <- fit$Intervals
    
for(k in 1:nrow(INT)){
               lower <- INT$Lower[k]
               upper <- INT$Upper[k]
      
          COVER[m,k] <- as.numeric(lower <= theta & theta <= upper)
      
         LENGTH[m,k] <- upper - lower
                     }
    
#####################################################################
##### MCMC DIAGNOSTICS ##############################################
#####################################################################
    
             IMH.ACC <- mean(fit$IMH$acceptance.rate)
             RWM.ACC <- mean(fit$RWM$acceptance.rate)
             IMH.ESS <- mean(fit$IMH$ESS)
             RWM.ESS <- mean(fit$RWM$ESS)
            IMH.RESS <- mean(fit$IMH$RelativeESS)
            RWM.RESS <- mean(fit$RWM$RelativeESS)
              IMH.GR <- fit$IMH$GelmanRubin$psrf[1,1]
              RWM.GR <- fit$RWM$GelmanRubin$psrf[1,1]
        MCMC.SUM[m,] <- c(IMH.ACC, IMH.ESS, IMH.RESS, IMH.GR,
                          RWM.ACC, RWM.ESS, RWM.RESS, RWM.GR)
    
###########
    
                pred <- Predict.XGamma(fit = fit, N = 20000)
    
               PMEAN <- pred$Summary$PredictiveMean
    
               LOWER <- pred$Summary$Lower
               UPPER <- pred$Summary$Upper
    
            T.future <- rXGamma(K.pred, theta)
    
             APE[m,] <- c(mean(abs(T.future - PMEAN[1])),
                          mean(abs(T.future - PMEAN[2])),
                          mean(abs(T.future - PMEAN[3])),
                          mean(abs(T.future - PMEAN[4])))
    
             SPE[m,] <- c(mean((T.future - PMEAN[1])^2),
                          mean((T.future - PMEAN[2])^2),
                          mean((T.future - PMEAN[3])^2),
                          mean((T.future - PMEAN[4])^2))
    
          P.COVER[m,] <- c(mean(as.numeric(T.future >= LOWER[1] & T.future <= UPPER[1])),
                           mean(as.numeric(T.future >= LOWER[2] & T.future <= UPPER[2])),
                           mean(as.numeric(T.future >= LOWER[3] & T.future <= UPPER[3])),
                           mean(as.numeric(T.future >= LOWER[4] & T.future <= UPPER[4])))
    
         P.LENGTH[m,] <- c(UPPER[1] - LOWER[1], UPPER[2] - LOWER[2], UPPER[3] - LOWER[3], UPPER[4] - LOWER[4])
    
              SUCCESS <- SUCCESS + 1
    
               pb$tick()
    
                    if(m %% 100 == 0)
                      {
                   cat(
               sprintf("PID %d | %d/%d (%.1f%%)\n",
            Sys.getpid(), m, MC, 100*m/MC))
         flush.console()
                      }
                      }
  
#####################################################################
##### TABLE 1 #######################################################
#####################################################################
  
                 TAB1 <- rbind(MeanEstimate = colMeans(EST, na.rm = TRUE),
                               MAE = colMeans(AE, na.rm = TRUE),
                               MeanRMSE.KM = colMeans(RMSE.KM, na.rm = TRUE))
       colnames(TAB1) <- EST.NAMES
  
#####################################################################
##### TABLE 2 #######################################################
#####################################################################
  
                 TAB2 <- rbind(Coverage = colMeans(COVER, na.rm = TRUE),
                               AvgLength = colMeans(LENGTH, na.rm = TRUE))
       colnames(TAB2) <- INT.NAMES
  
#####################################################################
##### TABLE 3 : PREDICTION ##########################################
#####################################################################
  
                 TAB3 <- rbind(MeanAPE = colMeans(APE, na.rm = TRUE),
                               MeanSPE = colMeans(SPE, na.rm = TRUE),
                               Coverage = colMeans(P.COVER, na.rm = TRUE),
                               AvgLength = colMeans(P.LENGTH, na.rm = TRUE))
       colnames(TAB3) <- PRED.NAMES
  
#####################################################################
##### TABLE 4 : MCMC DIAGNOSTICS ####################################
#####################################################################
  
##### SUMMARY FUNCTION ##############################################

               SUMFUN <- function(x)
                      {
                    x <- x[is.finite(x)]
                    if(length(x) == 0)
                      {
                return(c(Min = NA,
                         Q1 = NA,
                         Median = NA,
                         Mean   = NA,
                         Q3     = NA, 
                         Max    = NA))
                      }
                     c(Min    = min(x),
                       Q1 = quantile(x, 0.25),
                       Median = median(x),
                       Mean = mean(x),
                       Q3 = quantile(x, 0.75),
                       Max = max(x))
                      }
  
             ACC.IMH  <- SUMFUN(MCMC.SUM[, "IMH.ACC"])
             ACC.RWM  <- SUMFUN(MCMC.SUM[, "RWM.ACC"])
  
             ESS.IMH  <- SUMFUN(MCMC.SUM[, "IMH.ESS"])
             ESS.RWM  <- SUMFUN(MCMC.SUM[, "RWM.ESS"])
  
             RESS.IMH <- SUMFUN(MCMC.SUM[, "IMH.RESS"])
             RESS.RWM <- SUMFUN(MCMC.SUM[, "RWM.RESS"])
  
             GR.IMH   <- SUMFUN(MCMC.SUM[, "IMH.GR"])
             GR.RWM   <- SUMFUN(MCMC.SUM[, "RWM.GR"])
  
##### TABLE 4
             
                  TAB4 <- rbind(c(ACC.IMH,  ACC.RWM),
                                c(ESS.IMH,  ESS.RWM),
                                c(RESS.IMH, RESS.RWM),
                                c(GR.IMH,   GR.RWM))
  
        rownames(TAB4) <- c("Acceptance Rate", "ESS", "Relative ESS", "Rhat")
  
        colnames(TAB4) <- c("IMH.Min",
                            "IMH.Q1",
                            "IMH.Median",
                            "IMH.Mean",
                            "IMH.Q3",
                            "IMH.Max",
                            "RWM.Min",
                            "RWM.Q1",
                            "RWM.Median",
                            "RWM.Mean",
                            "RWM.Q3",
                            "RWM.Max")
  
                  TAB4 <- as.data.frame(round(TAB4, 4))
  
#####################################################################
##### OUTPUT ########################################################
#####################################################################
  
                 return(
                   list(Table.Estimates = round(TAB1,4),
                        Table.Intervals = round(TAB2,4),
                        Table.Prediction = round(TAB3,4),
                        Table.MCMC = TAB4,
                        Estimates = EST,
                        AbsoluteError = AE,
                        RMSE.KM = RMSE.KM,
                        Coverage = COVER,
                        Length = LENGTH,
                        APE = APE,
                        SPE = SPE,
                        PredictionCoverage = P.COVER,
                        PredictionLength = P.LENGTH,
                        SuccessfulReplications = SUCCESS,
                        MCMC.Diagnostics = MCMC.SUM,
                        Censored = NCENS))
                      }

