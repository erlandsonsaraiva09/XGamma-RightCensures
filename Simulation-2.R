###############################################################################
########## SIMULATION STUDY 2 DATA ############################################
###############################################################################

# Run Fit.functionsXGamma.R
# Run Monte-Carlo_functions.R

set.seed(19)


library(parallel)
library(coda)
library(survival)
library(progress)

NCORES <- 14

MC.total <- 1000

MC.split <- rep(floor(MC.total/NCORES), NCORES)

MC.split[NCORES] <- MC.total - sum(MC.split[-NCORES])

MC.split

sum(MC.split)

system.time(
           {
       OUT <- mclapply(1:NCORES,
              function(i){MonteCarlo.XGamma(theta = 0.6,
                                             n = 25,
                                             p.cens = 0.50,
                                             MC = MC.split[i],
                                             K.pred = 30)
                          },
   mc.cores = NCORES)})

##### MONTE CARLO

Merge.MonteCarlo <- function(OUT)
                 {
  
###############################################################
##### CENSORING ###############################################
###############################################################
  
           NCENS <- unlist(lapply(OUT, function(x) x$Censored))
  
###############################################################
##### ESTIMATES ###############################################
###############################################################
  
             EST <- do.call(rbind, lapply(OUT, function(x) x$Estimates))
              AE <- do.call(rbind,lapply(OUT, function(x) x$AbsoluteError))
  
         RMSE.KM <- do.call(rbind, lapply(OUT, function(x) x$RMSE.KM))
  
            TAB1 <- rbind(MeanEstimate = colMeans(EST, na.rm = TRUE),
                          MAE = colMeans(AE, na.rm = TRUE),
                          MeanRMSE.KM = colMeans(RMSE.KM, na.rm = TRUE))
  
###############################################################
##### INTERVALS ###############################################
###############################################################
  
           COVER <- do.call(rbind,lapply(OUT, function(x) x$Coverage))
  
          LENGTH <- do.call(rbind, lapply(OUT, function(x) x$Length))
  
            TAB2 <- rbind(Coverage = colMeans(COVER, na.rm = TRUE),
                          AvgLength = colMeans(LENGTH, na.rm = TRUE))
  
###############################################################
##### PREDICTION ##############################################
###############################################################
  
             APE <- do.call(rbind, lapply(OUT, function(x) x$APE))
  
             SPE <- do.call(rbind, lapply(OUT, function(x) x$SPE))
  
         P.COVER <- do.call(rbind, lapply(OUT, function(x) x$PredictionCoverage))
  
        P.LENGTH <- do.call(rbind, lapply(OUT, function(x) x$PredictionLength))
  
            TAB3 <- rbind(MeanAPE = colMeans(APE, na.rm = TRUE),
                          MeanSPE = colMeans(SPE, na.rm = TRUE),
                          Coverage = colMeans(P.COVER, na.rm = TRUE),
                          AvgLength = colMeans(P.LENGTH, na.rm = TRUE))
  
###############################################################
##### MCMC ####################################################
###############################################################
  
            MCMC <- do.call(rbind, lapply(OUT, function(x) x$MCMC.Diagnostics))
  
            TAB4 <- as.data.frame(round(apply(MCMC, 2, mean, na.rm = TRUE), 4))
  
###############################################################
##### OUTPUT ##################################################
###############################################################
  
           return(
             list(Table.Estimates = round(TAB1,4),
                  Table.Intervals = round(TAB2,4),
                  Table.Prediction = round(TAB3,4),
                  Table.MCMC = TAB4,
                  Mean.Censored = mean(NCENS, na.rm = TRUE),
                  Censored = NCENS))}

### RESULTS

              MC <- Merge.MonteCarlo(OUT)

             mean(MC$Censored)
             
          MC.EST <- MC$Table.Estimates
          MC.EST

    MC.INTERVALS <- MC$Table.Intervals
    MC.INTERVALS
    
              MT <- as.matrix(MC$Table.Prediction)
              MT <- MT[2:4,]
          MT[1,] <- round(sqrt(MT[1,]), 4)
              MT




