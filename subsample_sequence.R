# mean model fits for an increasing number of subsampling

# non-parametric bootstrap to assess stability of final version

# permutation test for ophi_code data
library(MASS)
library(plyr)
ophi_code <- read.csv("data/ophi_code.csv", colClasses = c(rep("numeric",2),rep("factor",98)))
# mod_code <- glm(day_falling ~., data = ophi_code)
# code_aic <- stepAIC(mod_code, direction = "backward")
# df_aic_code <- model.frame(code_aic)
# based on permutation ---------------------------------------------------------

subsamp_seq <- function(...) {
        
        bootstrap_ophi <- function(nsub, ...) {
                # randomization of allels according to probabilities
                ophi_code_temp <- (ophi_code[sample(nrow(ophi_code), size = round(nrow(ophi_code) * nsub)), ])
                
                check_contrasts <- function(snp) {
                        if (length(unique(snp)) == 1) {
                                snp <- 0
                                return(snp)
                        }
                        snp <- 1
                }
                
                f_levels <- unlist(lapply(ophi_code_temp, check_contrasts))
                
                ophi <- ophi_code_temp[f_levels != 0]
                
                # modelfit code
                mod_code <- glm(day_falling ~., data = ophi)
                
                # AIC backward stepwise 
                mod_aic_code <- stepAIC(mod_code, direction = "backward")
                
                # model frame
                df <- model.frame(mod_aic_code)
                
                # get final aic
                final_aic <- summary(mod_aic_code)$aic
                # get deviance
                null_dev <- summary(mod_aic_code)$null.deviance
                res_dev <- summary(mod_aic_code)$deviance
                # get final r2
                r2 <- ((summary(mod_aic_code)$null.deviance) - (summary(mod_aic_code)$deviance)) / (summary(mod_aic_code)$null.deviance)
                
                # get final adj r2
                npred <- length(summary(mod_aic_code)$contrasts)
                nobs <- nrow(ophi)
                adjr2 <- 1 - ((1 - r2^2) * (nobs - 1) / (nobs - npred - 1))
                
                # null
                mod_out <- list(nsub = nsub, aic = final_aic, null_dev = null_dev, res_dev = res_dev,
                                r2 = r2, adjr2 = adjr2, npred = npred, nobs = nobs)
                out <- data.frame(mod_out)
        }
        subsampled <- lapply(c(0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,0.97), bootstrap_ophi)
}

# parallel
library(parallel)

cl <- makeCluster(detectCores())
#get library support needed to run the code
clusterEvalQ(cl,library(MASS))

#put objects in place that might be needed for the code
clusterExport(cl, c("subsamp_seq", "ophi_code"))

all <- parLapply(cl, 1:2, subsamp_seq)
# combining into data.frame
results <- ldply(unlist(all, recursive = FALSE))

write.csv(results, "results.csv")

stopCluster(cl)
