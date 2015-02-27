# permutation test for ophi_code data
library(MASS)

ophi_code <- read.csv("data/ophi_code.csv", colClasses = c(rep("numeric",2),rep("factor",98)))
# mod_code <- glm(day_falling ~., data = ophi_code)
# code_aic <- stepAIC(mod_code, direction = "backward")
# df_aic_code <- model.frame(code_aic)
# based on permutation ---------------------------------------------------------
rand_aic <- function(...) {

        # randomization of allels according to probabilities
        rand_allelfreq <- function(x) {
                # x is a snp vector
                y <- sample(x)  
        }
        # create new df
        ophi_code_temp <- ophi_code
        ophi_code_temp[3:ncol(ophi_code_temp)] <- as.data.frame(apply(ophi_code[3:ncol(ophi_code)],2, rand_allelfreq))
        
        # modelfit code
        mod_code <- glm(day_falling ~., data = ophi_code_temp)
        
        # AIC backward stepwise 
        mod_aic_code <- stepAIC(mod_code, direction = "backward")
        
        # get final aic
        final_aic <- summary(mod_aic_code)$aic
        # get final r2
        r2 <- ((summary(mod_aic_code)$null.deviance) - (summary(mod_aic_code)$deviance)) / (summary(mod_aic_code)$null.deviance)
        # get final adj r2
        npred <- length(summary(mod_aic_code)$contrasts)
        nobs <- nrow(ophi_code)
        adjr2 <- 1 - ((1 - r2^2) * (nobs - 1) / (nobs - npred - 1))
        
        out <- list("aic" = final_aic, "r2" = r2, "adjr2" = adjr2)
        
}

# parallel
library(parallel)

cl <- makeCluster(detectCores()-1)
#get library support needed to run the code
clusterEvalQ(cl,library(MASS))

#put objects in place that might be needed for the code
clusterExport(cl, c("rand_aic", "ophi_code"))

all <- parLapply(cl, 1:5000, rand_aic)

all_aic <- sapply(all, function(x) x$aic)
all_r2 <- sapply(all, function(x) x$r2)
all_adjr2 <- sapply(all, function(x) x$adjr2)

write(all_aic, "aics.txt")
write(all_r2, "r2.txt")
write(all_adjr2, "adjr2.txt")

stopCluster(cl)
