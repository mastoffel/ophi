# non-parametric bootstrap to assess stability of final version

# permutation test for ophi_code data
library(MASS)

ophi_code <- read.csv("data/ophi_code.csv", colClasses = c(rep("numeric",2),rep("factor",98)))
# mod_code <- glm(day_falling ~., data = ophi_code)
# code_aic <- stepAIC(mod_code, direction = "backward")
# df_aic_code <- model.frame(code_aic)
# based on permutation ---------------------------------------------------------
bootstrap_ophi <- function(...) {
        # randomization of allels according to probabilities
        ophi_code_temp <- (ophi_code[sample(nrow(ophi_code), size = round(nrow(ophi_code) * 0.8)), ])
        
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
        
        final_names <- names(df)
}

# parallel
library(parallel)

cl <- makeCluster(detectCores()-1)
#get library support needed to run the code
clusterEvalQ(cl,library(MASS))

#put objects in place that might be needed for the code
clusterExport(cl, c("bootstrap_ophi", "ophi_code"))

all <- parLapply(cl, 1:10000, bootstrap_ophi)
write(unlist(all), "all_names.txt")

stopCluster(cl)
