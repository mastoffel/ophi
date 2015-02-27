# permutation test for ophi_code data
library(MASS)

ophi_code <- read.csv("data/ophi_code.csv", colClasses = c(rep("numeric",2),rep("factor",98)))
mod_code <- glm(day_falling ~., data = ophi_code)
code_aic <- stepAIC(mod_code, direction = "backward")
df_aic_code <- model.frame(code_aic)

# based on permutation ---------------------------------------------------------

perm_test <- function(df_aci_code) {
        
# construct snp data frame
snp_df <- df_aic_code[2:ncol(df_aic_code)]
# scramble genotypes
snp_df_temp <- as.data.frame(apply(snp_df,2, sample))
# bind day_falling to the df
df_aic_code_temp <- cbind(df_aic_code[1], snp_df_temp)
# fit model
mod_code_temp <- glm(day_falling ~., data = df_aic_code_temp)
# aic
aic <- summary(mod_code_temp)$aic
# get r2
r2 <- ((summary(mod_code_temp)$null.deviance) - (summary(mod_code_temp)$deviance)) / (summary(mod_code_temp)$null.deviance)

out <- list("aic" = aic, "r2" = r2)

}

all <- replicate(10000, perm_test(df_aic_code))


all_aic <- unlist(all[1, ])
all_r2 <- unlist(all[2, ])

# histograms
hist(all_aic, xlim = range(2400, 2650), main = "10000 permutations of 49 SNP genotypes",
     xlab = "AIC")
abline(v = 2447)

hist(all_r2, xlim = range(0, 0.6), main = "10000 permutations of 49 SNP genotypes", 
     xlab = "r2")
abline(v = 0.59)
