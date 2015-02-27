# analysis ophi data


# code model, genotype encoded as categorical variable -------------------------

# stepAIC 
library(MASS)
# load data
ophi_code <- read.csv("data/ophi_code.csv", colClasses = c(rep("numeric",2),rep("factor",98)))
ophi_add <- read.csv("data/ophi_add.csv")
ophi_het <- read.csv("data/ophi_het.csv")

# create empty complete table
complete_table <- data.frame(matrix(NA, nrow=100, ncol=5))
names(complete_table) <- c("Estimate", "Deviance Explained", "F", "P (F-test)", "P (Chisquared-test)")

# model selection code----------------------------------------------------------
# modelfit code
mod_code <- glm(day_falling ~., data = ophi_code)

## AIC backward stepwise ##
mod_aic_code <- stepAIC(mod_code, direction = "backward")
# final aic model variables data frame
df_aic_code <- model.frame(mod_aic_code)

# get deletion test table
library(minmodelr)
code_aic_table <- DelTestVar(df_aic_code)
# fill in complete SNP table
fulltab_code_aic <- complete_table
row.names(fulltab_code_aic) <- names(ophi_code)
fulltab_code_aic[row.names(code_aic_table), ] <- code_aic_table


## stepwise regression (crawley) ##
mod_crawley_code <- MinMod(ophi_code)
df_crawley_code <- mod_crawley_code[[1]]
code_crawley_table <- DelTestVar(df_crawley_code)

# fill in complete SNP table
fulltab_code_crawley <- complete_table
row.names(fulltab_code_crawley) <- names(ophi_code)
fulltab_code_crawley[row.names(code_crawley_table), ] <- code_crawley_table


# model selection additive------------------------------------------------------
# modelfit add
mod_add <- glm(day_falling ~., data = ophi_add)

## AIC backward stepwise ##
mod_aic_add <- stepAIC(mod_add, direction = "backward")
# final aic model variables data frame
df_aic_add <- model.frame(mod_aic_add)

# get deletion test table
library(minmodelr)
add_aic_table <- DelTestVar(df_aic_add)
# fill in complete SNP table
fulltab_add_aic <- data.frame(matrix(NA, nrow=100, ncol=5),
                               row.names = names(ophi_add))
names(fulltab_add_aic) <- names(add_aic_table)
fulltab_add_aic[row.names(add_aic_table), ] <- add_aic_table


## stepwise regression (crawley) ##
mod_crawley_add <- MinMod(ophi_add)
df_crawley_add <- mod_crawley_add[[1]]
add_crawley_table <- DelTestVar(df_crawley_add)

# fill in complete SNP table
fulltab_add_crawley <- complete_table
row.names(fulltab_add_crawley) <- names(ophi_add)
fulltab_add_crawley[row.names(add_crawley_table), ] <- add_crawley_table

# model selection het------------------------------------------------------
# modelfit het
mod_het <- glm(day_falling ~., data = ophi_het)

## AIC backward stepwise ##
mod_aic_het <- stepAIC(mod_het, direction = "backward")
# final aic model variables data frame
df_aic_het <- model.frame(mod_aic_het)

# get deletion test table
library(minmodelr)
het_aic_table <- DelTestVar(df_aic_het)
# fill in complete SNP table
fulltab_het_aic <- data.frame(matrix(NA, nrow=100, ncol=5),
                              row.names = names(ophi_het))
names(fulltab_het_aic) <- names(het_aic_table)
fulltab_het_aic[row.names(het_aic_table), ] <- het_aic_table


## stepwise regression (crawley) ##
mod_crawley_het <- MinMod(ophi_het)
df_crawley_het <- mod_crawley_het[[1]]
het_crawley_table <- DelTestVar(df_crawley_het)

# fill in complete SNP table
fulltab_het_crawley <- complete_table
row.names(fulltab_het_crawley) <- names(ophi_het)
fulltab_het_crawley[row.names(het_crawley_table), ] <- het_crawley_table




library(broom)
test <- tidy(mod_aic_code)


