# AIC backwards selection approach

# stepAIC backward-------------------------------------------------------------
ophi_code <- read.table("./data/Ophionotus_GLM_4.txt", 
                        colClasses = c(rep("numeric",2),rep("factor",98)), 
                        header = TRUE, row.names = NULL)
ophi_add <-  read.table("./data/Ophionotus_GLM_4.txt", 
                        colClasses = c(rep("numeric",100)), 
                        header = TRUE, row.names = NULL)
ophi_add[, 3:ncol(ophi_add)] <- ophi_add[, 3:ncol(ophi_add)] - 2
ophi_het <- read.csv("./data/joedata.csv",
                     colClasses = c(rep("numeric",2),rep("numeric",98)), 
                     header = TRUE)
ophi_het[, 2] <- ophi_code$disc_size
names(ophi_het)[2] <- "disc_size"

# databin <- read.csv("C:/Users/Martin/Studium/Ophionotus/files/joedata.csv", header = TRUE)

ophi <- ophi_het

numVar <- ncol(ophi)

# full model fit
fit <- glm(ophi[, 1] ~. , data = ophi[, 2:numVar])

# stepwise model selection AIC
library(MASS)
step <- stepAIC(fit, direction = "backward")

# extract final model
formula(step)

# extract variable names
modnames <- attr(terms(step), "term.labels")

# construct data frame with F and P values of each variable in the final Model
ophired <- cbind(ophi[1], ophi[, modnames])

# glmulti
library(glmulti)

mod <- glmulti(ophired[, 1] ~., data = ophired, method = "d", level = 1)








# deletion test every variable
vals <- DelTestVar(ophired)

# reduce data frame 
val_red <- vals[2:(length(modnames)+1), ]
val_red2 <- (val_red[, c("F","P (F-test)")])

# give the right names (not factor levels)
row.names(val_red2) <- modnames

# sort within all SNP´s
final_snp_df <- data.frame("F" = rep(NA, 98), "pval" = rep(NA, 98), row.names = names(ophi)[3:ncol(ophi)])

# fill in the final AIC model
final_snp_df[row.names(val_red2), ] <- val_red2
write.csv(final_snp_df, "stepAICbackHetFac.csv")



# pca
library(homals)

# random forests
library(randomForest)
library(AUCRF)
library(VSURF)
mod <- AUCRF(ophi[, 1] ~. , data = ophi[, 2:numVar])

mod <- VSURF(formula = ophi[, 1] ~. , data = ophi[, 2:numVar])
