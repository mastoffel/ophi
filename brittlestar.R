## checking MinMod with joes data
# extract coefficients from final models --> fitted values
# fit every single snp in a model
# + 


# load data---------------------------------------------------------------------
ophi_code <- read.table("./data/Ophionotus_GLM_4.txt", 
                      colClasses = c(rep("numeric",2),rep("factor",98)), 
                      header = TRUE, row.names = NULL)

# ophi as het, numeric
ophi_het <- read.csv("./data/joedata.csv",
                     colClasses = c(rep("numeric",2),rep("numeric",98)), 
                     header = TRUE)

ophi_het[, 2] <- ophi_code$disc_size
names(ophi_het)[2] <- "disc_size"

# ophi as numeric
ophi_add <-  read.table("./data/Ophionotus_GLM_4.txt", 
                               colClasses = c(rep("numeric",100)), 
                               header = TRUE, row.names = NULL)
ophi_add[, 3:ncol(ophi_add)] <- ophi_add[, 3:ncol(ophi_add)] - 2


# single variable deletion testing----------------------------------------------
library(minmodelr)
single_add <- DelTestVar(ophi_add)
single_het <- DelTestVar(ophi_het)


# methods for best model selection----------------------------------------------

library(leaps)
day_falling <- ophi_add$day_falling
inp_mat_add <- ophi_add[, -1]
inp_mat_het <- ophi_het[, -1]
inp_mat_code <- ophi_code[, -1]



# doesnt work for het because of factors, turn to numeric
inp_mat_het[] <- lapply(inp_mat_het, as.numeric)

regsub <- regsubsets(x = inp_mat_het, y = day_falling, 
                     nbest = 1, nvmax = 30, method = "backward",
                     really.big = T)


plot(regsub, scale = "adjr2", main = "Adjusted R^2 - best models for 1-30 predictors")
plot(regsub, scale = "r2", main = "R^2 - best models for 1-30 predictors")
plot(regsub, scale = "bic", main = "BIC - best models for 1-30 predictors")
plot(regsub, scale = "Cp", main = "Mallow큦 Cp - best models for 1-30 predictors")

# using bestglm ----------------------------------------------------------------

library(bestglm)

# additive model with AIC-------------------------------------------------------

Xy <- cbind(inp_mat_add, day_falling)
# AIC backwards algorithm
best_add <- bestglm(Xy, IC = "AIC", method = "backward")
# glm object containing final best model
bestmod_add <- best_add$BestModel
# preparing for full_table function with uses DelTestVar for
# single deletion testing and puts the final SNP큦 into the full SNP table (98)
varnames <- names(bestmod_add$coef)[-1]
bestdf_add <- cbind(day_falling, inp_mat_add[varnames])

source("full_table.R")
# reference table should just contain snp row names
table_add_AIC <- full_table(bestdf_add, inp_mat_add[-1])


# additive model with BIC (all vars deleted!!) ---------------------------------
Xy <- cbind(inp_mat_add, day_falling)
# AIC backwards algorithm
best_add <- bestglm(Xy, IC = "BIC", method = "backward")

# additive model with CV -------------------------------------------------------
Xy <- cbind(inp_mat_add, day_falling)
# AIC backwards algorithm
best_add <- bestglm(Xy, IC = "CV")




library(glmulti)

# LASSO -----------------------------------------------------------------
# paper "LASSO model selection with post-processing for a genome-wide association study data set"
library(glmnet)


mod <- glmnet(x = as.matrix(inp_mat), y = day_falling, family = "gaussian")

# 1) Holgers method: testing every single SNP in a model and then
# do a false discovery rate correction
# qvalue or benjamini and hochberg

ophi <- ophi_add

results <- data.frame(row.names = names(ophi[3:ncol(ophi)]))

for (i in 3:ncol(ophi)) {
        # fit model with disc size and one SNP as indepedent vars
        fit2 <- glm(ophi[, 1] ~ ophi[, 2] + ophi[, i])
        fit <- glm(ophi[, 1] ~ ophi[, 2])
        model_comp <- anova(fit, fit2, test = "F")
        
        # get p values in vector
        pval <- model_comp[2, 6]
        fval <- model_comp[2, 5]
        # extract pvalue from SNP
        # pval <- summary(fit)$coef[, "Pr(>|t|)"][3]
        
        results[i-2, "F"] <- fval
        results[i-2, "pval"] <- pval
}

# false discovery rate with qval
library("qvalue")

# lambda = 0 is Benjamini&Hochberg 1995; 
qval <- qvalue(results$pval,  pi0.method="smoother")

results$qval <- qval$qvalue
write.csv(results, "single_models_het.csv")

# p-value cut off for a given false discovery rate level
max(qval$pvalues[qval$qvalues <= 0.5])


# Joe: First try 98 models with disc size + each snp

fit <- glm(ophi[, 1])
fit2 <- glm(ophi[, 1] ~ ophi[, 2] + ophi[, i])



# stepAIC backward-------------------------------------------------------------
ophi <- ophi_code
ophi <- ophi_het
ophi <- ophi_add

# databin <- read.csv("C:/Users/Martin/Studium/Ophionotus/files/joedata.csv", header = TRUE)
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

# deletion test every variable
vals <- DelTestVar(ophired)

# reduce data frame 
val_red <- vals[2:(length(modnames)+1), ]
val_red2 <- (val_red[, c("F","P (F-test)")])

# give the right names (not factor levels)
row.names(val_red2) <- modnames

# sort within all SNP큦
final_snp_df <- data.frame("F" = rep(NA, 98), "pval" = rep(NA, 98), row.names = names(ophi)[3:ncol(ophi)])

# fill in the final AIC model
final_snp_df[row.names(val_red2), ] <- val_red2
write.csv(final_snp_df, "stepAICbackHetFac.csv")


# minmodelr------------------------------------------------------------------
library(minmodelr)
ophi_list <- list(ophi_add, ophi_het)

all_models <- lapply(ophi_list, MinMod)



# test crawley models again full models

# dependent var
day_falling <- ophi_add[, 1]

# additive
fullmod_add <- glm(day_falling ~., data=ophi_add[, 2:ncol(ophi_add)])
crawley_add <- all_models[[1]][[1]]
add_mod <- glm(day_falling ~., data = crawley_add[, 2:ncol(crawley_add)])
anova_add <- anova(fullmod_add,add_mod, test = "F")

# het
fullmod_het <- glm(day_falling ~., data = ophi_het[, 2:ncol(ophi_het)])
crawley_het <- all_models[[3]][[1]][-c(3:11)]
het_mod <- glm(day_falling ~., data = crawley_het[, 2:ncol(crawley_het)])
anova_het <- anova(fullmod_het, het_mod, test = "F")


# plot anova for every variable ------------------------------------------------

give.n <- function(x){
        return(c(y = mean(x), label = length(x)))
}

code_df <- all_models[[2]][[1]]
het_df <- all_models[[3]][[1]]
add_df <- all_models[[1]][[1]]

# change to factor for proper boxplots
add_df[, 2:ncol(add_df)] <- lapply((add_df[, 2:ncol(add_df)]), as.factor)
all_plots <- list()
for (i in names(add_df)[2:ncol(add_df)]) {

        # plot boxplots for 
        p <- ggplot(add_df, aes_string(y = "depVar", x = i)) + 
                geom_boxplot() + 
                stat_summary(fun.data = give.n, geom = "text", size = 10) +
                ylab("day falling") +
                theme_classic(base_size = 20)
        (assign(i, p))
}

for (i in names(add_df)[2:ncol(add_df)]) {
        
        # plot boxplots for 
        p <- ggplot(add_df, aes_string(y = "depVar", x = i)) + 
                geom_point() +
                stat_smooth(method = lm) + 
               stat_summary(fun.data = give.n, geom = "text", size = 10) +
                ylab("day falling") +
                theme_classic(base_size = 20)
        (assign(i, p))
}

source("multiplot.R")
multiplot(het_imp_14, het_imp_3, het_imp_33, het_imp_47, het_imp_53,
          het_imp_58, het_imp_66,het_imp_75,het_imp_79,het_imp_95, cols = 2)

multiplot(code_imp_1, code_imp_10, code_imp_59, code_imp_61, code_imp_88, code_imp_97, cols = 2)


# data frame 
add_df <- all_models[[1]][[1]]
add_mod <- all_models[[1]][[2]]
add_df2 <- data.frame(code_imp_1 = seq(min(het_df$code_imp_1), max(het_df$code_imp_1), 0.02),
                      code_imp_10 = mean(het_df$code_imp_10),
                      code_imp_59 = mean(het_df$code_imp_59),
                      code_imp_61 = mean(het_df$code_imp_61),
                      code_imp_88 = mean(het_df$code_imp_88),
                      code_imp_97 = mean(het_df$code_imp_97))

add_df2$code_imp_1_pred <- predict(add_mod, add_df2, type = "response")

                                   
bestmod <- MinMod(ophi) 
modnames <- attr(terms(bestmod[[2]]), "term.labels")


# get fitted values with predict
df <- bestmod[[1]]
df2 <- data.frame(het_imp_3 = seq(min(df$het_imp_3), max(df$het_imp_3, 0.02), het_imp_14 = mean(het_imp_14)))



# construct data frame with F and P values 
# of each variable in the final Model-----------------
ophired <- cbind(ophi[1], ophi[, modnames])

# deletion test every variable
vals <- DelTestVar(ophired)

# reduce data frame 
val_red <- vals[2:(length(modnames)+1), ]
val_red2 <- (val_red[, c("F","P (F-test)")])

# give the right names (not factor levels)
row.names(val_red2) <- modnames

# sort within all SNP큦
final_snp_df <- data.frame("F" = rep(NA, 98), "pval" = rep(NA, 98), row.names = names(ophi)[3:ncol(ophi)])

# fill in the final AIC model
final_snp_df[row.names(val_red2), ] <- val_red2

# het data
mod <- bestmod[[2]]
data <- bestmod[[1]]




# using granova package to visualize
library(granovaGG)
p <- granovagg.1w(data=data[, 1], 
        group = ((data[, 7])))
print(p)
write.csv(final_snp_df, "crawleyhet.csv")

bestmod2 <- minmodelr(ophi)

library(leaps)
leaps <- regsubsets(ophi[, 1] ~. , data = ophi[, 3:numVar], nbest = 2)

## num vars
numVar <- ncol(datafac)
datafac[, 3:numVar] <- lapply(data[, 3:numVar], factor)
model <- MinMod(datafac)
testmodel <- glm(datafac[, 1] ~. , data = datafac[, 2:numVar])
models <- glmulti(day_falling ~., data = subset(datafac, select = 1:20), method = "g", level=1)

## other best subset approaches
# glmulti, bestglm, leaps, step, stepAIC      (all based on AIC!)
library(glmulti)
mod <- glm(day_falling ~., data = ophi_add)
results <- glmulti(mod, method = "l", level = 1)

# lasso
library(glmnet)
inp_add <- as.matrix(ophi_het[2:ncol(ophi_het)])
y <- ophi_add[, 1]

grid <- 10^seq(10, -2, length = 100)

set.seed(1)
train <- sample(1:nrow(ophi_add), nrow(ophi_add)/2)
test <- (-train)
y.test <- y[test]

mod <- glmnet(x = inp_add[train, ], y = y[train], alpha = 1, lambda = grid)
plot(mod)

set.seed(1)
cv.out <- cv.glmnet(inp_add[train, ], y[train], alpha = 1)
plot(cv.out)
bestlam <- cv.out$lambda.min
lasso.pred <- predict(mod, s = bestlam, newx = inp_add[test, ])
mean((lasso.pred-y.test)^2)

# mean with just the intercept
mean((mean(y[train]) - y.test)^2)

# variable selection by lasso
out <- glmnet(inp_add, y, alpha=1, lambda = grid)
lasso.coef <- predict(out, type="coefficients", s= bestlam)
lasso.coef
