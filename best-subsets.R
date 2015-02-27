# best subsets selection

# ophi as numeric
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

# find best subsets for each size
library(leaps)
regfit.full <- regsubsets(day_falling ~., ophi_het, method = "backward", nvmax = 99)
reg.summary <- summary(regfit.full)

# plotting 
par(mfrow = c(2,2))
plot(reg.summary$rss, xlab = "Number of Variables", ylab = "RSS", type = "l")
plot(reg.summary$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq", type = "l")
which.max(reg.summary$adjr2)
points(39, reg.summary$adjr2[39], col = "red", cex = 2, pch = 20)
plot(reg.summary$cp, xlab = "Number of Variables", ylab = "Mallows Cp", type = "l")
which.min(reg.summary$cp)
points(16, reg.summary$cp[16], col = "red", cex=2, pch = 20)
which.min(reg.summary$bic)
plot(reg.summary$bic, xlab = "Number of Variables", ylab = "BIC", type = "l")
points(1, reg.summary$bic[1], col = "red", cex = 2, pch = 20)
# title("Best subsets for the additive model", outer = TRUE)

# cross validation

set.seed(1)
train <-  sample(c(TRUE,FALSE), nrow(ophi_het), rep = TRUE, prob = c(0.5,0.5))
test <- (!train)

# apply regsubsets to the training set

regfit.best <- regsubsets(day_falling ~. ,data = ophi_het[train, ], nvmax = 99,
                          method = "backward", really.big = TRUE)

test.mat <- model.matrix(day_falling ~. , data = ophi_het[test, ])

# loop over each best subset

val.errors <- rep(NA, 99)
for ( i in 1:99) {
        coefi <- coef(regfit.best, id = i)
        pred <- test.mat[, names(coefi)] %*% coefi
        val.errors[i] <- mean((ophi_het$day_falling[test]-pred)^2)
}
plot(val.errors)

# writing own predict method

predict.regsubsets <- function(object, newdata, id, ...) {
        form <- as.formula(object$call[[2]])
        mat <- model.matrix(form, newdata)
        coefi <- coef(object, id = id) 
        xvars <- names(coefi)
        mat[, xvars] %*% coefi
}


# regsubsets with full model
regfit.best <- regsubsets(day_falling ~. ,data = ophi_het, nvmax = 98, method = "backward")
coef(regfit.best, 4) # exactly the crawley model

k <- 10
set.seed(1)
folds <- sample(1:k, nrow(ophi_het), replace = TRUE)
cv.errors <- matrix(NA, k, 98, dimnames=list(NULL, paste(1:98)))

# cross validation
for(j in 1:k) {
        best.fit <- regsubsets(day_falling ~. , data = ophi_het[folds!=j, ], nvmax = 98, method = "backward")
        
        for (i in 1:98) {
                pred = predict.regsubsets(best.fit, ophi_het[folds == j, ], id = i)
                cv.errors[j, i] <- mean((ophi_het$day_falling[folds==j] - pred)^2)
        }
}

mean.cv.errors <- apply(cv.errors, 2, mean)
plot(mean.cv.errors, xlab="best subset size")
