# dimension reduction approach


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

# MCA multiple correspondance analysis on categorical data----------------------

library(FactoMineR)

new_ophi <- ophi_code[3:100]
mca1 <- MCA(new_ophi)
plot.MCA(mca1, invisible=c("ind"))
scores <- as.data.frame(mca1$ind$coord)

scores$day_falling <- ophi_code$day_falling

mod <- glm(day_falling ~. , data =scores)
summary(mod)



ophi_het <- read.csv("./data/joedata.csv",
                     colClasses = c(rep("numeric",2),rep("factor",98)), 
                     header = TRUE)

ophi_het[, 2] <- ophi_code$disc_size
names(ophi_het)[2] <- "disc_size"

library(FactoMineR)

new_ophi <- ophi_het[3:100]
mca1 <- MCA(new_ophi)
plot.MCA(mca1, invisible=c("ind"))
scores <- as.data.frame(mca1$ind$coord)

scores$day_falling <- ophi_code$day_falling

mod <- glm(day_falling ~. , data =scores)
summary(mod)


# explore dimensions

res = dimdesc(mca1, axes=2, proba=0.001)
res


# PCA on het data---------------------------------------------------------------
new_ophi <- ophi_add[3:100]
fit <- princomp(new_ophi, cor = TRUE)
summary(fit)

scores <- cbind(ophi_add["day_falling"], fit$scores[,1:10])

mod <- glm(day_falling ~. , data =scores)
summary(mod)

# PCA with rotation-------------------------------------------------------------
library(psych)

new_ophi <- ophi_het[3:100]

fit <- principal(new_ophi, nfactors=10, rotate="varimax")
scores <- cbind(ophi_add["day_falling"], fit$scores)

mod <- glm(day_falling ~. , data =scores)
summary(mod)

# FA ---------------------------------------------------------------------------
new_ophi <- ophi_add[3:100]

fit <- factanal(new_ophi, 10, rotation = "varimax")
# print(fit, digits=2, cutoff=.1, sort=TRUE)

scores <- matrix(rep(NA, 288*10), ncol = 10)
for (i in 1:10) {
scores[, i] <- as.vector(fit$loadings[,i]) %*% t(new_ophi)
}

scores <- cbind(ophi_add["day_falling"], scores)
mod <- glm(day_falling ~. , data =scores)
summary(mod)



cats = apply(new_ophi, 2, function(x) nlevels(as.factor(x)))

mca1_vars_df = data.frame(mca1$var$coord, Variable = rep(names(cats), cats))

# data frame with observation coordinates
mca1_obs_df = data.frame(mca1$ind$coord)


ggplot(data = mca1_obs_df, aes(x = Dim.1, y = Dim.2)) +
        geom_hline(yintercept = 0, colour = "gray70") +
        geom_vline(xintercept = 0, colour = "gray70") +
        geom_point(colour = "gray50", alpha = 0.7) +
        geom_density2d(colour = "gray80") +
        geom_text(data = mca1_vars_df, 
                  aes(x = Dim.1, y = Dim.2, 
                      label = rownames(mca1_vars_df), colour = Variable)) +
        ggtitle("MCA plot of variables using R package FactoMineR") +
        scale_colour_discrete(name = "Variable")