# joes table values


# ophi as numeric
ophi_add <-  read.table("./data/Ophionotus_GLM_4.txt", 
                        colClasses = c(rep("numeric",100)), 
                        header = TRUE, row.names = NULL)
ophi_add[, 3:ncol(ophi_add)] <- ophi_add[, 3:ncol(ophi_add)] - 2

# ophi as het, numeric
ophi_het <- read.csv("./data/joedata.csv",
                     colClasses = c(rep("numeric",2),rep("numeric",98)), 
                     header = TRUE)

ophi_het[, 2] <- ophi_add$disc_size
names(ophi_het)[2] <- "disc_size"


library(minmodelr)
ophi_list <- list(ophi_add, ophi_het)
all_models <- lapply(ophi_list, MinMod)

# extract df´s
ophi_crawley_add <- all_models[[1]][[1]]
ophi_crawley_het <- all_models[[2]][[1]]
names(ophi_crawley_het) <- names(all_models[[2]][[1]])

# extract models
mod_crawley_add <- all_models[[1]][[2]]
mod_crawley_het <- all_models[[2]][[2]]

table_add <- DelTestVar(ophi_crawley_add)
table_het <- DelTestVar(ophi_crawley_het)
row.names(table_het)[2:nrow(table_het)] <- names(ophi_crawley_het[2:ncol(ophi_crawley_het)])

# same for AIC-----------------------------------------------------------------------------------
ophi <- ophi_add

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
row.names(vals)[2:nrow(vals)] <- modnames

