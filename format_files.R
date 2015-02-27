# transform to ped format
library(stringr)
raw <- read.csv("raw.csv", row.names = 1, header = TRUE)[1:288, ]
snp_names <- read.csv("snp_names.csv", header = FALSE)
row.names(raw) <- str_replace_all(row.names(raw), "Ophi ", "")
row.names(raw)[14] <- 14
nobs <- length(row.names(raw))
nsnps <- length(names(raw))/2
# PED file
ped_temp <- data.frame("FID" = as.numeric(row.names(raw)), "IID" = as.numeric(row.names(raw)), "PatID" = rep(0, nobs),
                     "MatID" = rep(0, nobs), "Sex" = rep(0, nobs), "Phenotype" = ophi_code$day_falling)
ped <- cbind(ped_temp, raw)
write.table(ped, "ophi.ped", row.names = FALSE, col.names = FALSE)

# MAP file
map <- data.frame("Chr" = rep(0, nsnps), "FID" = snp_names, "gendist" = rep(0, nsnps), "snppos" = rep(9999, nsnps) )
write.table(map, "ophi.map", row.names = FALSE)

#phenotype
ophi_phen <- data.frame("FID" = as.numeric(row.names(raw)), "IID" = as.numeric(row.names(raw)), "Phenotype" = ophi_code$day_falling)
write.table(ophi_phen, "ophi.phen", row.names =FALSE, col.names = FALSE)

# canditate vs. non-candidate

ophi_code_can <- ophi_code[1:65]
ophi_code_nocan <- ophi_code[c(1, 66:100)]

mod1 <- glm(day_falling ~ ., data = ophi_code_can)
mod2 <- glm(day_falling ~ ., data = ophi_code_nocan)
summary(mod1)
summary(mod2)
summary(anova(mod1, mod2))


