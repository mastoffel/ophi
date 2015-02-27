# allel sharing analysis (as)

ophi_code_snp <- ophi_code[-c(1,2)]

ophi_as <- as.data.frame(matrix(nrow = nrow(ophi_code), ncol = nrow(ophi_code)))
row.names(ophi_as) <- 1:nrow(ophi_as)
names(ophi_as) <- 1:nrow(ophi_as)
str(ophi_as)

as_pair <- function(ind1, ind2) {
        df <- rbind(ind1, ind2)
        calc_points <- function(df) {
                if(df[1] == df[2]) {
                        point <- 1
                } else if (df[1] == 2 || df[2] == 2) {
                        point <- 0.5
                } else {
                        point <- 0
                }
        }
                        
        shared_sum <- sum(apply(df, 2, calc_points))
}

for (i in seq(1:288)) {
        for (k in seq(1:288)) {
                ophi_as[i, k] <- as_pair(ophi_code_snp[i, ], ophi_code_snp[k, ])
        }
}

ophi_as_rel <- ophi_as/98
ophi_as_rel[upper.tri(ophi_as_rel)] <- NA

288/4
