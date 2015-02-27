# idea for plot: allele frequency change with increasing number of deaths

ophi_code <- (read.csv("data/ophi_code.csv", colClasses = c(rep("numeric",100))))
ophi <- ophi_code[,-2]
ophi_aic_mod <- read.csv("data/df_aic_code.csv", colClasses = c(rep("numeric",50), row.names = TRUE))
ophi <- ophi[names(ophi_aic_mod)]

# calculate which allel per snp has the minor frequency
maf_ident <- function(snp, allel = -99) {
        all_alleles <- length(snp) * 2
        allel_freq <- function(snp_elem, allel) {
               count <- 0
               if (snp_elem == allel) count <- 2
               if (snp_elem == 2) count <- 1
               count
               }
        freq_1 <- sum(sapply(snp, allel_freq, allel = 1)) / all_alleles
        freq_2 <- sum(sapply(snp, allel_freq, allel = 3)) / all_alleles
        
        if (freq_1 > freq_2) {
                return(3)
        } else {
                return(1)
        }
}

# calculate the frequency of a given allel per snp
maf_calc <- function(snp, allel) {
        all_alleles <- length(snp) * 2
        allel_freq <- function(snp_elem, allel) {
                count <- 0
                if (snp_elem == allel) count <- 2
                if (snp_elem == 2) count <- 1
                count
        }
        freq <- sum(sapply(snp, allel_freq, allel)) / all_alleles
}

# identifies all maf alleles
all_maf_alleles <- apply(ophi[, 2:ncol(ophi)], 2, maf_ident)
# gets all frequencies of maf-alleles
all_maf <- lapply(1:(ncol(ophi)-1), function(x) maf <- maf_calc(snp = ophi[, x+1], allel = all_maf_alleles[x]))
all_maf <- unlist(all_maf)

# try out ggvis
maf_df <- data.frame(maf = all_maf)

maf_df %>%
        ggvis(~maf) %>%
        layer_densities(
                adjust = input_slider(.1, 2, value = 1, step = .1, label = "Bandwidth adjustment")
        )

# prep
obs <- nrow(ophi)
var <- ncol(ophi)

# calculate minor allele frequencies per quartile 
maf_0_25 <- unlist(lapply(1:(var-1), 
                            function(x) maf <- maf_calc(snp = ophi[1:(1 * obs/4), x+1], 
                                                        allel = all_maf_alleles[x])))
maf_25_50 <- unlist(lapply(1:(var-1), 
                          function(x) maf <- maf_calc(snp = ophi[(1 * obs/4 + 1):(2 * obs/4), x+1], 
                                                      allel = all_maf_alleles[x])))
maf_50_75 <- unlist(lapply(1:(var-1), 
                           function(x) maf <- maf_calc(snp = ophi[(2 * obs/4 + 1):(3 * obs/4), x+1], 
                                                       allel = all_maf_alleles[x])))
maf_75_100 <- unlist(lapply(1:(var-1), 
                           function(x) maf <- maf_calc(snp = ophi[(3 * obs/4 + 1):(4 * obs/4), x+1], 
                                                       allel = all_maf_alleles[x])))

maf_0_75 <- unlist(lapply(1:(var-1), 
                          function(x) maf <- maf_calc(snp = ophi[1:(3 * obs/4), x+1], 
                                                      allel = all_maf_alleles[x])))

plot(density(all_maf, adjust = 1))
lines(density(maf_0_25, adjust = 1), col = "blue")
lines(density(maf_25_50, adjust = 1), col = "red")
lines(density(maf_50_75, adjust = 1), col = "green")
lines(density(maf_75_100, adjust = 1), col = "purple")

lines(density(maf_0_75, adjust = 1), col = "green")

# KS test
ks.test(maf_0_25, maf_75_100)
