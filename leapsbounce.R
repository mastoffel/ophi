# regression subset selection with leaps and bounce algorithm-------------------------

## with code #
library(leaps)
regsubsets_out <- regsubsets(day_falling ~. , data = ophi_code,
                             nbest = 1, method = "backward", nvmax = NULL)
plot(regsubsets_out, scale = "adjr2", main = "Adjusted R^2")

# best model at each variable number
summary_out <- summary(regsubsets_out, matrix.logical=TRUE)
models_df <- as.data.frame(summary_out$outmat)

# exchange factor level names with variable names
library(stringr)
shorten <- function(x) {
        str_sub(x, 1, str_length(x) - 1)
}
original_names <- unlist(lapply(names(models_df), shorten))

# get indices of variables for all models
get_ind <- function(x){
        which(as.logical(x))
}
all_mod_ind <- apply(models_df, 1, get_ind)

number_snp <- vector()
for (i in 1:length(all_mod_ind)) {
        number_snp[i] <- length(unique(original_names[all_mod_ind[[i]]]))
}

# plot --------------------------------------------------------------------------
# number of snps per model
library(scales)
pdf("leapsbounce.pdf", width=7, height=5)

par(mar = c(6,6,5,5))
plot(summary_out$adjr2, 
     bty = "n", axes = F, xlab = "", ylab = "", pch=16, 
     ylim = c(0, 0.55))
axis(2, at=c(seq(from=0,to=0.5,by=0.1)),las=1, cex.axis = 1.2)
mtext(expression("Adjusted" ~ italic(r)^{2}), 2,line= 4, cex = 1.2)
axis(1, at=c(seq(1,181, by = 20),las=1), cex.axis = 1.2)
mtext("Contrasts", 1,line= 1 ,at=-22, cex = 1.2)


# axis(1,0:180,line=1,col="red",col.ticks="red",col.axis="red")
axis(1,at=c(seq(from=1,to=181,by=10)),
     labels=number_snp[seq(from = 1, to = 181, by = 10)],
     line=2.5, cex.axis = 1.2)
mtext("SNPs", 1,line= 3.5 ,at=-16, cex = 1.2)
# segments, arrows
arrows(89,summary_out$adjr2[89],89,0, lty = 1, lwd = 2, col = "black", length = 0.15)
points(89, summary_out$adjr2[89], col = "red", pch = 16, cex = 1)
# abline(v = 89, lty = 5, lwd = 2) # 89 factor levels
# text(132, 0.485,"Lowest RSS model \n(89 factor level comparisons,\n66 SNP´s)", col = "blue")
text(89, 0.525,"Lowest RSS model \n(66 SNPs, 89 Contrasts)", col = "black", cex = 1.2)

dev.off()
  

# ggvis
library(ggvis)
df_lb <- data.frame(lev = 1:length(summary_out$adjr2), snp = number_snp,
                    adjr2 = summary_out$adjr2)

df_lb %>% 
        ggvis(x = ~lev, y = ~adjr2) %>%
        layer_points() %>%
        add_axis("x") %>%
        add_axis("x", offset = 50, grid = FALSE, labels = listsnp)



library(broom)
test <- tidy(mod_aic_code)





# percentage of snps from aic model in leaps bounce best model
# snps in best leaps bounce model
snp_lb_code <- unique(original_names[unlist(all_mod_ind[89])])
# snps in AIC
snp_aic_code <- names(df_aic_code)[-1]

# percentage of snps from best AIC model in best lb model
sum(snp_aic_code %in% snp_lb_code)/49

which.max(summary_out$adjr2)
best_lp <- as.data.frame(summary_out$which[88,])
names(best_lp) <- "in"
best_lp <- best_lp[best_lp[1] == TRUE]

plot(summary.out$rss)
plot(summary.out$rsq)

## with additive #
regsubsets_out <- regsubsets(day_falling ~. , data = ophi_add,
                             nbest = 1, method = "backward", nvmax = NULL)
regsubsets_out
plot(regsubsets_out, scale = "adjr2", main = "Adjusted R^2")

# best model at each variable number
summary.out <- summary(regsubsets_out, matrix.logical=TRUE)
as.data.frame(summary.out$outmat)
plot(summary.out$adjr2)
plot(summary.out$rss)
plot(summary.out$rsq)
