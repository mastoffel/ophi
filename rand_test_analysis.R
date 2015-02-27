# IDEA: plot on a surface the relation of AIC and r2 of 10 000 randomized datasets.
# Then show that our model is outside the mass (aic = 2447, r2 = 0.59, adjr2 = 0.172)

library(ggvis)
library(dplyr)

# original
# loading data
all_aic <- unname(unlist(read.table("data/server/aic_10k.txt")))
all_r2 <- unname(unlist(read.table("data/server/r2_10k.txt")))
all_adjr2 <- unname(unlist(read.table("data/server/adjr2_10k.txt")))

# combine to df
results <- tbl_df(data.frame(aic = all_aic, r2 = all_r2, adjr2 = all_adjr2))

# what proportion is higher?
higher <- results[(results$aic < 2447) & (results$r2 > 0.59), ]
prop <- nrow(higher)/nrow(all_aic)

# ggvis
all_aic <- data.frame(aic = all_aic)
all_aic %>% ggvis(x = ~aic, fill := "#fff8dc", opacity := input_slider(0,1)) %>% 
        layer_histograms(width := input_slider(0,2,step=0.1, label = "width"), boundary = 0) %>%
        add_axis("x", title = "AIC") %>%
        add_axis("y", title = "frequency")


# getting grid for surface plot
# aic_grid <- cut(sort(results$aic), 1000)
# aic_grid <- seq(from = (min(results$aic)), to = (max(results$aic)), by = (max(results$aic)-min(results$aic))/1000)
# r2_grid <- seq(from = (min(results$r2)), to = (max(results$r2)), by = (max(results$r2)-min(results$r2))/1000)
# grid <- data.frame(aic = aic_grid, r2 = r2_grid)

# 2d density estimation
library(MASS)
# library("KernSmooth")

dens <- kde2d(results$aic, results$r2, n = 100)
# same with KernSmooth package
# dens_2 <- bkde2D(cbind(results$aic, results$r2), bandwidth = c(0.5,0.5),c(gridsize = c(100,100))
# terrain.colors(30)
surf.colors <- function(x, col = rev(gray(seq(0, 1, len = 40)))) {
        # First we drop the 'borders' and average the facet corners
        # we need (nx - 1)(ny - 1) facet colours!
        x.avg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] +
                          x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)]) / 4
        # Now we construct the actual colours matrix
        colors = col[cut(x.avg, breaks = length(col), include.lowest = T)]
        return(colors)
}

# border = NA
# col = surf.colors(dens$z)
res <- persp(dens$x, dens$y, dens$z, 
             col = "gray", #surf.colors(dens$z),
             phi = 30, theta = 240, box = T, # theta 240 phi 20
             shade = .6, ticktype = "detailed", xlab = "\n\n\naic", ylab = "\n\n\nr2", 
             zlab = "\n\n\ndensity", d=5, 
             nticks = 6, border = NA)

# arrows(trans3d(x = 2447, y = 0.59, z = c(0.06), pmat = res), trans3d(x = 2447, y = 0.59, z = c(0.02), pmat = res))
points(trans3d(x = 2447, y = 0.59, z = 0.03, pmat = res), col = "black", pch = 19, cex = 1.5)
lines(trans3d(x = 2447, y = 0.59, z = c(0,0.03), pmat = res),lwd = 2,lty = 5, col = "black")
lines(trans3d(x = 2447, y = 0.59, z = 0.001, pmat = res), col = "orange", pch = "0", cex = 1.5)

pdf("SavingExample.pdf", width=7, height=5)
filled.contour(dens,
               #color.palette = colorRampPalette(c('#f7f7f7', '#edf8e9', '#bae4b3','#74c476', '#31a354','#006d2c')),
               #color.palette = colorRampPalette(c('#f7f7f7', '#cccccc','#969696', '#636363','#252525')),
               # color.palette = colorRampPalette(c('white', '#a1dab4','#41b6c4', '#2c7fb8','#253494')),
               #color.palette = colorRampPalette(c('white', '#bdc9e1','#67a9cf', '#1c9099','#016c59')),
               #color.palette = colorRampPalette(c('white', '#bdbdbd','#969696','#737373','#525252', '#252525', '#000000')),
               #color.palette = colorRampPalette(c('white','lightblue', 'blue', 'red' ,'darkred')),
               # color.palette = colorRampPalette(c('white', '#dadaeb','#bcbddc', '#6a51a3', '#6a51a3' ,'#3f007d')),
               color.palette = colorRampPalette(c( 'white', '#d0d1e6','#74a9cf', '#0570b0','#023858')),
               xlim = c(2390, 2510),
               ylim = c(0.15, 0.7), 
               plot.axes     = { points(2447,0.59, pch = 21, col = "black", cex = 2, bg="red", lwd = 2)
                                 axis(1, seq(2380, 2520, by = 20))
                                 axis(2, seq(0.1, 0.8,   by = 0.1)) })
dev.off() 


library("scatterplot3d")
scatterplot3d(x = dens$x, y = dens$y, z = dens$z)
plot(dens)


plot(results$aic, 1-results$r2)
points(x = 2447, y = 0.41, col = "red")


library(rgl)
?hist3d
# get adjr2
ophi_code <- read.csv("data/ophi_code.csv", colClasses = c(rep("numeric",2),rep("factor",98)))
mod_code <- glm(day_falling ~., data = ophi_code)
# AIC backward stepwise 
library(MASS)
mod_aic_code <- stepAIC(mod_code, direction = "backward")
r2 <- ((summary(mod_aic_code)$null.deviance) - 
               (summary(mod_aic_code)$deviance)) / (summary(mod_aic_code)$null.deviance)
npred <- length(summary(mod_aic_code)$contrasts)
nobs <- nrow(ophi_code)
adjr2 <- 1 - ((1 - r2^2) * (nobs - 1) / (nobs - npred - 1))
