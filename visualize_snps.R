# visualizing the snps

source("multiplot.R")
library(ggplot2)
library(grid)

names_code_org <- read.csv("names_code_org.csv",header = FALSE, stringsAsFactors = FALSE)[1]
names(df_aic_code)[2:50] <- names_code_org[,1]
get_all_plots <- function(model_df) {     
        all_plots <- list()
        for(i in 2:ncol(model_df)) {       
        snp <- names(model_df)[i]
        
        p <- ggplot(model_df, aes_string(x = snp, y = "day_falling")) +
                theme_bw(base_size = 10) +
                geom_boxplot() +
                           # fill = NA, shape = 21) +
                stat_summary(fun.y = "mean", geom = "text", 
                             label="----", size= 5, color= "blue")  +
                geom_point(position = position_jitter(width = 0.2), 
                           alpha = 0.7, size = 0.5, color = "orange") +
                theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
                      axis.title.y = element_blank())
                
                #stat_summary(fun.y=mean, colour="darkred", geom="point", 
                           #shape=18, size=3,show_guide = FALSE)
        nam <- paste("plot_", i, sep = "")
        # assign(nam, p)
        all_plots[[i-1]] <- p
        }
        return(all_plots)
}

# transform dfs for het and add to categorical
df_add <- df_aic_add
df_add[2:ncol(df_add)] <- lapply(df_aic_add[2:ncol(df_add)], as.factor)
df_het <- df_aic_het
df_het[2:ncol(df_het)] <- lapply(df_aic_het[2:ncol(df_het)], as.factor)


all_plots <- get_all_plots(df_aic_code)

multiplot(plotlist = all_plots[1:49], cols = 7)

ggsave(filename = "Boxplots.jpg", width = 7, height = 5, units = "in")  
ggsave(filename = "Boxplots.pdf", width = 7, height = 5, units = "in")  

multiplot(plotlist = all_plots[(length(all_plots)/2+1):length(all_plots)], cols = 4)