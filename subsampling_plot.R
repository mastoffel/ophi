# non-parametric bootstrapping on the full data set with subsequent stepAIC
library(dplyr)
# load first 10 000
all_names <- tbl_df(read.table("data/server/all_names_08.txt",  stringsAsFactors = FALSE))
# second 10 000
all_names_2 <- tbl_df(read.table("data/server/all_names_08_2.txt", stringsAsFactors = FALSE))
all_names <- rbind(all_names, all_names_2)

# create vector
names_vec <- as.character(unlist(all_names))
# sort
tab_names <- sort(table(names_vec), decreasing = TRUE)
tab_names <- tab_names[2:length(tab_names)]
# delete day_falling
snps <- names(tab_names)
# logical for occuring in model
in_mod <- snps %in% names(df_aic_code)
# data frame
snp_df <- data.frame(names = (1:length(tab_names)), 
                     freq = unname(tab_names)/20000, row.names = 1:length(tab_names),
                     mod = as.numeric(in_mod),
                     stringsAsFactors = FALSE)
snp_df$mod <- factor(snp_df$mod, levels = rev(levels(factor(snp_df$mod))))

# theme
library(ggplot2)
library(grid)
theme_dens <- theme_minimal() +
        theme(strip.text.x = element_text(vjust=1,size = 18),
              #strip.background = element_rect(fill="white",linetype="blank"),
              axis.title.x = element_text(vjust= -2,size = 16),
              axis.title.y = element_text(vjust= 3,size = 16),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              plot.margin = unit(c(2,2,2,2), "cm"), # less margin on down side
              # axis.line = element_line(size=1),
              # panel.border = element_blank() ,
              panel.grid.major = element_blank(),
              # axis.line = element_line(color = 'black'),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour="black",fill=NA,size=1),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12, face = "plain"))

ggplot(data = snp_df, aes(x = names, y = freq, fill = factor(mod))) +
        geom_bar(color = "black", size = 0.1, width=1, stat="identity") + #(colour="white"
        theme_dens +
        scale_x_continuous(limits = c(0.5, 100),
                           expand = c(0, 0), 
                           breaks = c(seq(from = 0, to = 100, by = 10))) +
        scale_y_continuous(limits = c(0, 1),
                           expand = c(0, 0),
                           breaks = c(seq(from = 0, to = 1, by = 0.2))) +
        scale_fill_manual(name="Final Model",
                          values = c("black", "#f7f7f7"), ##636363
                          labels=c("Included", "Excluded")) +
        #theme(panel.border = element_blank(),
              #axis.line = element_line(color = 'black')) +
        theme(legend.position=c(0.85, 0.8)) +
                           # breaks = c(seq(from = 0, to = 0.05, by = 0.01))) +
        xlab("SNP") +
        ylab("Frequency")
ggsave(filename = "Subsamp.jpg", width = 7, height = 5, units = "in")  
ggsave(filename = "Subsamp.pdf", width = 7, height = 5, units = "in")       


library(ggvis)
snp_df %>% 
        ggvis(~names, ~freq) %>%
        layer_bars(fill = mod) %>%
        layer_bars(~freq[names == 1])

library(ggplot2)
qplot(x = names, y = freq, data=snp_df, geom="bar", fill=factor(mod))
