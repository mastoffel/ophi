# plotting

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


# main plot
ophi_code_2 <- ophi_code
ophi_code_2$day_falling <- ophi_code$day_falling + 71
ggplot(ophi_code_2, aes(x = day_falling)) +
        geom_histogram(color = "black", binwidth = 5, 
                       aes(y = ..density.., fill = ..count..)) + # fill = ..count..
        # scale_fill_brewer()
        scale_fill_gradient("Counts", low="#f7f7f7", high="#636363") +
        geom_density(size = 1, alpha=.2, adjust = 1) +
        # scale_fill_gradient2("Counts", low = "grey46", high = "grey33") +
        theme_dens +
        theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black')) +
        #legend.key = element_rect(colour = 'black'))+
        scale_x_continuous(limits = c(65, 160),
                           expand = c(0, 0),
                           breaks = c(seq(from = 70, to = 160, by = 10))) +
        scale_y_continuous(limits = c(0, 0.05),
                           expand = c(0, 0),
                          breaks = c(seq(from = 0, to = 0.05, by = 0.01))) +
        xlab("Survival day") +
        ylab("Density") 
  
ggsave(filename = "Survival_day.jpg", width = 7, height = 5, units = "in")  
ggsave(filename = "Survival_day.pdf", width = 7, height = 5, units = "in") 



# density plot 
source("multiplot.R")

bpl <- ggplot(ophi_code, aes(x = as.numeric(row.names(ophi_code)), y = day_falling)) +
        geom_boxplot(colour = I("#3366FF"), size = 1.3) +
        coord_flip() +
        theme_dens_2 +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())

d <- density(ophi_code$day_falling)     
DF <- data.frame(x = d$x, y = d$y)
gpl <- ggplot(DF, aes(x=x, y=y)) +
        geom_segment(aes(xend=x, yend=0, colour=abs(y)^0.7*sign(y),
                     lineend = "butt")) +
        geom_line() +
        # scale_colour_gradient(low="red", high="blue") +
        theme_dens_2 +
        xlab("Day falling") +
        ylab("Density") +
        theme(legend.position="none") 

multiplot(bpl, gpl, ncol = 1)


# plotting with ggplot (not really working)

plot_df <- data.frame(r = summary_out$adjr2, levels = 1:length(summary_out$adjr2), snps = number_snp)

p1 <- ggplot(data = plot_df, aes(x = levels, y = r)) +
        geom_point() +
        theme_dens +
        theme(panel.border = element_blank(),
              axis.line = element_line(color = 'black'))+
        scale_x_continuous(limits = c(0, 200),
                           expand = c(0, 0),
                           breaks = c(seq(from = 0, to = 200, by = 20))) +
        scale_y_continuous(limits = c(0, 0.5),
                           expand = c(0, 0),
                           breaks = c(seq(from = 0, to = 0.5, by = 0.1))) +
        xlab("Factor Levels") +
        ylab(expression("Adj." ~ italic(R)^{2}))

p2 <- ggplot(data = plot_df, aes(x = snps, y = r)) +
        geom_point() +
        theme_dens +
        theme(panel.border = element_blank(),
              axis.line = element_line(color = 'black'))+
        scale_x_continuous(limits = c(0, 100),
                           expand = c(0, 0),
                           breaks = c(seq(from = 0, to = 100, by = 10))) +
        scale_y_continuous(limits = c(0, 0.5),
                           expand = c(0, 0),
                           breaks = c(seq(from = 0, to = 0.5, by = 0.1))) +
        xlab("SNPs") +
        ylab(expression("Adj." ~ italic(R)^{2}))

