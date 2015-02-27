
# based on probability ----------------------------------------------------------------------------
rand_aic <- function(...) {
        
        ophi_code <- ophi_code
        # randomization of allels according to probabilities
        rand_allelfreq <- function(x) {
                # x is a snp vector
                df <- as.data.frame(table(x)) 
                df$prob  <-  df$Freq / sum(df$Freq)
                x_rand <- sample(df$x, length(x), replace = TRUE, prob = df$prob)      
        }
        # create new df
        ophi_code_temp <- ophi_code
        ophi_code_temp[3:ncol(ophi_code_temp)] <- as.data.frame(apply(ophi_code[3:ncol(ophi_code)],2, rand_allelfreq))
        ophi_code_temp
        
        # modelfit code
        mod_code <- glm(day_falling ~., data = ophi_code_temp)
        
        # AIC backward stepwise 
        mod_aic_code <- stepAIC(mod_code, direction = "backward")
        
        # get final aic
        final_aic <- summary(mod_aic_code)$aic
}