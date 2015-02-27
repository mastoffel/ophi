full_table <- function(df_mod, df_ref) {
        
        vals <- DelTestVar(df_mod)
        
        # reduce data frame by intercept 
        val_red <- vals[2:nrow(vals), ]
        
        # extract F and p value
        val_red2 <- (val_red[, c("F","P (F-test)")])
        
        # sort within all SNP´s
        final_df <- data.frame("F" = rep(NA, 98), "pval" = rep(NA, 98),
                               row.names = names(df_ref))
        
        final_df[row.names(val_red2), ] <- val_red2
        
        final_df

}