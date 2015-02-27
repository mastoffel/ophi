all <- read.csv("gene_names_full.csv")
nonadd <- read.csv("nonadd_names.csv")
add <- read.csv("add_names.csv")
code <- read.csv("code_names.csv")

nonadd_names <- nonadd$variable
nonadd$variable <- str_replace(nonadd_names, "het", "code")

nonadd_c <- read.csv("nonadd_c_names.csv")
add_c <- read.csv("add_c_names.csv")
code_c <- read.csv("code_c_names.csv")

nonaddc_names <- nonadd_c$variable
nonadd_c$variable <- str_replace(nonaddc_names, "het", "code")


add_out <- all[!(all$variable %in% add$variable), ] 
nonadd_out <- all[!(all$variable %in% nonadd$variable), ] 
code_out <- all[!(all$variable %in% code$variable), ] 

write.csv(add_out, "add_aic.csv")
write.csv(nonadd_out, "nonadd_aic.csv")
write.csv(code_out, "code_aic.csv")

add_out_c <- all[!(all$variable %in% add_c$variable), ] 
nonadd_out_c <- all[!(all$variable %in% nonadd_c$variable), ] 
code_out_c <- all[!(all$variable %in% code_c$variable), ] 

write.csv(add_out_c, "add_c.csv")
write.csv(nonadd_out_c, "nonadd_c.csv")
write.csv(code_out_c, "code_c.csv")
