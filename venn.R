# venn diagram different approaches
library(stringr)
snps_code <- names(df_aic_code)[-1]
snps_add <- names(df_aic_add)[-1]
snps_het <- names(df_aic_het)[-1]
snps_addhet <- unique(c(snps_add, snps_het))
snps_subsamp <- snps[1:50]
snps_leaps <- snp_lb_code

snps_add <- str_replace_all(snps_add, "add", "code")
snps_het <- str_replace_all(snps_het, "het", "code")

# example
length(intersect(snps_code, snps_add))


library(VennDiagram)
draw.triple.venn((area1 = length(snps_add)),
                 area2 = length(snps_het),
                 area3 = length(snps_code),
                 n12 = length(intersect(snps_add, snps_het)),
                 n23 = length(intersect(snps_code, snps_het)),
                 n13 = length(intersect(snps_code, snps_add)),
                 n123 = length(intersect(snps_code,intersect(snps_add, snps_het))),
                 category = c("additive", "non-additive", "Final model"),
                 fill = c("red", "green", "blue"),
                 lty = 1,
                 cex = 2,
                 cat.cex = 2,
                 cat.col = c("red", "green", "blue"),
                 margin = 0.13
)
# "dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"

par(mar = c(3,3,3,3))
draw.quad.venn(area1 = length(snps_code),
               area2 = length(snps_addhet),
               area3 = length(snps_subsamp),
               area4 = length(snps_leaps),
               n12 = length(intersect(snps_code, snps_addhet)),
               n13 = length(intersect(snps_code, snps_subsamp)),
               n14 = length(intersect(snps_code, snps_leaps)),
               n23 = length(intersect(snps_addhet, snps_subsamp)),
               n24 = length(intersect(snps_addhet, snps_leaps)),
               n34 = length(intersect(snps_subsamp, snps_leaps)),
               n123 = length(intersect(snps_code,intersect(snps_addhet, snps_subsamp))),
               n124 = length(intersect(snps_code,intersect(snps_addhet, snps_leaps))),
               n134 = length(intersect(snps_code,intersect(snps_subsamp, snps_leaps))),
               n234 = length(intersect(snps_addhet,intersect(snps_subsamp, snps_leaps))),
               n1234 = length(intersect(snps_code,
                       (intersect(snps_addhet, intersect(snps_subsamp, snps_leaps))))),
                       category = c("Final", "Add & non-add", "Sub-sampling", "Leaps & bounds"),
                       fill = c("orange", "red", "green", "blue"),
                       lty = 1,
                       cex = 2,
                       cat.cex = 2,
                       cat.col = c("orange", "red", "green", "blue"),
               margin = 0.08
               )
               
               
draw.quintuple.venn(area1 = 49,
                    area2 = 27,
                    area3 = 36,
                    area4 = 50,
                    area5 = 66,
                    n12   = 19,
                    n13   = 24,
                    n14   = 42,
                    n15   = 47,
                    n23   = 15,
                    n24   = 19,
                    n25   = 21,
                    n34   = 27,
                    n35   = 31,
                    n45   = 49,
                    n123  = 13,
                    n124  = 18,
                    n125  = 19,
                    n134  = 22,
                    n135  = 23,
                    n145  = 41,
                    n234  = 13,
                    n235  = 14,
                    n245  = 19,
                    n345  = 27,
                    n1234 = 13,
                    n1235 = 13,
                    n1245 = 18,
                    n1345 = 22,
                    n2345 = 13,
                    n12345 = 13,
                    category = c("Final", "Add", "Non-add", "Subsampled", "LeapsBounds"),
                    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                    cat.cex = 2,
                    margin = 0.05,
                    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
                            1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
                    ind = TRUE
)


                    