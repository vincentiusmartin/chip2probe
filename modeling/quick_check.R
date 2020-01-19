library(dplyr)
library("Hmisc")

df <- read.table("/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/modeling/all_features_ori_one_hot.tsv",
                 header = TRUE,sep="\t")

df_x <- df %>% select(-label)#(-label,-sequence,-linker,-site1,-site2)
#cor(df_x)
rcorr(as.matrix(df_x))

coop <- df %>% filter(label == "cooperative")
add <- df %>% filter(label == "additive")
anti <- df %>% filter(label == "anticoop")

cor.test(df$A, df$G, method=c("pearson"))

shapiro.test(coop$A)$p.value
shapiro.test(coop$C)$p.value
shapiro.test(coop$G)$p.value
shapiro.test(coop$T)$p.value

t.test(add$GC_content, coop$GC_content, "greater")$p.value
wilcox.test(add$C, coop$C, "greater")$p.value
wilcox.test(add$G, coop$G, "greater")$p.value
wilcox.test(add$T, coop$T, "greater")$p.value


