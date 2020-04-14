library(dplyr)
library(tidyr)
library("Hmisc")

df <- read.table("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/modeling/train_features.tsv",
                 header = TRUE,sep="\t")

df_x <- df %>% select(-label)#(-label,-sequence,-linker,-site1,-site2)

df_x2 <- df_x[ , !grepl( 'flankseq|ori|max|min' , names( df_x ) )]
cr <- cor(df_x2)
which(mat<=sort(cr), arr.ind = TRUE)
cor_df <- as.data.frame(as.table(cr))
sub_cor <- subset(cor_df, abs(Freq) > 0.5 & abs(Freq) < 1)
arr_cor <- sub_cor %>% arrange(desc(abs(Freq)))
# RCORR
write.csv(arr_cor,file="arrcortmp.csv")

df_x3 <-  df_x2[ , grepl( 'prot|roll|mgw|helt' , names( df_x2 ) )]
df_x3$label <- df$label

cmp <- df_x3 %>% 
  gather(key = variable, value = value, -label) %>%
  group_by(label, variable) %>%
  summarise(value = list(value)) %>%
  spread(label, value) %>%
  group_by(variable) %>%
  mutate(p_val_wilcox = wilcox.test(unlist(additive), unlist(cooperative))$p.value,
         p_val_t = t.test(unlist(additive), unlist(cooperative))$p.value,) %>%
  select(variable, p_val_wilcox, p_val_t) %>%
  arrange(p_val_wilcox,p_val_t)
write.csv(cmp,file="cmp.csv")

coop <- df_x3 %>% filter(label == "cooperative")
add <- df_x3 %>% filter(label == "additive")

x1 <- coop %>% filter(mgw_str_pos.2 > -500)
x2 <- add %>% filter(mgw_str_pos.2 > -500)
boxplot(x1$mgw_str_pos.2, x2$mgw_str_pos.2)
wilcox.test(coop$mgw_str_pos.2, add$mgw_str_pos.2, "less")$p.value

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


