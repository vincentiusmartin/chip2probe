library(dplyr)
library(ggplot2)

setwd("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/main_nar")

df <- read.csv("output/Runx1Ets1/label_pr/both_ori_plt_runx1_ets1.csv")
df$log_indiv = log(df$intensity_x)
df$log_two = log(df$intensity_y)
plot(df$log_indiv , df$log_two, cex=0.5)
identify(df$log_indiv , df$log_two)
df[1060,]


# df <- read.csv("output/Ets1Runx1/label_pr/both_ori_plt_ets1_runx1.csv") %>%
#     rename("Ets1"=intensity_x, "Ets1+Runx1"=intensity_y) %>%
#     select(c("Name","Ets1","Ets1+Runx1")) %>%
#     pivot_longer(., cols = c("Ets1","Ets1+Runx1"), names_to = "Experiment", valu  es_to = "intensity")

# df <- read.csv("output/Runx1Ets1/label_pr/both_ori_plt_runx1_ets1.csv") %>%
#     rename("Runx1"=intensity_x, "Runx1+Ets1"=intensity_y) %>%
#     select(c("Name","Runx1","Runx1+Ets1")) %>%
#     pivot_longer(., cols = c("Runx1","Runx1+Ets1"), names_to = "Experiment", values_to = "intensity")

df <- read.csv("output/Ets1Ets1/label_pr/m1m2m3wt.csv") %>%
  select(c("Sequence","m3","m1","m2","wt")) %>%
  pivot_longer(., cols = c("m3","m1","m2","wt"), names_to = "#sites", values_to = "intensity")
df$`#sites`[df$`#sites` == "m1"] <- "one_site (m1/m2)"
df$`#sites`[df$`#sites` == "m2"] <- "one_site (m1/m2)"
df$`#sites`[df$`#sites` == "m3"] <- "no_site (m3)"
df$`#sites`[df$`#sites` == "wt"] <- "two_sites (wt)"

# x=Experiment/`#sites`, y=intensity) red ["#b22222","#FFA07A"] blue #75bbfd","#0343df
p <- ggplot(df, aes(x=`#sites`, y=intensity)) + 
  geom_violin(aes(fill= `#sites`),scale=3) + 
  geom_boxplot(width=0.1) +
  #scale_fill_manual(values=c("#b22222","#FFA07A")) + 
  theme(legend.position = "none") +
  scale_y_continuous(trans = 'log', breaks = c(250,1000,3000,10000,25000,60000, 150000))
ggsave("er_signal.pdf", width=8) # 5.92 x 6.04



#=========

df <- read.csv("output/Ets1Runx1/training/train_ets1_runx1.tsv",sep="\t")
nrow(df[(df$distance == 5) & (df$label == "cooperative"),])
dist <- df[,(c("distance","label"))]
fisher.test(dist)

# 336/528
# 995/2202

dfdist <- df[df$distance == 11,]
ndist <- nrow(dfdist)
dist_coop <- nrow(dfdist[dfdist["label"] == "cooperative",]) 
dist_not_coop <- ndist - dist_coop
notdist_coop <- nrow(df[df["label"] == "cooperative",]) - dist_coop
notdist_not_coop <- nrow(df[df["label"] != "cooperative",]) - dist_not_coop
# dist == 5
dat <- data.frame(
  "dist" = c(dist_coop, notdist_coop),
  "not_dist" = c(dist_not_coop, notdist_not_coop),
  row.names = c("coop", "not_coop"),
  stringsAsFactors = FALSE
)
dat
chisq.test(dat)$expected
fisher.test(dat)

wilcox.test(c(0.82,
              0.859,
              0.872,
              0.854,
              0.871,
              0.822,
              0.823,
              0.861,
              0.866,
              0.87),
            c(0.818,
              0.862,
              0.869,
              0.858,
              0.872,
              0.812,
              0.824,
              0.858,
              0.861,
              0.865),
            'greater')



