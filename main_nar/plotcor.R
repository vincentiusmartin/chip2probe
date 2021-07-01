library(dplyr)
library(ggplot2)

setwd("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/main_nar")


df <- read.csv("output/Ets1Runx1/label_pr/med_m1m2.csv")
#df <- read.csv("output/Runx1Ets1/label_pr/med_m1m2.csv")
df.lm <- lm(intensity_y ~ intensity_x, df)
rsq <- summary(df.lm)$r.squared
df.lm
rsq

df <- rename(df, "Ets1"=intensity_x, "Ets1+Runx1"=intensity_y)

# geom_abline(slope = coef(df.lm)[[2]], intercept = coef(df.lm)[[1]], color="red", linetype = "dashed") 
p <- ggplot(df, aes(x=Ets1, y=`Ets1+Runx1`)) +
    geom_point(color='blue', size = 0.8) +
    scale_x_continuous(trans = 'log',  breaks = c(250,1000,3000,10000,25000,60000, 150000)) +
    scale_y_continuous(trans = 'log',  breaks = c(250,1000,3000,10000,25000,60000, 150000)) +
    geom_abline(color="black", size=0.8) +
    geom_smooth(method='lm', formula= y~x, color="red", linetype = "dashed") +
    theme_bw()+ 
    theme(text = element_text(size=22), axis.text.x = element_text(angle = 30, hjust = 1)) +
    ggsave("er_signal.pdf", width=8) 

