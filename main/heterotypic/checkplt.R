path <- "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/ch1_ch2/coop_ch1_vs_ch2/tables/both_ori_plt.csv"
df <- read.csv(path)

plot(log(df$Alexa488Adjusted_x), log(df$Alexa488Adjusted_y), cex=0.3)
identify(log(df$Alexa488Adjusted_x), log(df$Alexa488Adjusted_y))
log(df[972,]$Alexa488Adjusted_x)
df[1864,]
