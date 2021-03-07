#path <- "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/ch1_ch2/coop_ch1_vs_ch2/tables/both_ori_plt.csv"
path <- "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/homotypic/training/lbled_o1_selected.csv"
df <- read.csv(path)

plot(log(df$indiv_median), log(df$two_median), cex=0.3)
identify(log(df$indiv_median), log(df$two_median))
log(df[972,]$Alexa488Adjusted_x)
df[1877,]
