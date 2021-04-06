#path <- "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/ch1_ch2/coop_ch1_vs_ch2/tables/both_ori_plt.csv"
path <- "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/main_nar/output/Runx1Ets1/label_pr/both_ori_plt_runx1_ets1.csv"
df <- read.csv(path)

axis(side=1)
plot(log(df$indiv_median), log(df$two_median), cex=0.3)
identify(log(df$indiv_median), log(df$two_median))
df[100,]
