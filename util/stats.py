import rpy2.robjects as robjects

wilcox_r = robjects.r['wilcox.test']
shapiro_r = robjects.r['shapiro.test']
numeric = robjects.r['as.numeric']
t_test = robjects.r['t.test']

def wilcox(x,y,alternative="two.sided"):
    """
    Return wilcox of x and y. From https://www.rdocumentation.org/packages/stats/versions/3.6.1/topics/wilcox.test
    :param x,y
    :param alternative: "two.sided", "less", "greater"
    """
    x_num = numeric(x)
    y_num = numeric(y)
    return wilcox_r(x_num,y_num,alternative=alternative).rx("p.value")[0][0]
