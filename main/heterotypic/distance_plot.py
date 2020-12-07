import seaborn as sns
import matplotlib.pyplot as plt
from  matplotlib.ticker import FuncFormatter

dd = [int(x) for x in pd.read_csv("peak_dist.csv")["dist"].tolist()]
plt.figure(figsize=(2,4))
plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: int(y)))
sns.boxplot(data=dd, width=.2, color="gray")
plt.title("")
plt.ylabel("Distance")
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.tight_layout(pad=0, w_pad=0)
plt.savefig("a.png")
