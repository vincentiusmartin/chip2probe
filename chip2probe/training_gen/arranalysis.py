
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_chamber_corr(dfx, dfy,
                namecol="Name", valcol="Alexa488Adjusted", repcols="rep", extrajoincols=[],
                title="",
                median=False, log=False,
                shownames=False):
    """
    Get correlation between

    Args:
        dfx: data frame from chamber 1
        dfy: data frame from chamber 2
        namecol: column name for the probe names
        valcol: column name for the array intensity
        extrajoincols: column name for the extras such as orientation
        title: plot title
        median: get median between replicates
        log: return result in log scale
        shownames: show names of the probe on the plot, useful for debugging
    Return:
        Show plot result for correlation between chamber
    """
    # correlation plot of negative controls in chamber 1 vs chamber 2

    joincols = [namecol, repcols] + extrajoincols
    dfcombined = dfx.merge(dfy, on=joincols)[joincols + ["%s_x"%valcol, "%s_y"%valcol]]

    if median:
        dfcombined = dfcombined.groupby([namecol] + extrajoincols)[["%s_x"%valcol, "%s_y"%valcol]].median().reset_index()
    if log:
        dfcombined["%s_x"%valcol] = np.log(dfcombined["%s_x"%valcol])
        dfcombined["%s_y"%valcol] = np.log(dfcombined["%s_y"%valcol])

    x = dfcombined["%s_x"%valcol].values
    y = dfcombined["%s_y"%valcol].values

    slope, intercept = np.polyfit(x, y, 1)
    # Create a list of values in the best fit line
    abline_values = [slope * i + intercept for i in x]
    # plot diagonal
    plt.plot([min(min(x),min(y)),max(max(x),max(y))], [min(min(x),min(y)),max(max(x),max(y))], color='black', label='diagonal', linewidth=0.5)
    # plot best fit line
    plt.plot(x, abline_values, color='red', label='best fit', linewidth=0.5)
    r_squared = np.corrcoef(x,y)[0,1]**2
    plt.scatter(x, y, s=2)

    if shownames:
        names = dfcombined[namecol].values
        for i in range(len(names)):
            plt.text(x[i], y[i], names[i], fontsize=3)

    plt.xlabel('Chamber 1')
    plt.ylabel('Chamber 2')
    plt.title(title)
    plt.text(min(x), max(y), 'R² = %0.2f' % r_squared, color='red')
    plt.text(min(x), 0.9*max(y), 'best fit: y=%0.2fx + %0.2f' % (slope,intercept), color='red')
    plt.legend(loc='lower right')
    plt.show()

def plot_classified_labels(df):
    # HARDCODED - FIX
    # cooperative
    x = np.log(df[df['label']=='cooperative']['Alexa488Adjusted_x'].values)
    y = np.log(df[df['label']=='cooperative']['Alexa488Adjusted_y'].values)
    plt.scatter(x, y, color='blue', s=3)

    #additive
    x = np.log(df[df['label']=='additive']['Alexa488Adjusted_x'].values)
    y = np.log(df[df['label']=='additive']['Alexa488Adjusted_y'].values)
    plt.scatter(x, y, color='gray', s=3)

    # anti-coop
    x = np.log(df[df['label']=='anticoop']['Alexa488Adjusted_x'].values)
    y = np.log(df[df['label']=='anticoop']['Alexa488Adjusted_y'].values)
    plt.scatter(x, y, color='red', s=3)

    plt.xlabel('Chamber 1')
    plt.ylabel('Chamber 2')
    x = np.log(df['Alexa488Adjusted_x'].values)
    y = np.log(df['Alexa488Adjusted_y'].values)
    plt.plot([min(min(x),min(y)),max(max(x),max(y))], [min(min(x),min(y)),max(max(x),max(y))], color='black')
    plt.title("")
    plt.show()
