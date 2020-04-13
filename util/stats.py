"""This file contains functions to perfrom hypothesis testing using Python."""
from scipy.stats import mannwhitneyu, ttest_ind
import numpy as np


def wilcox(x, y, alternative="two-sided"):
    """
    Return p-value using mannwhitney u test of x and y.

    :param x,y
    :param alternative: "two.sided", "less", "greater"
    """
    x = np.array(x).reshape((-1, 1)).flatten()
    y = np.array(y).reshape((-1, 1)).flatten()
    statistic, p_value = mannwhitneyu(x, y, alternative=alternative)
    return p_value


def t_test(x, y):
    """
    Return p-value using t-test of x and y.

    :param x,y
    :param alternative: "two.sided", "less", "greater"
    """
    statistic, p_value = ttest_ind(x, y)
    return p_value
