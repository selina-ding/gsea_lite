"""
FDR-based p-value adjuster.
https://github.com/dceoy/fdra
"""

import numpy as np
import pandas as pd
from pathlib import Path


def bh_fdr_correction(pvalue):
    """
    Benjamini-Hochberg FDR correction

    Parameters
    ----------
    pvalue : array-like
        DEG结果表中的 pvalue 列

    Returns
    -------
    numpy.ndarray
        q-value (FDR)
    """

    pvalue = np.asarray(pvalue, dtype=float)

    if np.any((pvalue < 0) | (pvalue > 1)):
        raise ValueError("Invalid p-values")

    # 总检验数
    m = len(pvalue)

    # 按 pvalue 从小到大排序
    # np.argsort():返回排好序的元素在原数组中的索引（位置编号）
    sorted_index = np.argsort(pvalue) #返回排好序的元素在原数组pvalue中的位置编号
    sorted_pvalue = pvalue[sorted_index]

    # BH校正
    qvalue_sorted = sorted_pvalue * m / np.arange(1, m + 1)

    # 保证FDR单调递增
    # 一个基因的qvalue不能大于排在它后面的任何一个基因的qvalue。如果大于了，就得降到和后面一样低
    # np.minimum.accumulate()：沿给定轴计算元素的累积最小值，返回一个数组，其中每个位置包含到目前为止遇到的最小值。
    qvalue_sorted = np.minimum.accumulate(
        qvalue_sorted[::-1] #先反转，pvalue变成从大到小排
    )[::-1] #最后再反转回从小到大

    # 限制在[0,1]
    # np.clip(): 强制将所有算出来的qvalue限制在0到1的范围内。如果数组里有大于1的数，强制变成1。
    qvalue_sorted = np.clip(qvalue_sorted, 0, 1)

    # 恢复原顺序
    qvalue = np.empty(m)
    qvalue[sorted_index] = qvalue_sorted #qvalue_sorted是从小到大排好序的，sorted_index是原pvalue的索引

    return qvalue
