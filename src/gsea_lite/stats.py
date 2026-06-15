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
    pvalue 

    Returns
    -------
    numpy.ndarray
        q-value (FDR)
    """

    pvalue = np.asarray(pvalue, dtype=float)

    if np.any((pvalue < 0) | (pvalue > 1)):
        raise ValueError("Invalid p-values")

    # Total number of statistical tests
    m = len(pvalue)

    # sorting pvalue in ascending order
    sorted_index = np.argsort(pvalue) # Return the position number of the sorted element in the original array
    sorted_pvalue = pvalue[sorted_index]

    # Benjamini–Hochberg correction
    qvalue_sorted = sorted_pvalue * m / np.arange(1, m + 1)

    qvalue_sorted = np.minimum.accumulate(qvalue_sorted[::-1])[::-1]

    # Restrict all q-values to the interval [0, 1]
    qvalue_sorted = np.clip(qvalue_sorted, 0, 1)

    # Restore the q-values to the original p-value order
    qvalue = np.empty(m)
    qvalue[sorted_index] = qvalue_sorted

    return qvalue
