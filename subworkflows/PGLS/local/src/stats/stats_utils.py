#!/usr/bin/env python3
# ============================ src/stats/stats_utils.py ===============================
# Author: Miguel Ramón (IBE-UPF)
# Date: 2025
# =====================================================================================
# This module provides statistical utility functions for data analysis in comparative studies.
# It includes functions to summarize numeric arrays, apply Benjamini-Hochberg FDR correction to p-values,
# perform FDR correction within groups in a DataFrame, and standardize arrays (z-scores).
# These utilities support robust statistical workflows for biological and evolutionary data.
#
# Main Features:
#   - summarize_series: Computes summary statistics (mean, median, sd, var, skewness, normality p-value, min, max) for a numeric array.
#   - bh_fdr: Applies Benjamini-Hochberg FDR correction to a vector of p-values.
#   - bh_fdr_per_group: Applies BH FDR correction within each group of a DataFrame, adds q-value and significance columns.
#   - standardize: Standardizes an array (z-scores), mean-centering if sd is zero or not finite.
#
# Dependencies: numpy, pandas, scipy.stats
# =====================================================================================

from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import skew, normaltest, kurtosis

def summarize_series(x: np.ndarray) -> dict:
    """
    Compute summary statistics for a numeric array.

    Parameters:
        x (np.ndarray): Input array of numbers.

    Returns:
        dict: Dictionary containing count, mean, median, standard deviation,
              variance, skewness, normality test p-value, min, and max.
              If no finite values, all statistics are NaN except n=0.
    """
    x = np.asarray(x, float)  # Ensure input is a float array
    finite = x[np.isfinite(x)]  # Filter out non-finite values (NaN, inf)
    if finite.size == 0:
        # Return all NaN if no finite values
        return dict(n=0, mean=np.nan, median=np.nan, sd=np.nan, var=np.nan,
                    skew=np.nan, normal_p=np.nan, min=np.nan, max=np.nan)
    # Compute statistics
    m = float(np.mean(finite))
    med = float(np.median(finite))
    sd = float(np.std(finite, ddof=1))  # Sample standard deviation
    var = float(np.var(finite, ddof=1)) # Sample variance
    sk = float(skew(finite)) if finite.size >= 3 else np.nan  # Skewness if enough data
    kurt = float(kurtosis(finite)) if finite.size >= 4 else np.nan  # Kurtosis if enough data
    try:
        # Normality test if enough data (>=8)
        k2, pnorm = normaltest(finite) if finite.size >= 8 else (np.nan, np.nan)
    except Exception:
        pnorm = np.nan
    # Return all statistics in a dictionary
    return dict(
        n=finite.size,
        mean=m,
        median=med,
        sd=sd,
        var=var,
        skew=sk,
        kurt=kurt,
        normal_p=float(pnorm),
        min=float(np.min(finite)),
        max=float(np.max(finite))
    )

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR correction for a vector of p-values.

    Parameters:
        pvals (np.ndarray): Array of p-values.

    Returns:
        np.ndarray: Array of FDR-adjusted q-values (same shape as input).
                    Non-finite p-values remain NaN.
    """
    p = np.asarray(pvals, float)  # Ensure input is float array
    q = np.full_like(p, np.nan, dtype=float)  # Output array initialized to NaN
    ok = np.isfinite(p)  # Mask of finite p-values
    m = int(ok.sum())    # Number of finite p-values
    if m == 0:
        return q  # Return all NaN if no finite p-values
    p_ok = p[ok]  # Extract finite p-values
    order = np.argsort(p_ok)  # Sort p-values
    ranks = np.arange(1, m + 1, dtype=float)  # Ranks for BH procedure
    q_ok = p_ok[order] * (m / ranks)  # Compute BH q-values
    q_ok = np.minimum.accumulate(q_ok[::-1])[::-1]  # Ensure monotonicity
    inv = np.empty_like(order)
    inv[order] = np.arange(m)  # Inverse permutation to restore original order
    q[ok] = q_ok[inv]  # Assign q-values to output array
    return np.minimum(q, 1.0)  # Cap q-values at 1.0

def bh_fdr_per_group(df: pd.DataFrame, pcol: str, group_col: str, alpha: float = 0.05) -> pd.DataFrame:
    """
    Apply Benjamini-Hochberg FDR correction within each group in a DataFrame.

    Parameters:
        df (pd.DataFrame): Input DataFrame.
        pcol (str): Name of the column containing p-values.
        group_col (str): Name of the column to group by.
        alpha (float): Significance threshold for q-values.

    Returns:
        pd.DataFrame: DataFrame with added columns:
            - 'q_{pcol}': FDR-adjusted q-values.
            - 'sig_{pcol}': Boolean indicating significance (q <= alpha).
    """
    df = df.copy()  # Work on a copy to avoid modifying original
    qcol = f'q_{pcol}'      # Name for q-value column
    sigcol = f'sig_{pcol}'  # Name for significance column
    df[qcol] = np.nan       # Initialize q-value column
    df[sigcol] = False      # Initialize significance column

    # Apply BH FDR correction within each group
    for _, group in df.groupby(group_col):
        idx = group.index
        pvals = df.loc[idx, pcol].values
        qvals = bh_fdr(pvals)
        df.loc[idx, qcol] = qvals
        df.loc[idx, sigcol] = qvals <= alpha

    return df

def standardize(x: np.ndarray) -> np.ndarray:
    """
    Standardize an array: subtract mean and divide by standard deviation.

    Parameters:
        x (np.ndarray): Input array.

    Returns:
        np.ndarray: Standardized array (z-scores).
                    If standard deviation is zero or not finite, only mean is subtracted.
    """
    x = np.asarray(x, float)  # Ensure input is float array
    mu = np.nanmean(x)        # Mean, ignoring NaNs
    sd = np.nanstd(x, ddof=1) # Sample standard deviation, ignoring NaNs
    # Standardize if sd is finite and > 0, else just mean-center
    return (x - mu) / sd if np.isfinite(sd) and sd > 0 else (x - mu)

def compute_weights(w_raw: np.ndarray, func: str, y: np.ndarray | None = None) -> np.ndarray:
    """
    Compute weights for GLS regression based on raw weights and function type.

    Parameters:
        w_raw (np.ndarray): Raw weight values (must be positive).
        func (str): Weight function ('linear', 'log', 'binomial').
        y (np.ndarray, optional): Response values for binomial weights.

    Returns:
        np.ndarray: Computed weights.

    Raises:
        ValueError: If weights are not positive or if y is required but not provided.
    """
    w = np.asarray(w_raw, float)
    if not np.all(w > 0):
        raise ValueError("Weights must be positive")
    if func == 'linear':
        return w
    if func == 'log':
        return np.log(w + 1e-6)
    if func == 'binomial':
        if y is None:
            raise ValueError("y required for binomial weights")
        y = np.asarray(y, float)
        return w * y * (1 - y)
    raise ValueError(f"Unknown weight function: {func}")
