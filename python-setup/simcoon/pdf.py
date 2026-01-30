"""
Probability Distribution Function Module
========================================

Tools for working with Probability Distribution Functions (PDF).
Replaces the legacy PDF executable.

This module provides tools for:
- Common distribution functions (log-normal, Weibull, etc.)
- PDF discretization for homogenization
- Distribution fitting from data
- Visualization utilities

Typical applications:
- Grain size distributions
- Fiber diameter distributions
- Particle size distributions
- Porosity distributions

Examples
--------
>>> import numpy as np
>>> from simcoon.pdf import lognormal, discretize_pdf, fit_distribution
>>>
>>> # Create log-normal grain size distribution
>>> grain_sizes = np.linspace(1, 100, 1000)
>>> pdf_values = lognormal(grain_sizes, mu=3.0, sigma=0.5)
>>>
>>> # Discretize for homogenization
>>> bins = discretize_pdf(
...     lambda x: lognormal(x, 3.0, 0.5),
...     x_min=1, x_max=100, n_bins=10
... )
"""

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Callable, Tuple, Union
import numpy as np

try:
    from scipy import stats
    from scipy.integrate import quad
    from scipy.optimize import minimize_scalar, curve_fit
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# =============================================================================
# Distribution Functions
# =============================================================================

def lognormal(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    """
    Log-normal probability density function.

    Common for grain sizes in polycrystalline materials.

    Parameters
    ----------
    x : ndarray
        Values at which to evaluate PDF (must be > 0)
    mu : float
        Mean of the underlying normal distribution (log-scale)
    sigma : float
        Standard deviation of the underlying normal distribution

    Returns
    -------
    ndarray
        PDF values

    Notes
    -----
    The mean grain size is exp(mu + sigma^2/2)
    The mode is exp(mu - sigma^2)
    """
    if HAS_SCIPY:
        return stats.lognorm.pdf(x, s=sigma, scale=np.exp(mu))
    else:
        x = np.asarray(x)
        return np.where(
            x > 0,
            np.exp(-0.5 * ((np.log(x) - mu) / sigma)**2) / (x * sigma * np.sqrt(2 * np.pi)),
            0.0
        )


def weibull(x: np.ndarray, shape: float, scale: float) -> np.ndarray:
    """
    Weibull probability density function.

    Common for fiber strength and fatigue life distributions.

    Parameters
    ----------
    x : ndarray
        Values at which to evaluate PDF (must be >= 0)
    shape : float
        Shape parameter (k > 0)
    scale : float
        Scale parameter (lambda > 0)

    Returns
    -------
    ndarray
        PDF values
    """
    if HAS_SCIPY:
        return stats.weibull_min.pdf(x, c=shape, scale=scale)
    else:
        x = np.asarray(x)
        k, lam = shape, scale
        return np.where(
            x >= 0,
            (k / lam) * (x / lam)**(k - 1) * np.exp(-(x / lam)**k),
            0.0
        )


def normal(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    """
    Normal (Gaussian) probability density function.

    Parameters
    ----------
    x : ndarray
        Values at which to evaluate PDF
    mu : float
        Mean
    sigma : float
        Standard deviation

    Returns
    -------
    ndarray
        PDF values
    """
    if HAS_SCIPY:
        return stats.norm.pdf(x, loc=mu, scale=sigma)
    else:
        x = np.asarray(x)
        return np.exp(-0.5 * ((x - mu) / sigma)**2) / (sigma * np.sqrt(2 * np.pi))


def gamma_dist(x: np.ndarray, shape: float, scale: float) -> np.ndarray:
    """
    Gamma probability density function.

    Parameters
    ----------
    x : ndarray
        Values at which to evaluate PDF (must be >= 0)
    shape : float
        Shape parameter (k > 0)
    scale : float
        Scale parameter (theta > 0)

    Returns
    -------
    ndarray
        PDF values
    """
    if HAS_SCIPY:
        return stats.gamma.pdf(x, a=shape, scale=scale)
    else:
        from math import gamma as gamma_func
        x = np.asarray(x)
        k, theta = shape, scale
        return np.where(
            x >= 0,
            x**(k - 1) * np.exp(-x / theta) / (theta**k * gamma_func(k)),
            0.0
        )


def uniform(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """
    Uniform probability density function.

    Parameters
    ----------
    x : ndarray
        Values at which to evaluate PDF
    a : float
        Lower bound
    b : float
        Upper bound

    Returns
    -------
    ndarray
        PDF values
    """
    x = np.asarray(x)
    return np.where((x >= a) & (x <= b), 1.0 / (b - a), 0.0)


def bimodal(x: np.ndarray, mu1: float, sigma1: float, mu2: float, sigma2: float,
            weight1: float = 0.5) -> np.ndarray:
    """
    Bimodal (mixture of two normals) probability density function.

    Parameters
    ----------
    x : ndarray
        Values at which to evaluate PDF
    mu1, sigma1 : float
        Parameters for first mode
    mu2, sigma2 : float
        Parameters for second mode
    weight1 : float
        Weight of first mode (0 < weight1 < 1)

    Returns
    -------
    ndarray
        PDF values
    """
    return weight1 * normal(x, mu1, sigma1) + (1 - weight1) * normal(x, mu2, sigma2)


def gaussian_mixture(x: np.ndarray, means: np.ndarray,
                     stds: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """
    Gaussian mixture model PDF.

    Parameters
    ----------
    x : ndarray
        Values at which to evaluate PDF
    means : ndarray
        Mean of each component
    stds : ndarray
        Standard deviation of each component
    weights : ndarray
        Weight of each component (will be normalized)

    Returns
    -------
    ndarray
        PDF values
    """
    means = np.asarray(means)
    stds = np.asarray(stds)
    weights = np.asarray(weights)
    weights = weights / weights.sum()

    pdf = np.zeros_like(np.asarray(x), dtype=float)
    for mu, sigma, w in zip(means, stds, weights):
        pdf += w * normal(x, mu, sigma)
    return pdf


# =============================================================================
# Discretization
# =============================================================================

def discretize_pdf(pdf_func: Callable[[np.ndarray], np.ndarray],
                   x_min: float, x_max: float, n_bins: int,
                   total_volume: float = 1.0,
                   method: str = 'integrate') -> Dict[str, np.ndarray]:
    """
    Discretize a PDF into bins for homogenization.

    Parameters
    ----------
    pdf_func : callable
        PDF function f(x) -> density
    x_min : float
        Minimum value
    x_max : float
        Maximum value
    n_bins : int
        Number of discrete bins
    total_volume : float
        Total volume fraction (default 1.0)
    method : str
        'integrate' (accurate) or 'midpoint' (fast)

    Returns
    -------
    dict
        'centers': bin center values
        'edges': bin edge values
        'widths': bin widths
        'volume_fractions': volume fraction in each bin

    Examples
    --------
    >>> # Discretize log-normal grain size distribution
    >>> bins = discretize_pdf(
    ...     lambda x: lognormal(x, mu=3.0, sigma=0.5),
    ...     x_min=1, x_max=100, n_bins=10
    ... )
    >>> print(bins['centers'])
    >>> print(bins['volume_fractions'])
    """
    edges = np.linspace(x_min, x_max, n_bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    widths = edges[1:] - edges[:-1]

    if method == 'integrate' and HAS_SCIPY:
        # Integrate PDF over each bin
        fractions = np.array([
            quad(pdf_func, e1, e2)[0]
            for e1, e2 in zip(edges[:-1], edges[1:])
        ])
    else:
        # Midpoint approximation
        fractions = pdf_func(centers) * widths

    # Normalize
    fractions = fractions / fractions.sum() * total_volume

    return {
        'centers': centers,
        'edges': edges,
        'widths': widths,
        'volume_fractions': fractions
    }


def discretize_to_phases(pdf_func: Callable, x_min: float, x_max: float,
                         n_phases: int, property_name: str = 'size') -> List[Dict]:
    """
    Discretize PDF into phase definitions for RVE.

    Parameters
    ----------
    pdf_func : callable
        PDF function
    x_min, x_max : float
        Value range
    n_phases : int
        Number of phases
    property_name : str
        Name of the varying property

    Returns
    -------
    list of dict
        Phase definitions with 'volume_fraction' and property value
    """
    bins = discretize_pdf(pdf_func, x_min, x_max, n_phases)

    phases = []
    for i in range(n_phases):
        phases.append({
            'phase_number': i,
            property_name: bins['centers'][i],
            'volume_fraction': bins['volume_fractions'][i]
        })

    return phases


# =============================================================================
# Distribution Fitting
# =============================================================================

def fit_lognormal(data: np.ndarray) -> Dict[str, float]:
    """
    Fit log-normal distribution to data.

    Parameters
    ----------
    data : ndarray
        Sample data (must be positive)

    Returns
    -------
    dict
        'mu': log-scale mean
        'sigma': log-scale standard deviation
        'mean': actual mean
        'median': actual median
    """
    if not HAS_SCIPY:
        # Simple moment estimation
        log_data = np.log(data[data > 0])
        mu = np.mean(log_data)
        sigma = np.std(log_data)
    else:
        sigma, _, scale = stats.lognorm.fit(data, floc=0)
        mu = np.log(scale)

    return {
        'mu': mu,
        'sigma': sigma,
        'mean': np.exp(mu + sigma**2 / 2),
        'median': np.exp(mu)
    }


def fit_weibull(data: np.ndarray) -> Dict[str, float]:
    """
    Fit Weibull distribution to data.

    Parameters
    ----------
    data : ndarray
        Sample data (must be non-negative)

    Returns
    -------
    dict
        'shape': shape parameter (k)
        'scale': scale parameter (lambda)
        'mean': distribution mean
    """
    if not HAS_SCIPY:
        raise ImportError("scipy is required for Weibull fitting")

    shape, _, scale = stats.weibull_min.fit(data, floc=0)

    from math import gamma as gamma_func
    mean = scale * gamma_func(1 + 1/shape)

    return {
        'shape': shape,
        'scale': scale,
        'mean': mean
    }


def fit_normal(data: np.ndarray) -> Dict[str, float]:
    """
    Fit normal distribution to data.

    Parameters
    ----------
    data : ndarray
        Sample data

    Returns
    -------
    dict
        'mu': mean
        'sigma': standard deviation
    """
    mu = np.mean(data)
    sigma = np.std(data)
    return {'mu': mu, 'sigma': sigma}


def fit_distribution(data: np.ndarray, dist_type: str = 'auto') -> Dict[str, any]:
    """
    Fit a distribution to data.

    Parameters
    ----------
    data : ndarray
        Sample data
    dist_type : str
        Distribution type: 'normal', 'lognormal', 'weibull', 'gamma', or 'auto'
        If 'auto', tries multiple distributions and returns best fit.

    Returns
    -------
    dict
        'type': distribution type
        'params': fitted parameters
        'aic': Akaike Information Criterion (lower is better)
    """
    if not HAS_SCIPY:
        raise ImportError("scipy is required for distribution fitting")

    if dist_type == 'auto':
        # Try multiple distributions
        results = []
        for dt in ['normal', 'lognormal', 'weibull', 'gamma']:
            try:
                result = fit_distribution(data, dt)
                results.append(result)
            except:
                pass

        # Return best fit (lowest AIC)
        return min(results, key=lambda x: x['aic'])

    data = np.asarray(data)
    n = len(data)

    if dist_type == 'normal':
        params = fit_normal(data)
        log_lik = np.sum(stats.norm.logpdf(data, params['mu'], params['sigma']))
        k = 2

    elif dist_type == 'lognormal':
        params = fit_lognormal(data)
        log_lik = np.sum(stats.lognorm.logpdf(data, params['sigma'], scale=np.exp(params['mu'])))
        k = 2

    elif dist_type == 'weibull':
        params = fit_weibull(data)
        log_lik = np.sum(stats.weibull_min.logpdf(data, params['shape'], scale=params['scale']))
        k = 2

    elif dist_type == 'gamma':
        shape, _, scale = stats.gamma.fit(data, floc=0)
        params = {'shape': shape, 'scale': scale}
        log_lik = np.sum(stats.gamma.logpdf(data, shape, scale=scale))
        k = 2

    else:
        raise ValueError(f"Unknown distribution type: {dist_type}")

    # AIC = 2k - 2*log(L)
    aic = 2 * k - 2 * log_lik

    return {
        'type': dist_type,
        'params': params,
        'aic': aic,
        'log_likelihood': log_lik
    }


# =============================================================================
# Statistics
# =============================================================================

def pdf_statistics(pdf_func: Callable, x_min: float, x_max: float) -> Dict[str, float]:
    """
    Compute statistics of a PDF.

    Parameters
    ----------
    pdf_func : callable
        PDF function
    x_min, x_max : float
        Integration range

    Returns
    -------
    dict
        'mean': expected value
        'variance': variance
        'std': standard deviation
        'mode': mode (approximate)
    """
    if not HAS_SCIPY:
        raise ImportError("scipy is required for PDF statistics")

    # Mean
    mean, _ = quad(lambda x: x * pdf_func(x), x_min, x_max)

    # Variance
    moment2, _ = quad(lambda x: x**2 * pdf_func(x), x_min, x_max)
    variance = moment2 - mean**2

    # Mode (find maximum)
    result = minimize_scalar(lambda x: -pdf_func(x), bounds=(x_min, x_max), method='bounded')
    mode = result.x

    return {
        'mean': mean,
        'variance': variance,
        'std': np.sqrt(variance),
        'mode': mode
    }


# =============================================================================
# Visualization
# =============================================================================

def plot_pdf(x: np.ndarray, pdf_values: np.ndarray, ax=None,
             histogram_data: np.ndarray = None, n_bins: int = 30, **kwargs):
    """
    Plot a PDF with optional histogram comparison.

    Parameters
    ----------
    x : ndarray
        x values
    pdf_values : ndarray
        PDF values
    ax : matplotlib.axes.Axes, optional
        Axes to plot on
    histogram_data : ndarray, optional
        Data to overlay as histogram
    n_bins : int
        Number of histogram bins
    **kwargs
        Additional arguments passed to plot

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(x, pdf_values, **kwargs)

    if histogram_data is not None:
        ax.hist(histogram_data, bins=n_bins, density=True, alpha=0.3, label='Data')

    ax.set_xlabel('Value')
    ax.set_ylabel('Probability Density')

    return ax


def plot_discretization(bins: Dict, ax=None, **kwargs):
    """
    Plot discretized PDF as bar chart.

    Parameters
    ----------
    bins : dict
        Output from discretize_pdf()
    ax : matplotlib.axes.Axes, optional
        Axes to plot on
    **kwargs
        Additional arguments passed to bar

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()

    ax.bar(bins['centers'], bins['volume_fractions'],
           width=bins['widths'] * 0.9, **kwargs)

    ax.set_xlabel('Value')
    ax.set_ylabel('Volume Fraction')

    return ax


# =============================================================================
# Convenience Exports
# =============================================================================

__all__ = [
    # Distribution functions
    'lognormal',
    'weibull',
    'normal',
    'gamma_dist',
    'uniform',
    'bimodal',
    'gaussian_mixture',
    # Discretization
    'discretize_pdf',
    'discretize_to_phases',
    # Fitting
    'fit_lognormal',
    'fit_weibull',
    'fit_normal',
    'fit_distribution',
    # Statistics
    'pdf_statistics',
    # Visualization
    'plot_pdf',
    'plot_discretization',
]
