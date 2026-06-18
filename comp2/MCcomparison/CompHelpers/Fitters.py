from iminuit import Minuit
import numpy as np

# ## example usage
# res = fit_quadratic_minimum(lyList, redChiOnly)
# print(res['success'])

# quadratic: y = a*(x - x0)^2 + y0

def quadratic(x, a, x0, y0):
    return a * (x - x0)**2 + y0

def fit_quadratic_minimum(xdata, ydata, option=None):
    """Fit a quadratic y = a*(x - x0)^2 + y0 near the minimum of a curve.
    
    Parameters
    ----------
    xdata, ydata : array-like
        Full x and y arrays (e.g. lyList, redChiOnly).
    option : any
        Unused option flag (kept for API compatibility).
    
    Returns
    -------
    result : dict with keys
        'x0'     : fitted minimum x  ± error
        'y0'     : fitted minimum y  ± error
        'a'      : curvature         ± error
        'x0_err', 'y0_err', 'a_err' : 1-sigma errors
        'minuit'  : the Minuit object (for further inspection)
        'success' : bool, whether the fit converged
        'mask'    : boolean mask of points used
    """
    xdata = np.asarray(xdata, dtype=float)
    ydata = np.asarray(ydata, dtype=float)
    if xdata.shape != ydata.shape:
        raise ValueError("xdata and ydata must have the same shape.")
    if xdata.size < 3:
        raise ValueError("Need at least 3 points for a quadratic fit.")
    
    order = np.argsort(xdata)
    x_sorted = xdata[order]
    y_sorted = ydata[order]

    imin = int(np.argmin(y_sorted))
    if imin == 0 or imin == len(y_sorted) - 1:
        raise ValueError("Minimum is at edge; need one point on each side for fit.")

    idx_sorted = np.array([imin - 1, imin, imin + 1], dtype=int)
    xfit = x_sorted[idx_sorted]
    yfit = y_sorted[idx_sorted]

    mask_sorted = np.zeros_like(x_sorted, dtype=bool)
    mask_sorted[idx_sorted] = True
    mask = np.zeros_like(mask_sorted)
    mask[order] = mask_sorted
    x_min_raw = float(x_sorted[imin])
    
    def chi2(*params):
        a, x0, y0 = params
        return np.sum((yfit - quadratic(xfit, a, x0, y0))**2)

    m = Minuit(chi2, 10.0, x_min_raw, y_sorted[imin], name=("a", "x0", "y0"))
    m.limits['x0'] = (xfit.min(), xfit.max())
    m.limits['y0'] = (yfit.min(), yfit.max())
    
    m.errordef = 1.0  # least-squares
    m.migrad()
    m.hesse()
    
    result = {
        'x0':      m.values['x0'],
        'y0':      m.values['y0'],
        'a':       m.values['a'],
        'x0_err':  m.errors['x0'],
        'y0_err':  m.errors['y0'],
        'a_err':   m.errors['a'],
        'success': m.valid,
        'minuit':  m,
        'mask':    mask,
    }

    outstr = "Fit results: x0 = {x0:.3f} ± {x0_err:.3f}, y0 = {y0:.3f} ± {y0_err:.3f}, a = {a:.3f} ± {a_err:.3f}, success = {success}".format(**result)
    
    return result, outstr

def quadHalfWidth(a):
    """Calculate the width of a quadratic at y = y0 + 1 (i.e. where chi2 increases by 1)"""
    if a <= 0:
        raise ValueError("Curvature a must be positive to define a width.")
    return np.sqrt(1/a)