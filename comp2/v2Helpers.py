import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares

def gaussian(x, amplitude, mean, sigma):
  """Gaussian function for fitting"""
  return amplitude * np.exp(-0.5 * ((x - mean) / sigma)**2)

def fit_histogram_gaussian(bin_centers, counts, errors=None, initial_guess=None):
  """
  Fit a histogram with a Gaussian using iminuit
  
  Parameters:
  bin_centers: array of bin center values
  counts: array of histogram counts
  errors: array of errors (if None, uses sqrt(counts))
  initial_guess: dict with keys 'amplitude', 'mean', 'sigma'
  
  Returns:
  m: fitted Minuit object
  """
  
  # Calculate errors if not provided
  if errors is None:
      errors = np.sqrt(counts)
      # Avoid zero errors
      errors[errors == 0] = 1
  
  # Initial parameter guesses if not provided
  if initial_guess is None:
    amplitude_guess = np.max(counts)
    mean_guess = bin_centers[np.argmax(counts)]
    sigma_guess = np.std(bin_centers[counts > np.max(counts)/2])
    
    initial_guess = {
        'amplitude': amplitude_guess,
        'mean': mean_guess, 
        'sigma': sigma_guess
    }
  
  # Create cost function
  cost = LeastSquares(bin_centers, counts, errors, gaussian)
  
  # Create Minuit object
  m = Minuit(cost, **initial_guess)
  
  # Set reasonable limits
  m.limits['amplitude'] = (0, None)
  m.limits['sigma'] = (0, None)
  
  # Perform the fit
  m.migrad()
  return m

def fit_mc_histogram(mc_obj, scale_factor=1.0, rangeLow = 0, rangeHigh = 100):
  """Simple wrapper to fit MC histogram with Gaussian"""
  # Get the MC histogram data
  mc_bin_centers = mc_obj.mBinX
  mc_counts = mc_obj.m_scaledCounts * scale_factor
  
  # Only fit in the range where we have reasonable statistics
  fit_mask = (mc_bin_centers >= rangeLow) & (mc_bin_centers <= rangeHigh) & (mc_counts > 0)
  mc_centers_fit = mc_bin_centers[fit_mask]
  mc_counts_fit = mc_counts[fit_mask]
  
  # Perform the fit
  fit_result = fit_histogram_gaussian(mc_centers_fit, mc_counts_fit)
  return fit_result