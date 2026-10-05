#done by Ana
import argparse
import os
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def parse_arguments():
    parser = argparse.ArgumentParser(description="Gaussian and continuum fitting of the spectrum.")
    parser.add_argument('spectrum_file', type=str, help="Path of the file text of the spectrum (ej. spectrum.txt)")
    return parser.parse_args()

import os
import pandas as pd

def read_spectrum(filename):
    wave_unit, flux_unit = "Å", "ADU"
    
    clean_filename = filename.strip("'\"")
    
    if os.path.isabs(clean_filename):
        file_path = clean_filename
    else:
        file_path = os.path.abspath(clean_filename)
    
    # Verify if the file exists 
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"We couldn't find the file in this path: {file_path}")

    #Find wavelength column
    #with open(file_path, 'r') as f:
        #lineas = f.readlines()
        #start_line = next(i for i, line in enumerate(lineas) if 'WAVELENGTH' in line)
    #read it 
    obs = pd.read_csv(
        file_path, 
        skiprows=27,
        sep=r'\s+|,', 
        engine='python',
        on_bad_lines='skip'
    )
    
    obs.columns = obs.columns.str.strip()

    # Make columns into numpy arrays 
    wave_ang = obs.iloc[:, 0].to_numpy(dtype=float)
    flux_adu = obs.iloc[:, 1].to_numpy(dtype=float)

    red = 0.0  # Redshift
    wave = wave_ang / (1 + red)
    flux = flux_adu

    return wave, flux, wave_unit, flux_unit

#MODELING WITH scipy.optimize.curve_fit

def polynomial_baseline(x, c0, c1):
    """1st garde polynomial for the continuum."""
    return c0 + c1 * x

def gaussian_with_continuum(x, A, mu, sigma, c0, c1):
    """Complete model: Gaussian+ Continuum."""
    return (c0 + c1 * x) + A * np.exp(- (x - mu)**2 / (2 * sigma**2))

def fit_spectrum(x, y):
    # define the peak region and continuum regions
    line_mask = (x >= 6680.0) & (x <= 6695.0)
    cont_mask = ~line_mask
    
    # Continuum configuration
    popt_cont, pcov_cont = curve_fit(polynomial_baseline, x[cont_mask], y[cont_mask]) #parameter optimized y covarianza
    c0_init, c1_init = popt_cont
    perr_cont = np.sqrt(np.diag(pcov_cont)) #parameter error
    
    # initial estimations for errors parameters 
    x_line = x[line_mask]
    y_line = y[line_mask]
    
    max_idx = np.argmax(y_line)
    mu_init = x_line[max_idx]
    A_init = y_line[max_idx] - polynomial_baseline(mu_init, c0_init, c1_init)
    sigma_init = 2.0
    
    p0 = [A_init, mu_init, sigma_init, c0_init, c1_init]
    
    # non-lineal fitting with curve_fit
    popt_full, pcov_full = curve_fit(gaussian_with_continuum, x, y, p0=p0, maxfev=10000)
    perr_full = np.sqrt(np.diag(pcov_full))
    
    return popt_full, perr_full, popt_cont, perr_cont, cont_mask

def main():
    args = parse_arguments()
    x, y, wave_unit, flux_unit = read_spectrum(args.spectrum_file)
    
    popt, perr, popt_cont, perr_cont, cont_mask = fit_spectrum(x, y)
    
    A, mu, sigma, c0, c1 = popt
    A_err, mu_err, sigma_err, c0_err, c1_err = perr
    
    # FWHM and its uncertainties 
    fwhm_factor = 2.35482004503
    fwhm = fwhm_factor * sigma
    fwhm_err = fwhm_factor * sigma_err
    
    # Evalution of the models 
    poly_fit = polynomial_baseline(x, c0, c1)
    full_fit = gaussian_with_continuum(x, *popt)
    
    print("\n" + "="*55)
    print(" BEST-FIT PARAMETERS AND UNCERTAINTIES ")
    print("="*55)
    print("1. CONTINUUM POLYNOMIAL:")
    print(f"   - Baseline Offset (c0):     {c0:.3f} ± {c0_err:.3f} [{flux_unit}]")
    print(f"   - Baseline Slope (c1):      {c1:.5f} ± {c1_err:.5f} [{flux_unit}/{wave_unit}]")
    print("-" * 55)
    print("2. GAUSSIAN EMISSION LINE:")
    print(f"   - Amplitude (A):            {A:.3f} ± {A_err:.3f} [{flux_unit}]")
    print(f"   - Central Wavelength (μ):   {mu:.3f} ± {mu_err:.3f} [{wave_unit}]")
    print(f"   - Standard Deviation (σ):   {sigma:.3f} ± {sigma_err:.3f} [{wave_unit}]")
    print(f"   - FWHM:                     {fwhm:.3f} ± {fwhm_err:.3f} [{wave_unit}]")
    print("="*55 + "\n")
    
    #Graphics
    plt.figure(figsize=(8, 4))
    plt.plot(x, y, label="Full Spectrum", color="black", lw=1)
    plt.xlabel(f"Wavelength [{wave_unit}]")
    plt.ylabel(f"Flux [{flux_unit}]")
    plt.title("1. Full Spectrum")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 4))
    plt.plot(x, y, label="Spectrum", color="gray", alpha=0.7)
    plt.plot(x, poly_fit, color="red", lw=2, linestyle="--", label="Polynomial Continuum Fit")
    plt.xlabel(f"Wavelength [{wave_unit}]")
    plt.ylabel(f"Flux [{flux_unit}]")
    plt.title("2. Full Spectrum with Polynomial Baseline Overlaid")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 4))
    plt.plot(x, y, label="Spectrum", color="gray", alpha=0.7)
    plt.plot(x, poly_fit, color="red", lw=1.5, linestyle="--", label="Continuum")
    plt.plot(x, full_fit, color="blue", lw=2, label="Gaussian + Polynomial Model")
    plt.xlabel(f"Wavelength [{wave_unit}]")
    plt.ylabel(f"Flux [{flux_unit}]")
    plt.title("3. Full Spectrum with Polynomial and Gaussian Overlaid")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()