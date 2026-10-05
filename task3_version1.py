import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from astropy.modeling import models, fitting
from functions_ajuste import *
from functions_own import *
import sys 


obs = pd.read_csv(
    'spectrum.txt',
    skiprows=27, # It skips first 27 lines 
    sep=','
)

obs.columns = obs.columns.str.strip()

wave_ang = obs.iloc[:, 0]
flux_adu = obs.iloc[:, 1]

red = 0.0 #reshift is 0
wave = wave_ang / (1 + red)
flux = flux_adu

plt.plot(wave, flux, label='Spectrum')
plt.xlabel('Wavelenght  [$\AA$]')
plt.ylabel('Flux [ADU]')
plt.legend()

lam_um = wave       
flambda = flux

#Configuration for the continuum
#Continuum regions
wl1, wl2 = 6670, 6680
wl3, wl4 = 6698, 6705

flux_nocont, lam_cont, flux_cont, poly_cont = ajuste_continuo (wl1, wl2, wl3, wl4, lam_um, flambda) #we use the function ajuste_continuo to get it 

resid_cont = flux_cont - poly_cont(lam_cont) #residuals
noise = np.std(resid_cont)

#GRAPHICS

plt.figure(figsize=(8,5))
plt.title('UVES Spectrum (Units: Å and ADU)')
plt.plot(lam_um, flambda, label="Original spectrum")
plt.plot(lam_um, poly_cont(lam_um), label = 'Continuum')
plt.plot(lam_um, flux_nocont, label="Spectrum without continuum")
plt.xlabel("Wavelength [$\AA$]")
plt.ylabel("Flux [ADU]")
plt.legend()
# We adjust the limits of the x axis to a real range of the spectrum in Angstroms
plt.xlim(6660, 6720)
plt.show()

#Selection of a specific region of the spectrum
w1 = 6670.0
w2 = 6720.0
mask = (wave_ang > w1) & (wave_ang < w2) 

flux_reg = flambda[mask]
wave_reg = wave_ang[mask]
cont_reg_model = poly_cont(wave_reg)

flux_reg_nocont = flux_reg - cont_reg_model

flux_reg_nocont=flux_reg - cont_reg_model

#It adjusts the rest value
wave_rest = 6685.66  # [Å]

#Estimates the initial width for the gaussian in velocity (~400 km/s)
stddev_n_line = (400.0 / 2.9979e5) * wave_rest / 2.35  

# Without normalisation because we are working in ADU
flux_norm = 1.0
flux_reg_nocont_norm = flux_reg_nocont / flux_norm

anchura_ins = 0.0  #instrumental resolution [km/s]  if you know it 

# GAUSSIAN MODEL DEFINITION (1COMPONENT) 
# Approximated initial amplitude (~45 ADU if we look at the detected peak)
line_n = models.Gaussian1D(
    amplitude=45.0, 
    mean=wave_rest, 
    stddev=stddev_n_line, 
    name='line-n'
)
line_n.amplitude.min = 0.0

combo_3c = line_n

# CALCULTATION AND ERRORS

combo_3c_fit_error = np.zeros([3])
combo_3c_fit, combo_3c_fit_error[:] = calc_bestfit_parameters(
    combo_3c, 
    wave_reg, 
    flux_reg_nocont_norm, 
    Noise=2 * noise / flux_norm, 
    Nsimul=20
)

combo_3c_mod = combo_3c_fit(wave_reg)

curspec = flux_reg_nocont_norm
curfit = combo_3c_fit
curfit_err = combo_3c_fit_error[:]
curmod = combo_3c_mod

# Rebuild the fitted model
line_n_fit = models.Gaussian1D(
    amplitude=curfit.amplitude, 
    mean=curfit.mean, 
    stddev=curfit.stddev
)
line_n_fit_err = assign_errors(curfit_err)

# PHYSICAL PARAMETERS CALCULus (FWHM, VELOCITY AND FLUX)
# FWHM in velocity(km/s)
FWHMr_line_n, FWHMr_line_n_err = calc_gaussian_fwhm(
    line_n_fit, velocity=True, Errpars=line_n_fit_err
)

fwhm_corr = np.sqrt(np.maximum(0, FWHMr_line_n**2 - anchura_ins**2))
print(f'FWHM_n (corrected) = {fwhm_corr:.2f} +/- {FWHMr_line_n_err:.2f} km/s')

# Velocity with respect to the wave_rest
vel_line_n, vel_line_n_err = calc_vel(
    line_n_fit, Errpars=line_n_fit_err, wave=wave_rest
)
print(f'Velocity_n = {vel_line_n:.2f} +/- {vel_line_n_err:.2f} km/s')

#Integrated flux(area under the gaussian in  ADU * Å)
flux_line_n, fluxr_line_n_err = calc_gaussian_flux(
    line_n_fit, flux_norm=flux_norm, Errpars=line_n_fit_err
)
print(f'Flux_n = {flux_line_n:.2f} +/- {fluxr_line_n_err:.2f} ADU*Å')

#GRAPHICS
plt.figure(figsize=(8, 5))

#Observed spectrum
plt.step(
    wave_reg, 
    flux_reg_nocont + cont_reg_model, 
    label='Obs', 
    color='royalblue', 
    alpha=0.8
)

# Gussian fitted component + Continuum
plt.plot(
    wave_reg, 
    line_n_fit(wave_reg) + cont_reg_model, 
    label='Gaussian Narrow component', 
    color='yellow', 
    ls='--', 
    alpha=0.9
)

#Total model
plt.plot(
    wave_reg, 
    curmod + cont_reg_model, 
    label='Model', 
    color='black', 
    alpha=0.8
)

#CONTINUUM
plt.plot(
    wave_reg, 
    poly_cont(wave_reg), 
    label = 'Continuum',
    color = 'red'
)

plt.ylabel("Flux [ADU]", size=12)
plt.xlabel(r"Wavelength [$\AA$]", size=12)
plt.legend()
plt.title('Spectral line fitting')
plt.tight_layout()
plt.show()

line_n_fit_err = assign_errors(combo_3c_fit_error)

# Extraction of the central parameters 
A_fit = line_n_fit.amplitude.value
mu_fit = line_n_fit.mean.value
sigma_fit = line_n_fit.stddev.value

A_err = line_n_fit_err['amplitude_err'][()]
mu_err = line_n_fit_err['mean_err'][()]
sigma_err = line_n_fit_err['stddev_err'][()]

# FWHM in Ångstroms 
fwhm_factor = np.sqrt(8 * np.log(2.0))  # ~2.35482
fwhm_ang = fwhm_factor * sigma_fit
fwhm_ang_err = fwhm_factor * sigma_err

# Calculations made with its own functions
FWHMr_line_n, FWHMr_line_n_err = calc_gaussian_fwhm(
    line_n_fit, velocity=True, Errpars=line_n_fit_err
)
fwhm_corr = np.sqrt(np.maximum(0, FWHMr_line_n**2 - anchura_ins**2))

vel_line_n, vel_line_n_err = calc_vel(
    line_n_fit, Errpars=line_n_fit_err, wave=wave_rest
)

flux_line_n, fluxr_line_n_err = calc_gaussian_flux(
    line_n_fit, flux_norm=flux_norm, Errpars=line_n_fit_err
)

#Parameters of the continuum
if hasattr(poly_cont, 'parameters'):
    params = poly_cont.parameters
    c0_fit = params[0] if len(params) > 0 else 0.0
    c1_fit = params[1] if len(params) > 1 else 0.0
elif hasattr(poly_cont, '__getitem__'):
    c0_fit = poly_cont[1] if len(poly_cont) > 1 else poly_cont[0]
    c1_fit = poly_cont[0] if len(poly_cont) > 1 else 0.0
else:
    c0_fit, c1_fit = 0.0, 0.0

c0_err = noise

print("\n" + "="*60)
print(" BEST-FIT PARAMETERS AND UNCERTAINTIES ")
print("="*60)
print("1. CONTINUUM FIT (Polynomial Baseline):")
print(f"   - Baseline Offset (c0):     {c0_fit:.3f} ± {c0_err:.3f} [ADU]")
print(f"   - Baseline Slope (c1):      {c1_fit:.5f} [ADU/Å]")
print(f"   - Continuum RMS Noise:      {noise:.3f} [ADU]")
print("-" * 60)
print("2. GAUSSIAN EMISSION LINE FIT:")
print(f"   - Amplitude (A):            {A_fit:.3f} ± {A_err:.3f} [ADU]")
print(f"   - Central Wavelength (μ):   {mu_fit:.3f} ± {mu_err:.3f} [Å]")
print(f"   - Standard Deviation (σ):   {sigma_fit:.3f} ± {sigma_err:.3f} [Å]")
print(f"   - FWHM (Wavelength):        {fwhm_ang:.3f} ± {fwhm_ang_err:.3f} [Å]")
print("-" * 60)
print("3. DERIVED PHYSICAL PARAMETERS:")
print(f"   - FWHM (Velocity):          {fwhm_corr:.2f} ± {FWHMr_line_n_err:.2f} [km/s]")
print(f"   - Centroid Velocity:        {vel_line_n:.2f} ± {vel_line_n_err:.2f} [km/s]")
print(f"   - Integrated Line Flux:     {flux_line_n:.2f} ± {fluxr_line_n_err:.2f} [ADU*Å]")
print("="*60 + "\n")