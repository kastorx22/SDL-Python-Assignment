#done by Ana
import numpy as np
from numpy import random as rnd
from astropy.modeling import models, fitting

#It calculates the continuum
def ajuste_continuo (wl1, wl2, wl3, wl4, lam_um, flambda):
        mask_cont = (
            ((lam_um >= wl1) & (lam_um <= wl2)) |
            ((lam_um >= wl3) & (lam_um <= wl4)) )
        #We fitt them usisng the mask
        lam_cont = lam_um[mask_cont]
        flux_cont = flambda[mask_cont]
        #polinomials
        linfitter = fitting.LinearLSQFitter()
        poly_cont = linfitter(models.Polynomial1D(1), lam_cont, flux_cont)
        #we substract it
        flux_nocont = flambda - poly_cont(lam_um) #flux with no continuum
        return flux_nocont, lam_cont, flux_cont, poly_cont

#To calculate the gaussian 
import numpy as np
from lmfit import Parameters, minimize


def gaussian_model(wave, vel, FWHM, fluxL, lambda0):
    c = 299792.458  # km/s

    mu = lambda0 * (1 + vel / c)

    sigma_v = FWHM / 2.3548
    sigma_l = mu * sigma_v / c

    amp = fluxL / (sigma_l * np.sqrt(2 * np.pi))

    return amp * np.exp(-(wave - mu)**2 / (2 * sigma_l**2)) 


def fit_gaussian_minimize(wave, flux, flux_err,
                          vel_init, FWHM_init, flux_init,
                          lambda0):

    c = 299792.458  # km/s

    def residual(params, wave, data, err):
        vel   = params['vel']
        FWHM  = params['FWHM']
        fluxL = params['flux']
        

        # Modelo usando la función definida
        model = gaussian_model(wave, vel, FWHM, fluxL, lambda0)

        if err is None:
            return data - model
        return (data - model) / err

    params = Parameters()
    params.add('vel',  value=vel_init)
    params.add('FWHM', value=FWHM_init, min=10)
    params.add('flux', value=flux_init, min=0)
    

    #minimizing
    result = minimize(residual, params, args=(wave, flux, flux_err))

    return result



       
