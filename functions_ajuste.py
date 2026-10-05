#done by Ana
import numpy as np
from numpy import random as rnd
from astropy.modeling import models, fitting

cspeed = 2.999e5


def indexes_gauss_comp(compoundmodel,icomp):
	param_names = ['amplitude_'+str(icomp),'mean_'+str(icomp),'stddev_'+str(icomp)]
	index_names =[]
	for name in param_names:
		index_names.append(compoundmodel.param_names.index(name))
	return index_names

def assign_errors(gauss_err):
#    err_ndarray= np.zeros((gauss_err[0],gauss_err[1],gauss_err[2]),\
#		dtype=[('amplitude_err','<f8'),('mean_err','<f8'),('stddev_err','<f8')])
	err_ndarray= np.array((gauss_err[0],gauss_err[1],gauss_err[2]),\
		dtype=[('amplitude_err','<f8'),('mean_err','<f8'),('stddev_err','<f8')])
	return err_ndarray

def calc_gaussian_fwhm(gaussian,velocity=True,Errpars=None):
	#print(gaussian)
	fwhm=np.sqrt(8*np.log(2.))*gaussian.stddev
	if velocity:
		fwhm *= (cspeed/gaussian.mean)
	if Errpars:
		fwhm_err = 	np.sqrt(8*np.log(2.))*Errpars['stddev_err']
		if velocity:
			fwhm_err *= (cspeed/gaussian.mean)
		return fwhm, fwhm_err
	return fwhm

def calc_vel(gaussian, Errpars=None, wave = None):
	#print(gaussian)
	vel = (gaussian.mean / wave -1) * cspeed #5486.961467393204
	if Errpars:
		vel_err = 	cspeed/wave*Errpars['mean_err']
		return vel, vel_err
	return vel

def calc_gaussian_flux(gaussian,flux_norm=1,Errpars=None):
	flux = np.sqrt(2*np.pi)*gaussian.stddev*flux_norm*gaussian.amplitude
	if Errpars:
		flux_err = np.sqrt(2*np.pi)*flux_norm*\
			np.sqrt(gaussian.stddev**2*Errpars['amplitude_err']**2+\
			gaussian.amplitude**2*Errpars['stddev_err']**2)
		return flux, flux_err   
	return flux

def calc_bestfit_parameters(fit_model,x,y,Noise=None,Nsimul=50):
	# compute best fit parameters and determine errors if noise is specified	

	fitter = fitting.TRFLSQFitter() #fitting.SLSQPLSQFitter()

	bestfit_pars = fitter(fit_model, x,y)
	bestfit_mod = bestfit_pars(x) 
	if (Noise != None):
		fit_model_init = bestfit_pars
		nparam = bestfit_pars.param_sets.size
		#print(' params_bestfit ',bestfit_pars.parameters)
		params_simul = np.zeros([nparam,Nsimul])
		# generate noisy spectra and repeat the fit
		
		for isimul in range(Nsimul):
			#print(' -------- ')
			#print(' executing simulation ',isimul)
			y_sim = rnd.normal(y,Noise)
			bestfit_pars_sim = fitter(fit_model_init,x,y_sim)
			params_simul[:,isimul]=bestfit_pars_sim.parameters
			#print(' params_simul ',bestfit_pars_sim.parameters)

		error_params=np.std(params_simul,axis=1) 
		return bestfit_pars, error_params
