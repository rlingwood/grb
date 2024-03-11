import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

wavelengths_ang = np.array([21705.5, 16468.4, 12399.2, 8989.6, 7640.7, 6219.8, 5468, 4586.9, 4392, 3465, 2600, 2246, 1928])
wavelength = wavelengths_ang*10**(-10) # wavelengths in meters

class parameters: # define a set of parameters for the extinction curve
    def __init__(self, a, b, lam, n):
        self.a = a
        self.b = b
        self.lam = lam
        self.n = n
        
class galaxy:
    def __init__(self, bkg, fuv, fi, fii, fiii, fir):
        self.bkg = bkg
        self.fuv = fuv
        self.fi = fi
        self.fii = fii
        self.fiii = fiii
        self.fir = fir
        
# milky way parameters
bkg = parameters(165.0, 90.0, 0.047, 2.0) # parameters for BKG
fuv = parameters(14.0, 4.0, 0.08, 6.5) # parameters for FUV
fi = parameters(0.045, -1.95, 0.22, 2.0) # parameters for 2175 angstrom feature
fii = parameters(0.002, -1.95, 9.7, 2.0) # parameters for 9.7 micrometer feature
fiii = parameters(0.002, -1.8, 18.0, 2.0) # parameters for 18 micrometer feature
fir = parameters(0.012, 0.0, 25.0, 2.0) # parameters for FIR
mw = galaxy(bkg, fuv, fi, fii, fiii, fir)

def term(wavelen, parameters): # a single term in the extinction curve equation
    microns = wavelen*10**(6) # convert meter values to microns
    a = parameters.a
    b = parameters.b
    lam = parameters.lam
    n = parameters.n
    return a/((microns/lam)**n + (lam/microns)**n + b)

# add all terms of the extinction curve equation
def eta(wavelen, galaxy):
    return term(wavelen, galaxy.bkg) + term(wavelen, galaxy.fuv) + term(wavelen, galaxy.fi) + term(wavelen, galaxy.fii)  + term(wavelen, galaxy.fiii) + term(wavelen, galaxy.fir)

v = 5500*10**(-10)
rv = 3.08
err = 0.1
beta = 1

def flux(wavelength, ext, beta):
    av = rv*ext
    alambda = (eta(wavelength, mw)/eta(v, mw))*av
    tau = alambda/1.086
    initialflux = wavelength**beta
    flux = initialflux*(np.e**(-1*tau))
    noise = np.random.normal(0, err*flux)
    fluxwnoise = flux + noise
    return fluxwnoise

def model(wavelength, ebv, beta):
    modelav = rv*ebv
    modelambda = (eta(wavelength, mw)/eta(v, mw))*modelav
    modeltau = modelambda/1.086
    initialflux = wavelength**beta
    modelflux = initialflux*(np.e**(-1*modeltau))
    return modelflux

def fit(param):
    fittedav = rv*param[0]
    fittedalambda = (eta(wavelength, mw)/eta(v, mw))*fittedav
    fittedtau = fittedalambda/1.086
    initialflux = wavelength**param[1]
    fittedflux = initialflux*(np.e**(-1*fittedtau))
    return fittedflux

fig, ax = plt.subplots(figsize=(8, 6))
norm = plt.Normalize(np.log10(1/(7500*10**(-10))), np.log10(1/(3800*10**(-10))))
ax.set(title = "Flux for Milky Way",
       xlabel = "log[1/λ] (m)",
       ylabel = "log[Flux]")
plt.rcParams['font.family'] = "serif"

# reduced chi squared
def chi(model, data, unc, var):
    sum = 0
    n = 0
    while n < len(model):
        sum = sum + ((model[n]-data[n])/unc[n])**2
        n = n + 1
    dof = len(a) - var # flux.__code__.co_argcount
    return sum/dof

x = np.log10(1/wavelength)
inputebv = np.zeros(11)
outputebv = np.zeros(11)
outputbeta = np.zeros(11)
for i in range(0, 11):
    j = i/10
    np.put(inputebv, [i], [j])
    y = np.log10(flux(wavelength, j, beta)) # plots given ebv and beta, not fitted
    param, param_cov = curve_fit(model, wavelength, flux(wavelength, j, beta))
    np.put(outputebv, [i], [param[0]])
    np.put(outputbeta, [i], [param[1]])
    plt.scatter(x, y, s=10, c=x, cmap='nipy_spectral_r', norm=norm)
    plt.plot(x, np.log10(fit(param)), c = 'k', lw = 0.1)
    model = fit(param)
    data = flux(wavelength, j, beta)
    unc = np.random.normal(0, err*flux)
    var = 3
    fitqual = chi(model, data, unc, var)
plt.show()

# plot fitted ebv against fitted beta
fig2, ax = plt.subplots(figsize=(8, 6))
plt.scatter(outputebv, outputbeta, s=10, c='k')
a, b = np.polyfit(outputebv, outputbeta, 1)
plt.plot(outputebv, a*outputebv+b, c = 'k', lw = 0.2)  




plt.show()
# this seems to only show an actual degeneracy if the inputted beta is very large