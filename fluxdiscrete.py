import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# wavelengths
wavelength = np.array([4640, 6580, 8060, 9000, 12200, 16300, 21900, 5510, 4450, 3650, 2910, 2310, 2120])
class parameters: # define a set of parameters for the curve
    def __init__(self, a, b, lam, n):
        self.a = a
        self.b = b
        self.lam = lam
        self.n = n

# milky way parameters
bkg = parameters(165.0, 90.0, 0.047, 2.0) # parameters for BKG
fuv = parameters(14.0, 4.0, 0.08, 6.5) # parameters for FUV
fi = parameters(0.045, -1.95, 0.22, 2.0) # parameters for 2175 angstrom feature
fii = parameters(0.002, -1.95, 9.7, 2.0) # parameters for 9.7 micrometer feature
fiii = parameters(0.002, -1.8, 18.0, 2.0) # parameters for 18 micrometer feature
fir = parameters(0.012, 0.0, 25.0, 2.0) # parameters for FIR

def term(wavelen, parameters): # a single term in the extinction curve equation
    microns = wavelen*10**(-4) # convert angstrom values to microns
    a = parameters.a
    b = parameters.b
    lam = parameters.lam
    n = parameters.n
    return a/((microns/lam)**n + (lam/microns)**n + b)

# add all terms of the extinction curve equation
def eta(wavelen):
    return term(wavelen, bkg) + term(wavelen, fuv) + term(wavelen, fi) + term(wavelen, fii)  + term(wavelen, fiii) + term(wavelen, fir)

v = 5500
rv = 3.08

def flux(wavelength, ext):
    av = rv*ext
    alambda = (eta(wavelength)/eta(v))*av
    tau = alambda/1.086
    initialflux = wavelength
    flux = initialflux*(np.e**(-1*tau))
    return flux

def addnoise(wavelength, ext):
    flx = flux(wavelength, ext)
    noise = np.random.normal(0, np.sqrt(flx))
    fluxwnoise = flx + noise
    fluxwnoisejy = fluxwnoise*((wavelength)**2)*3.34*(10**(4))
    return fluxwnoisejy

x = np.log10(1/wavelength)
ynine = np.log10(addnoise(wavelength, 0.9))
yeight = np.log10(addnoise(wavelength, 0.8))
yseven = np.log10(addnoise(wavelength, 0.7))
ysix = np.log10(addnoise(wavelength, 0.6))
yfive = np.log10(addnoise(wavelength, 0.5))
yfour = np.log10(addnoise(wavelength, 0.4))
ythree = np.log10(addnoise(wavelength, 0.3))
ytwo = np.log10(addnoise(wavelength, 0.2))
yone = np.log10(addnoise(wavelength, 0.1))
yzero = np.log10(addnoise(wavelength, 0.0))

def model(wavelength, ext, initialflux):
    modelav = rv*ext
    modelambda = (eta(wavelength)/eta(v))*modelav
    modeltau = modelambda/1.086
    modelflux = initialflux*wavelength*(np.e**(-1*modeltau))
    modelfluxjansky = modelflux*((wavelength)**2)*3.34*(10**(4))
    return modelfluxjansky

# curve fitting
paramnine, param_cov = curve_fit(model, wavelength, addnoise(wavelength, 0.9))
print("Fitted value for E(B-V): ", round(paramnine[0],4)) 
print("Fitted value for Initial Flux: ", round(paramnine[1],4)) 
print("") 
parameight, param_cov = curve_fit(model, wavelength, addnoise(wavelength, 0.8))
print("Fitted value for E(B-V): ", round(parameight[0],4)) 
print("Fitted value for Initial Flux: ", round(parameight[1],4)) 
print("") 
paramseven, param_cov = curve_fit(model, wavelength, addnoise(wavelength, 0.7))
print("Fitted value for E(B-V): ", round(paramseven[0],4)) 
print("Fitted value for Initial Flux: ", round(paramseven[1],4)) 
print("") 
paramsix, param_cov = curve_fit(model, wavelength, addnoise(wavelength, 0.6))
print("Fitted value for E(B-V): ", round(paramsix[0],4)) 
print("Fitted value for Initial Flux: ", round(paramsix[1],4)) 
print("") 
paramfive, param_cov = curve_fit(model, wavelength, addnoise(wavelength, 0.5))
print("Fitted value for E(B-V): ", round(paramfive[0],4)) 
print("Fitted value for Initial Flux: ", round(paramfive[1],4)) 
print("") 
paramfour, param_cov = curve_fit(model, wavelength, addnoise(wavelength, 0.4))
print("Fitted value for E(B-V): ", round(paramfour[0],4)) 
print("Fitted value for Initial Flux: ", round(paramfour[1],4)) 
print("")
paramthree, param_cov = curve_fit(model, wavelength, addnoise(wavelength, 0.3))
print("Fitted value for E(B-V): ", round(paramthree[0],4)) 
print("Fitted value for Initial Flux: ", round(paramthree[1],4))
print("") 
paramtwo, param_cov = curve_fit(model, wavelength, addnoise(wavelength, 0.2))
print("Fitted value for E(B-V): ", round(paramtwo[0],4)) 
print("Fitted value for Initial Flux: ", round(paramtwo[1],4)) 
print("")
paramone, param_cov = curve_fit(model, wavelength, addnoise(wavelength, 0.1))
print("Fitted value for E(B-V): ", round(paramone[0],4)) 
print("Fitted value for Initial Flux: ", round(paramone[1],4)) 
print("")
paramzero, param_cov = curve_fit(model, wavelength, addnoise(wavelength, 0.0))
print("Fitted value for E(B-V): ", round(paramzero[0],4)) 
print("Fitted value for Initial Flux: ", round(paramzero[1],4)) 

## plot
fig, ax = plt.subplots(figsize=(8, 6))
norm = plt.Normalize(np.log10(1/(7500)), np.log10(1/(3800)))
plt.scatter(x, ynine, s=10, c=x, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, yeight, s=10, c=x, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, yseven, s=10, c=x, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, ysix, s=10, c=x, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, yfive, s=10, c=x, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, yfour, s=10, c=x, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, ythree, s=10, c=x, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, ytwo, s=10, c=x, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, yone, s=10, c=x, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, yzero, s=10, c=x, cmap='nipy_spectral_r', norm=norm)

def fit(param):
    fittedextinction = param[0]
    fittedav = rv*fittedextinction
    fittedalambda = (eta(wavelength)/eta(v))*fittedav
    fittedtau = fittedalambda/1.086
    initialflux = param[1]*wavelength
    fittedflux = initialflux*(np.e**(-1*fittedtau))
    fittedfluxjansky = fittedflux*((wavelength)**2)*3.34*(10**(4))
    return fittedfluxjansky

plt.errorbar(x, np.log10(fit(paramnine)), yerr = np.log10(np.sqrt(ynine)), c = 'k', lw = 0.2)
plt.errorbar(x, np.log10(fit(parameight)), yerr = np.log10(np.sqrt(yeight)), c = 'k', lw = 0.2)
plt.errorbar(x, np.log10(fit(paramseven)), yerr = np.log10(np.sqrt(yseven)), c = 'k', lw = 0.2)
plt.errorbar(x, np.log10(fit(paramsix)), yerr = np.log10(np.sqrt(ysix)), c = 'k', lw = 0.2)
plt.errorbar(x, np.log10(fit(paramfive)), yerr = np.log10(np.sqrt(yfive)), c = 'k', lw = 0.2)
plt.errorbar(x, np.log10(fit(paramfour)), yerr = np.log10(np.sqrt(yfour)), c = 'k', lw = 0.2)
plt.errorbar(x, np.log10(fit(paramthree)), yerr = np.log10(np.sqrt(ythree)), c = 'k', lw = 0.2)
plt.errorbar(x, np.log10(fit(paramtwo)), yerr = np.log10(np.sqrt(ytwo)), c = 'k', lw = 0.2)
plt.errorbar(x, np.log10(fit(paramone)), yerr = np.log10(np.sqrt(yone)), c = 'k', lw = 0.2)
plt.errorbar(x, np.log10(fit(paramzero)), yerr = np.log10(np.sqrt(yzero)), c = 'k', lw = 0.2)

# axis labels and formatting
ax.set(title = "Flux for Milky Way",
       xlabel = "log[1/λ] (Å)",
       ylabel = "log[Flux] (Jy)")
plt.rcParams['font.family'] = "serif"

array = flux(wavelength, 0.0)
array2 = addnoise(wavelength,0.0)

# show plot
plt.show()
