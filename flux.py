import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# wavelengths
wavelength = np.arange(1500*(10**(-10)), 10000*(10**(-10)), 200*(10**(-10))) # range of wavelengths (microns)
c = 3*(10**8) # speed of light (metres per second)
freq = c/wavelength # frequency of light (hertz)

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
    microns = wavelen*10**(6)
    a = parameters.a
    b = parameters.b
    lam = parameters.lam
    n = parameters.n
    return a/((microns/lam)**n + (lam/microns)**n + b)

# add all terms of the extinction curve equation
def eta(wavelen):
    return term(wavelen, bkg) + term(wavelen, fuv) + term(wavelen, fi) + term(wavelen, fii)  + term(wavelen, fiii) + term(wavelen, fir)

v = 5500*(10**(-10))
rv = 3.08

def flux(wavelength, ext):
    av = rv*ext
    alambda = (eta(wavelength)/eta(v))*av
    tau = alambda/1.086
    initialflux = wavelength
    flux = initialflux*(np.e**(-1*tau))
    fluxjansky = flux*((wavelength)**2)*3.34*(10**(4))
    return fluxjansky

def addnoise(wavelength, ext):
    noise = np.random.normal(0, np.sqrt(flux(wavelength, ext))*10**(-9))
    return flux(wavelength, ext) + noise

x = np.log10(1/wavelength)
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
    return initialflux*wavelength*(np.e**(-1*modeltau))*((wavelength)**2)*3.34*(10**(4))

# curve fitting
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

# plot
fig, ax = plt.subplots(figsize=(8, 6))
t = x
norm = plt.Normalize(np.log10(1/(7500*(10**(-10)))), np.log10(1/(3800*(10**(-10)))))
plt.scatter(x, yfive, s=10, c=t, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, yfour, s=10, c=t, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, ythree, s=10, c=t, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, ytwo, s=10, c=t, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, yone, s=10, c=t, cmap='nipy_spectral_r', norm=norm)
plt.scatter(x, yzero, s=10, c=t, cmap='nipy_spectral_r', norm=norm)

def fit(param):
    fittedextinction = param[0]
    fittedav = rv*fittedextinction
    fittedalambda = (eta(wavelength)/eta(v))*fittedav
    fittedtau = fittedalambda/1.086
    initialflux = param[1]*wavelength
    return initialflux*(np.e**(-1*fittedtau))*((wavelength)**2)*3.34*(10**(4))

plt.plot(x, np.log10(fit(paramfive)), c = 'k', lw = 0.4)
plt.plot(x, np.log10(fit(paramfour)), c = 'k', lw = 0.4)
plt.plot(x, np.log10(fit(paramthree)), c = 'k', lw = 0.4)
plt.plot(x, np.log10(fit(paramtwo)), c = 'k', lw = 0.4)
plt.plot(x, np.log10(fit(paramone)), c = 'k', lw = 0.4)
plt.plot(x, np.log10(fit(paramzero)), c = 'k', lw = 0.4)

# axis labels and formatting
ax.set(title = "Flux for Milky Way",
       xlabel = "log[1/λ] (m)",
       ylabel = "log[Flux] (Jy)")
plt.rcParams['font.family'] = "serif"

array = flux(wavelength, 0.0)
array2 = addnoise(wavelength,0.5)-flux(wavelength, 0.0)

# show plot
plt.show()
