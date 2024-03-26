import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# constants
c = 3*10**8
h = 6.63*10**(-34)
ev = 1.6*10**(-19)

beta = 1 # slope of intrinsic grb spectrum
nh = 0.01
rv = 3.08
v = 5500*10**(-10)
err = 0.1 # level of noise to add


xray_kev = np.arange(0.5, 10, 0.1) # xray energies in kev
optical_ang = np.array([21705.5, 16468.4, 12399.2, 8989.6, 7640.7, 6219.8, 5468, 4586.9, 4392, 3465, 2600, 2246, 1928]) # optical wavelengths in angstrom

def conv(e):
    # convert xray energies in kev to wavelengths in m
    m = e*1000*ev # convert kev to j
    return (h*c)/m

xray = conv(xray_kev) # get xray wavelengths in meters
optical = optical_ang*10**(-10) # get optical wavelengths in meters

# join xray + otpical arrays
wavelengths = np.concatenate((optical, xray), axis=0) # all wavelengths in meters



# !!! optical extinction !!!

class parameters: # define a set of parameters for the extinction curve
    def __init__(self, a, b, lam, n):
        self.a = a
        self.b = b
        self.lam = lam
        self.n = n
        
class galaxy: # define a set of sets of parameters for a galaxy
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

def flux(wavelength, ext, beta):
    av = rv*ext
    alambda = (eta(wavelength, mw)/eta(v, mw))*av
    tau = alambda/1.086
    initialflux = wavelength**beta
    flux = initialflux*(np.e**(-1*tau))
    noise = np.random.normal(0, err*flux)
    fluxwnoise = flux + noise
    return fluxwnoise



# !!! xray extinction !!!

xsections = fits.getdata('IsmabsAtomicData_reduced.fits',1)
ISMabund = {'H':12.0,'He':10.99,'C':8.38,'N':7.88,'O':8.69,'Fe':7.43,'Ne':7.94,'Mg':7.40,'Si':7.27,'S':7.09,'Ar':6.41,'Ca':6.20}
# beta = {'H':1.0,'He':1.0,'C':0.5,'N':1.0,'O':0.6,'Fe':0.3,'Ne':1.0,'Mg':0.2,'Si':0.1,'S':0.6,'Ar':1.0,'Ca':0.003}
H = xsections['H']*10.**(ISMabund['H']-12)#*beta['H']
He = xsections['HeI']*10.**(ISMabund['He']-12)#*beta['He']
C = xsections['CI']*10.**(ISMabund['C']-12)#*beta['C']
N = xsections['NI']*10.**(ISMabund['N']-12)#*beta['N']
O = xsections['OI']*10.**(ISMabund['O']-12)#%*beta['O']
Fe = xsections['FeI']*10.**(ISMabund['Fe']-12)#*beta['Fe']
Ne = xsections['NeI']*10.**(ISMabund['Ne']-12)#*beta['Ne']
Mg = xsections['MgI']*10.**(ISMabund['Mg']-12)#*beta['Mg']
Si = xsections['SiI']*10.**(ISMabund['Si']-12)#*beta['Si']
S = xsections['SI']*10.**(ISMabund['S']-12)#*beta['S']
Ar = xsections['ArI']*10.**(ISMabund['Ar']-12)#*beta['Ar']
Ca = xsections['CaI']*10.**(ISMabund['Ca']-12)#*beta['Ca']
keV = xsections['Energy']/1000.
totsigma = (H+He+C+O+Fe+Ne+Mg+Si+S+Ar+Ca)*1e24

dE = 0.1
sigma = np.zeros(95)

def attenuation(e):
     index = np.where((keV<e+dE) & (keV>e))
     s = totsigma[index[0]]
     return s[0]

for n, e in enumerate(xray_kev):
    sigma[n] = attenuation(e)
    
def xrayflux(xray, sigma): # this might not work?
    initialflux = xray**beta
    flux = initialflux*np.e**(-1*sigma*nh)
    noise = np.random.normal(0, err*flux)
    fluxwnoise = flux + noise
    return fluxwnoise



# !!! plot formatting !!!
fig, ax = plt.subplots(figsize=(10, 6))
norm = plt.Normalize(np.log10(1/(7500*10**(-10))), np.log10(1/(3800*10**(-10))))
ax.set(title = "Flux for Milky Way",
       xlabel = "log[1/λ] (m)",
       ylabel = "log[Flux]")
plt.rcParams['font.family'] = "serif"



# !!! curve fitting !!!

def model(wavelength, ebv, beta):
    modelav = rv*ebv
    modelambda = (eta(wavelength, mw)/eta(v, mw))*modelav
    modeltau = modelambda/1.086
    initialflux = wavelength**beta
    modelflux = initialflux*(np.e**(-1*modeltau))
    return modelflux

def fit(param):
    fittedav = rv*param[0]
    fittedalambda = (eta(optical, mw)/eta(v, mw))*fittedav
    fittedtau = fittedalambda/1.086
    initialflux = optical**param[1]
    fittedflux = initialflux*(np.e**(-1*fittedtau))
    return fittedflux



# !!! combining the two !!!
x = np.log10(1/wavelengths)

for i in range(0, 11): # range of ebvs
    j = i/10
    yopt = np.log10(flux(optical, j, beta)) # plots optical
    param, param_cov = curve_fit(model, optical, flux(optical, j, beta)) # optical curvefitting
    yxray = np.log10(xrayflux(xray, sigma)) # plots xray
    y = np.concatenate((yopt, yxray), axis=0) # combines the two
    plt.scatter(x, y, s=8, c=x, cmap='nipy_spectral_r', norm=norm)
