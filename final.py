import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit

# constants
c = 3*10**8 # speed of light
h = 6.63*10**(-34) # planck constant
ev = 1.6*10**(-19) # electron volt

beta = 1 # instrinsic slope of flux
zhost = 2 # redshift of host galaxy
zint = 0.5 # redshift of intervener

xray = np.arange(0.3, 10, 0.1) # xray energies in kev
optical_ang = np.array([21705.5, 16468.4, 12399.2, 8989.6, 7640.7, 6219.8, 5468, 4586.9, 4392, 3465]) # optical wavelengths in angstrom
optical = optical_ang*10**(-10) # optical wavelengths in m

def remove(energy, z):
    shiftedenergy = energy*(1+z)
    energy = energy[shiftedenergy < 10]
    return energy

xray = remove(xray, zhost) # remove xray energies that will be out of bounds once shifted

def conv(e):
    # convert xray energies in kev to wavelengths in m
    m = e*1000*ev # convert kev to j
    return (h*c)/m

class parameters: # define a set of parameters for optical extinction
    def __init__(self, a, b, lam, n):
        self.a = a
        self.b = b
        self.lam = lam
        self.n = n
        
class galaxy: # define a set of sets of parameters to describe a galaxy's full optical extinction curve
    def __init__(self, bkg, fuv, fi, fii, fiii, fir, rv, ratio):
        self.bkg = bkg
        self.fuv = fuv
        self.fi = fi
        self.fii = fii
        self.fiii = fiii
        self.fir = fir
        self.rv = rv
        self.ratio = ratio
        
# milky way parameters
bkg = parameters(165.0, 90.0, 0.047, 2.0) # parameters for BKG
fuv = parameters(14.0, 4.0, 0.08, 6.5) # parameters for FUV
fi = parameters(0.045, -1.95, 0.22, 2.0) # parameters for 2175 angstrom feature
fii = parameters(0.002, -1.95, 9.7, 2.0) # parameters for 9.7 micrometer feature
fiii = parameters(0.002, -1.8, 18.0, 2.0) # parameters for 18 micrometer feature
fir = parameters(0.012, 0.0, 25.0, 2.0) # parameters for FIR
mwrv = 3.08
mwratio = 0.18*10**(22)
mw = galaxy(bkg, fuv, fi, fii, fiii, fir, mwrv, mwratio)

# lmc parameters
lmcbkg = parameters(175.0, 90.0, 0.046, 2.0) # parameters for BKG
lmcfuv = parameters(19.0, 5.5, 0.08, 4.5) # parameters for FUV
lmcfi = parameters(0.023, -1.95, 0.22, 2.0) # parameters for 2175 angstrom feature
lmcfii = parameters(0.005, -1.95, 9.7, 2.0) # parameters for 9.7 micrometer feature
lmcfiii = parameters(0.006, -1.8, 18.0, 2.0) # parameters for 18 micrometer feature
lmcfir = parameters(0.02, 0.0, 25.0, 2.0) # parameters for FIR
lmcrv = 3.16
lmcratio = 0.7*10**(22)
lmc = galaxy(lmcbkg, lmcfuv, lmcfi, lmcfii, lmcfiii, lmcfir, lmcrv, lmcratio)

# smc parameters
smcbkg = parameters(185.0, 90.0, 0.042, 2.0) # parameters for BKG
smcfuv = parameters(27.0, 5.5, 0.08, 4.0) # parameters for FUV
smcfi = parameters(0.005, -1.95, 0.22, 2.0) # parameters for 2175 angstrom feature
smcfii = parameters(0.01, -1.95, 9.7, 2.0) # parameters for 9.7 micrometer feature
smcfiii = parameters(0.012, -1.8, 18.0, 2.0) # parameters for 18 micrometer feature
smcfir = parameters(0.03, 0.0, 25.0, 2.0) # parameters for FIR
smcrv = 2.93
smcratio = 1.6*10**(22)
smc = galaxy(smcbkg, smcfuv, smcfi, smcfii, smcfiii, smcfir, smcrv, smcratio)

def term(wl, parameters): # a single term in the extinction curve equation for optical data
    microns = wl*10**(6) # convert meter values to microns
    a = parameters.a
    b = parameters.b
    lam = parameters.lam
    n = parameters.n
    return a/((microns/lam)**n + (lam/microns)**n + b)

# add all terms of the extinction curve equation for optical data
def eta(wl, galaxy):
    return term(wl, galaxy.bkg) + term(wl, galaxy.fuv) + term(wl, galaxy.fi) + term(wl, galaxy.fii)  + term(wl, galaxy.fiii) + term(wl, galaxy.fir)

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
totsigma = (H+He+C+O+Fe+Ne+Mg+Si+S+Ar+Ca)

def sigma(e): # find sigma for a given xray energy
     dE = 0.1
     index = np.where((keV<e+dE) & (keV>e))
     s = totsigma[index[0]]
     return s[0]

def chi(model, data):
    sum = 0
    n = 0
    while n < len(data):
        sum = sum + (data[n]-model[n])**2/model[n]
        n = n + 1
    dof = len(model) - 1
    return sum/dof

fig, ax = plt.subplots(figsize=(10, 6))
norm = plt.Normalize(np.log10(1/(7500*10**(-10))), np.log10(1/(3800*10**(-10))))
ax.set(title = "Flux",
       xlabel = "log[1/λ] (m)",
       ylabel = "log[Flux]")
plt.rcParams['font.family'] = "serif"

opticalwl = optical
xraywl = (conv(xray))
wl = np.concatenate((opticalwl, xraywl), axis=0)
normalization = 10**10
initialflux = normalization*(wl**(beta))

def opticalfunc(wl, gal, z):
    v = 5500*10**(-10)
    shiftedwavelength = wl/(1+z)
    pei = (eta(shiftedwavelength, gal)/eta(v, gal))
    return (gal.rv/1.086)*pei

def xrayfunc(energy, gal, z):
    shiftedenergy = energy*(1+z)
    length = len(energy)
    s = np.zeros(length)
    for n, e in enumerate(shiftedenergy): 
        s[n] = sigma(e)
    return s*gal.ratio*gal.rv

def host(op, xr, gal, z):
    copthost = opticalfunc(op, gal, z)
    cxrayhost = xrayfunc(xr, gal, z)
    return np.concatenate((copthost, cxrayhost), axis=0)

chost = host(optical, xray, smc, zhost)

def attenhost(ebv):
    return np.e**(-1*chost*ebv)

def intervener(op, xr, gal, z):
    coptint = opticalfunc(op, gal, z)
    cxrayint = xrayfunc(xr, gal, z)
    return np.concatenate((coptint, cxrayint), axis=0)

cint = intervener(optical, xray, smc, zint)

def attenint(ebv):
    return np.e**(-1*cint*ebv)

flux = initialflux*attenhost(1)*attenint(0.5)
fluxwn = flux + np.random.normal(0, 3*np.sqrt(flux))
fluxwn[fluxwn < 0] = 0

def model(wl, ebv):
    init = normalization*(wl**(beta))
    return init*np.e**(-1*chost*ebv)

param, param_cov = curve_fit(model, wl, fluxwn)
fittedflux = initialflux*np.e**(-1*chost*param[0])

x = np.log10(1/wl)
y = np.log10(fluxwn)
plt.scatter(x, y, s=3, c=x, cmap='nipy_spectral_r', norm=norm)
plt.plot(x, np.log10(fittedflux), c = 'k', lw = 0.1)

outputchi = chi(fittedflux, fluxwn)
