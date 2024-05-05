import numpy as np
import matplotlib.pyplot as plt

# wavelengths
wavelength = np.arange((10**(-3)), (10**(2)), (10**(-3))) # range of wavelengths to plot curve for
frac = 1/wavelength

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

# LMG parameters
lmgbkg = parameters(175.0, 90.0, 0.046, 2.0) # parameters for BKG
lmgfuv = parameters(19.0, 5.5, 0.08, 4.5) # parameters for FUV
lmgfi = parameters(0.023, -1.95, 0.22, 2.0) # parameters for 2175 angstrom feature
lmgfii = parameters(0.005, -1.95, 9.7, 2.0) # parameters for 9.7 micrometer feature
lmgfiii = parameters(0.006, -1.8, 18.0, 2.0) # parameters for 18 micrometer feature
lmgfir = parameters(0.02, 0.0, 25.0, 2.0) # parameters for FIR


# SMG parameters
smgbkg = parameters(185.0, 90.0, 0.042, 2.0) # parameters for BKG
smgfuv = parameters(27.0, 5.5, 0.08, 4.0) # parameters for FUV
smgfi = parameters(0.005, -1.95, 0.22, 2.0) # parameters for 2175 angstrom feature
smgfii = parameters(0.01, -1.95, 9.7, 2.0) # parameters for 9.7 micrometer feature
smgfiii = parameters(0.012, -1.8, 18.0, 2.0) # parameters for 18 micrometer feature
smgfir = parameters(0.03, 0.0, 25.0, 2.0) # parameters for FIR

def term(wavelength, parameters): # a single term in the extinction curve equation
    a = parameters.a
    b = parameters.b
    lam = parameters.lam
    n = parameters.n
    return a/((wavelength/lam)**n + (lam/wavelength)**n + b)

# add all terms of the extinction curve equation
def xi(wavelength):
    return term(wavelength, bkg) + term(wavelength, fuv) + term(wavelength, fi) + term(wavelength, fii)  + term(wavelength, fiii) + term(wavelength, fir)

def xilmg(wavelength):
    return term(wavelength, lmgbkg) + term(wavelength, lmgfuv) + term(wavelength, lmgfi) + term(wavelength, lmgfii)  + term(wavelength, lmgfiii) + term(wavelength, lmgfir)

def xismg(wavelength):
    return term(wavelength, smgbkg) + term(wavelength, smgfuv) + term(wavelength, smgfi) + term(wavelength, smgfii)  + term(wavelength, smgfiii) + term(wavelength, smgfir)

# plotting
fig, ax = plt.subplots(figsize=(8, 6))
x = np.log(frac)

# plot milky way curve
ymw = np.log(xi(wavelength))
plt.plot(x, ymw,
        color = 'black',
        label="MW",
        linewidth=1)

# plot lmg curve
ylmg = np.log(xilmg(wavelength))
plt.plot(x, ylmg,
        color = 'darkslateblue',
        label="LMG",
        linewidth=1)

# plot smg curve
ysmg = np.log(xismg(wavelength))
plt.plot(x, ysmg,
        color = 'cornflowerblue',
        label="SMG",
        linewidth=1)

# axis labels and formatting
plt.legend(loc="best", edgecolor = '0')
ax.set(title = "Extinction Curves for Milky Way, LMG and SMG",
       xlabel = "log[1/λ] (1/μm)",
       ylabel = "log[ξ(λ)]")
plt.rcParams['font.family'] = "serif"

# show plot
plt.show()
