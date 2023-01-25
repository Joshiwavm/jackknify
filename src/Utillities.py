import numpy as np
import scipy

from astropy import units as u
import astropy.constants as const
from astropy.cosmology import FlatLambdaCDM

global Tcmb; Tcmb = 2.7255
global mec2; mec2 = ((const.m_e*const.c*const.c).to(u.keV)).value
global clight; clight = const.c.value
global kboltz; kboltz = const.k_B.value
global hplanck; hplanck = const.h.value

class getinfo:
      def __init__(self,
                   reffreq=1e11,
                   parline=[-5,5,100],
                   limdist=np.inf,
                   limepsr=1.00E-06,
                   fb=0.175,
                   mu=0.590,
                   mue=1.140,
                   cosmo=None):

        self.reffreq = reffreq
        self.limdist = limdist
        self.limepsr = limepsr

        self.linmesh = np.array([])

        if len(parline):
            if len(parline)!=3: 
                printError('Provide correct numbers of line parameters [log scale]')
            
        self.linmesh = np.append(0.0,np.logspace(parline[0],parline[1],parline[2],dtype=np.float64))

        self.ysznorm = const.sigma_T/const.m_e/const.c**2
        self.ysznorm = self.ysznorm.to(u.cm**3/u.keV/u.Mpc)

        self.cosmo = FlatLambdaCDM(H0=70.00,Om0=0.30) if cosmo is None else cosmo

        self.fb = fb
        self.mu = mu
        self.mue = mue
        
# Adimensional frequency
# ----------------------------------------------------------------------
def getx(freq):
    factor = const.h*freq*u.Hz/const.k_B/(Tcmb*u.Kelvin)
    return factor.to(u.dimensionless_unscaled).value
        
        
# CMB surface brightness
# ----------------------------------------------------------------------
def getJynorm():
    factor  = 2e26
    factor *= (const.k_B*Tcmb*u.Kelvin)**3 # (kboltz*Tcmb)**3.0
    factor /= (const.h*const.c)**2         # (hplanck*clight)**2.0
    return factor.value


# Compton y to Jy/pixel
# ----------------------------------------------------------------------
def comptonToJyPix(freq,ipix,jpix):
    x = getx(freq)
    factor  = getJynorm()
    factor *= -4.0+x/np.tanh(0.5*x)
    factor *= (x**4)*np.exp(x)/(np.expm1(x)**2)
    factor *= np.abs(ipix*jpix)*(np.pi/1.8e2)*(np.pi/1.8e2)
    return factor

def comptonRelativ(freq,order=1):
    x = getx(freq)
    xt = x/np.tanh(0.5*x)
    st = x/np.sinh(0.5*x)
    Y0 = -4.0+xt
    
    if (order==0):
        return 1.0
    if (order==1):
        Y1 = -10.+((47./2.)+(-(42./5.)+(7./10.)*xt)*xt)*xt+st*st*(-(21./5.)+(7./5.)*xt)
        return np.divide(Y1,Y0)
    if (order==2):
        Y2 = (-15./2.)+((1023./8.)+((-868./5.)+((329./5.)+((-44./5.)+(11./30.)*xt)*xt)*xt)*xt)*xt+ \
             ((-434./5.)+((658./5.)+((-242./5.)+(143./30.)*xt)*xt)*xt+(-(44./5.)+(187./60.)*xt)*(st*st))*st*st
        return np.divide(Y2,Y0)
    if (order==3):
        Y3 = (15./2.)+((2505./8.)+((-7098./5.)+((14253./10.)+((-18594./35.)+((12059./140.)+((-128./21.)+(16./105.)*xt)*xt)*xt)*xt)*xt)*xt)*xt+ \
             (((-7098./10.)+((14253./5.)+((-102267./35.)+((156767./140.)+((-1216./7.)+(64./7.)*xt)*xt)*xt)*xt)*xt) +
             (((-18594./35.)+((205003./280.)+((-1920./7.)+(1024./35.)*xt)*xt)*xt) +((-544./21.)+(992./105.)*xt)*st*st)*st*st)*st*st
        return np.divide(Y3,Y0)
    if (order==4):
        Y4 = (-135./32.)+((30375./128.)+((-62391./10.)+((614727./40.)+((-124389./10.)+((355703./80.)+((-16568./21.)+((7516./105.)+((-22./7.)+(11./210.)*xt)*xt)*xt)*xt)*xt)*xt)*xt)*xt)*xt + \
             ((-62391./20.)+((614727./20.)+((-1368279./20.)+((4624139./80.)+((-157396./7.)+((30064./7.)+((-2717./7.)+(2761./210.)*xt)*xt)*xt)*xt)*xt)*xt)*xt + \
             ((-124389./10.)+((6046951./160.)+((-248520./7.)+((481024./35.)+((-15972./7.)+(18689./140.)*xt)*xt)*xt)*xt)*xt +\
             ((-70414./21.)+((465992./105.)+((-11792./7.)+(19778./105.)*xt)*xt)*xt+((-682./7.)+(7601./210.)*xt)*st*st)*st*st)*st*st)*st*st
        return np.divide(Y4,Y0)
        
# Relativistic corrections
def comptonCorrect(y, Te=0.0,limsize=np.inf):
    if hasattr(Te,'__len__') and (Te.size>=limsize):
        local_dict = {'t': ne.evaluate('Te/mec2')}
        for i in range(len(y)):
            local_dict.update({'y{0}'.format(i): y[i]})
        
        if (not np.shape(y)): return y
        elif y.ndim==1:       return np.full(np.shape(Te),y[0])
        elif y.ndim==2:       return ne.evaluate('y0+y1*t',local_dict=local_dict)
        elif y.ndim==3:       return ne.evaluate('y0+y1*t+y2*(t**2)',local_dict=local_dict)
        elif y.ndim==4:       return ne.evaluate('y0+y1*t+y2*(t**2)+y3*(t**3)',local_dict=local_dict)
        elif y.ndim==5:       return ne.evaluate('y0+y1*t+y2*(t**2)+y3*(t**3)+y4*(t**4)',local_dict=local_dict)
    else:
        ycorr = 0.0
        if (not np.shape(y)): y = [y]
        for i in range(len(y)-1,0,-1): ycorr = (y[i]+ycorr)*Te/mec2
        return ycorr+y[0]
    
def yszCorrect(freq,cdelt,order):
    return comptonToJyPix(freq,cdelt[0],cdelt[1])* comptonRelativ(freq,order)

# Provide SZ relativistic correction terms
def computeFlatCompton(freq,cdelt,order):
    if (freq[0]!=freq[1]):
        return scipy.integrate.quad(yszCorrect,freq[0],freq[1],args=([cdelt[0],cdelt[1]],order))[0]/(freq[1]-freq[0])
    else: 
        return yszCorrect(freq[0],cdelt,order)
