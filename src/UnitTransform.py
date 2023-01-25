from src.Utillities import *

class TransformInput():
    
    def __init__(self, 
                 popt, 
                 typ, 
                 c500=1.00, 
                 mass=1.00,
                 limdist=np.inf,
                 epsrel=1.00E-06,
                 freeLS=None):
        
        self.popt = popt
        self.type = typ
        self.c500 = c500
        self.mass = mass
        self.limdist = limdist
        self.epsrel = epsrel
        self.freeLS = freeLS
        
    def _A10Pressure(self, info):
        param_input = {
            "offset":self.popt['Offset'],
            "amp":self.popt['log10'],
            "major":self.popt['c500'],
            "e":self.popt['e'],
            "alpha":self.popt['Alpha'],
            "beta":self.popt['Beta'],
            "gamma":self.popt['Gamma'],
            "ap":self.popt['Alpha_p'],
            "c500":self.c500,
            "mass":self.mass,
            "limdist":self.limdist,
            "epsrel":self.epsrel,
            "freeLS":self.freeLS
        }
        
        Hz = info.cosmo.H(self.popt['z'])

        c500 = self.popt['c500']
        m500 = self.popt['log10']
        bias = self.popt['bias']

        r500 = ((3.00/4.00/np.pi/500.00/info.cosmo.critical_density(self.popt['z']))*(1.00-bias)*(10**m500)*u.solMass)**(1.00/3.00)

        param_input['major'] = r500.to(u.Mpc).value/c500
        param_input['amp']   = self.popt['P0']*(3.00/8.00/np.pi)*(info.fb*info.mu/info.mue)
        param_input['amp']  *= (((((2.5e2*Hz*Hz)**2.00)*((1.00-bias)*(10**(m500-15.00))*u.solMass)/(const.G**(0.50)))**(2.00/3.00)).to(u.keV/u.cm**3)).value
        param_input['amp' ] *= 1e10

        mass = ((1.00-bias)*(10**(m500-14.00))*(info.cosmo.H0.value/70.00)/3.00) 

        limdist = info.limdist*c500 if np.isfinite(info.limdist) else np.inf
        freeLS = self.popt['depth'] if self.type =='A10PressureLS' else None

        param_input['c500'] = c500
        param_input['mass'] = mass
        param_input['limdist'] = limdist
        param_input['freeLS'] = freeLS
        return param_input
    
    def _gnfwPressure(self, info):
        param_input = {
            "offset":self.popt['Offset'],
            "amp":self.popt['Amplitude'],
            "major":np.deg2rad(self.popt['Major']),
            "e":self.popt['e'],
            "alpha":self.popt['Alpha'],
            "beta":self.popt['Beta'],
            "gamma":self.popt['Gamma'],
            "limdist":self.limdist,
            "epsrel":self.epsrel,
            "freeLS":self.freeLS
        }
               
        param_input['major'] *= info.cosmo.angular_diameter_distance(self.popt['z'])        
        return param_input
    
    def _betaPressure(self, info): #done
        param_input = {
            "offset": self.popt['Offset'] ,
            "amp": self.popt['Amplitude'] ,
            "major": np.deg2rad(self.popt['Major']),
            "e": self.popt['e'] ,
            "beta": self.popt['Beta'],
            "limdist": self.limdist,
            "epsrel": self.epsrel,
            "freeLS": self.freeLS
        }
        
        param_input['major'] *= info.cosmo.angular_diameter_distance(self.popt['z']).value
        return param_input
    
    def generate(self):     
        
        if self.type.split('_')[0] == 'betaPressure': 
            return self._betaPressure(getinfo())
        
        elif self.type.split('_')[0] == 'gnfwPressure': 
            return self._gnfwPressure(getinfo())
         
        elif (self.type.split('_')[0] == 'A10Pressure') or (self.type.split('_')[0] == 'A10PressureLS'):
            return self._A10Pressure(getinfo())
