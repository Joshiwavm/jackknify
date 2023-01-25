from scipy.integrate import quad
import numpy as np

def printError(message):
    print(message)
    raise ValueError(42)

def elos(e): return e/np.sqrt(2.00-e*e)

# 3D A10 model profile
# ----------------------------------------------------------------------
def a10ProfileIntegrand(x,xi,alpha,beta,gamma,ap,c500,mass): 
    return (x**(-gamma))*((1.00+(x**alpha))**((gamma-beta)/alpha))*(x/((x*x-xi*xi)**0.50))*(mass**((ap+0.10)/(1.00+(2.00*x/c500)**3.00)))

# A10 model integrale
# ----------------------------------------------------------------------
def _a10ProfileIntegral(x,alpha,beta,gamma,ap,c500,mass,limdist=np.inf,epsrel=1.00E-06): 
    return 2.00*quad(a10ProfileIntegrand,x,limdist,args=(x,alpha,beta,gamma,ap,c500,mass),epsrel=epsrel)[0]
a10ProfileIntegral = np.vectorize(_a10ProfileIntegral)

# Integrated elliptical A10 model profile
# ----------------------------------------------------------------------
def a10Profile(grid,offset,amp,major,e,alpha,beta,gamma,ap,c500,mass, limdist=np.inf,epsrel=1.00E-06,freeLS=None): 
    integral = np.zeros_like(grid,dtype=np.float64)
    integral[grid<=limdist] = a10ProfileIntegral(grid[grid<=limdist],alpha,beta,gamma,ap,c500,mass,limdist,epsrel)
    ellipse = np.sqrt(1.00-elos(e)**2) if freeLS is None else freeLS
    return offset+amp*major*ellipse*integral




# 3D gNFW model profile
# ----------------------------------------------------------------------
def gnfwProfileIntegrand(x,xi,alpha,beta,gamma): 
    return (x**(-gamma))*((1.00+(x**alpha))**((gamma-beta)/alpha))*x/((x*x-xi*xi)**0.50)

# gNFW model integral
# ----------------------------------------------------------------------
def _gnfwProfileIntegral(x,alpha,beta,gamma,limdist=np.inf,epsrel=1.00E-06): 
    return 2.00*quad(gnfwProfileIntegrand,x,limdist,args=(x,alpha,beta,gamma),epsrel=epsrel)[0]
gnfwProfileIntegral = np.vectorize(_gnfwProfileIntegral)

# Integrated elliptical gNFW model profile
# ----------------------------------------------------------------------
def gnfwProfile(grid,offset,amp,major,e,alpha,beta,gamma,limdist=np.inf,epsrel=1.00E-06,freeLS=None): 
    integral = np.zeros_like(grid,dtype=np.float64)
    integral[grid<=limdist] = gnfwProfileIntegral(grid[grid<=limdist],alpha,beta,gamma,limdist,epsrel)
    ellipse = np.sqrt(1.00-elos(e)**2) if freeLS is None else freeLS
    return offset+amp*major*ellipse*integral




# 3D beta model profile
# ----------------------------------------------------------------------
def betaProfileIntegrand(x,xi,alpha,beta,gamma): 
    return ((1.00+(x**2.00))**(-1.50*beta))*x/((x*x-xi*xi)**0.50)

# Beta model integral
# ----------------------------------------------------------------------
def _betaProfileIntegral(x,alpha,beta,gamma,limdist=np.inf,epsrel=1.00E-06): 
    return 2.00*quad(betaProfileIntegrand,x,limdist,args=(x,beta),epsrel=epsrel)[0]
betaProfileIntegral = np.vectorize(_betaProfileIntegral)

# Integrated elliptical beta model profile 
# ----------------------------------------------------------------------
def betaProfile(grid,offset,amp,major,e,beta,limdist=np.inf,epsrel=1.00E-06,freeLS=None):
    integral = np.zeros_like(grid,dtype=np.float64)
    integral[grid<=limdist] = betaProfileIntegral(grid[grid<=limdist],beta,limdist,epsrel)
    ellipse = np.sqrt(1.00-elos(e)**2) if freeLS is None else freeLS
    return offset+amp*major*ellipse*integral