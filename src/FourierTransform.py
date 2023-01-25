from src.Utillities import *
from src.Settings import *

import scipy.integrate as scintegr
import scipy.interpolate as scinterp
from galario.double import sampleImage

class FT():
    def __init__(self, uvdata, pb_head, pb_im, info):
        self.uvdata   = uvdata
        self.pb_head  = pb_head
        self.pb_im    = pb_im
        self.info     = info
        
    def _uvpoint(self, dx, dy, u, v, off, amp):
        return off+amp*np.exp(2.00*np.pi*1j*(u*dx+v*dy))  
    
    def _uvgauss(self, dx, dy, u, v, off, amp, scale, e, theta):
        sint = np.sin(theta*np.pi/1.8e2)
        cost = np.cos(theta*np.pi/1.8e2)

        ur = np.pi*(-u*sint-v*cost)*(scale*np.pi/1.8e2)
        vr = np.pi*( u*cost-v*sint)*(scale*np.pi/1.8e2)*np.sqrt(1.00-e**2)

        factor = amp*np.exp(-2.00*(ur**2+vr**2))
        return off+factor*np.exp(2.00*np.pi*1j*(u*dx+v*dy))
    
    ################################################################################
        
    def _spw_SimpleSource(self, model):
        spwRA   = self.pb_head['CRPIX1'] 
        spwRA  += (model['RA']-self.pb_head['CRVAL1'])*np.cos(np.deg2rad(self.pb_head['CRVAL2']))/self.pb_head['CDELT1'] 
        spwDec  = self.pb_head['CRPIX2']+(model['Dec']-self.pb_head['CRVAL2'])/self.pb_head['CDELT2']
        spwairy = scinterp.RectBivariateSpline(np.linspace(0,self.pb_head['NAXIS2'],self.pb_head['NAXIS2']),
                                               np.linspace(0,self.pb_head['NAXIS1'],self.pb_head['NAXIS1']),
                                               self.pb_im,kx=1,ky=1,s=0)
        
        spwbeam = spwairy.ev(spwDec,spwRA)
        return spwbeam

        
    def _spw_powerLaw(self, model, spectrum, reffreq, spwbeam):
        if list(model.keys())[2] == 'Amplitude':
            spwpoint  = spwbeam*model['Amplitude']
        else:
            spwpoint  = spwbeam*model['log10']
            
        index     = spectrum['SpecIndex']
        spwpoint *= ((self.uvdata[2]/reffreq)**(index))+0j
        return spwpoint
    
    def _spw_powerLawMod(self, model, spectrum, reffreq, spwbeam):
        if list(model.keys())[2] == 'Amplitude':
            spwpoint  = spwbeam*model['Amplitude']
        else:
            spwpoint  = spwbeam*model['log10']

        index     = spectrum['SpecIndex']+spectrum['SpecCurv']*np.log(self.uvdata[2]/reffreq)
        spwpoint *= ((self.uvdata[2]/reffreq)**(index))+0j
        return spwpoint
    
    def _spw_powerDust(self, model, spectrum, reffreq, spwbeam): #Haven't tested this yet
        frqfactor = 0.50*(self.uvdata[2][0]+self.uvdata[2][1])
        
        spwfactor  = (1.00+spectrum['z'])*((spectrum['kappa0']*u.m**2/u.kg).to(u.Mpc**2/u.solMass))
        spwfactor /= self.info.cosmo.luminosity_distance(spectrum['z'])**2
        spwfactor *= (10**spectrum['SpecCurv'])*u.solMass
        spwfactor *= ((1.00+spectrum['z'])*frqfactor/spectrum['nu0'])**spectrum['beta'] 
        
        spwfactor *= 2.00*const.h*(((1.00+spectrum['z'])*frqfactor*u.Hz)**3)/const.c**2
        spwfactor /= np.expm1((const.h*(1.00+spectrum['z'][4])*frqfactor*u.Hz/const.k_B/spectrum['Temp']/u.Kelvin).to(u.dimensionless_unscaled).value)

        if list(model.keys())[2] == 'Amplitude':
            spwfactor  = ((frqfactor/self.info.reffreq)**spectrum['SpecIndex'])+spwfactor.to(u.Jy).value/model['Amplitude']
            spwpoint = spwfactor*spwbeam*model['Amplitude'] 
        else:
            spwfactor  = ((frqfactor/self.info.reffreq)**spectrum['SpecIndex'])+spwfactor.to(u.Jy).value/model['log10']
            spwpoint = spwfactor*spwbeam*model['log10'] 
            
        return spwpoint
    
    def _spw_tSZ(self, model, osz = 4):
        freqs = np.unique(self.uvdata[2])
        freq1 = freqs[0]
        freq2 = freqs[-1]
        
        cdelt1 = np.abs(self.pb_head[list(filter(lambda x: x in self.pb_head,['CDELT1','CD1_1']))[0]])
        cdelt2 = np.abs(self.pb_head[list(filter(lambda x: x in self.pb_head,['CDELT2','CD2_2']))[0]])

        spwconv = np.array([computeFlatCompton([freq1, freq2], [cdelt1, cdelt2], order) for order in range(1+osz)])
        return comptonCorrect(spwconv, model['Temperature'])
    
    def _spw_correct(self, model, spectrum, spectrum_type, parttouv_check = True):

        spectrum_type = spectrum_type.split('_')[0]
        reffreq = 1e11 #QUICKFIX
        
        if spectrum_type == 'powerLaw':
            if parttouv_check: spwbeam = self._spw_SimpleSource(model) 
            else:              spwbeam = 1.
            spwpoint = self._spw_powerLaw(model, spectrum, reffreq, spwbeam)
        elif spectrum_type == 'powerLawMod':
            if parttouv_check: spwbeam = self._spw_SimpleSource(model) 
            else:              spwbeam = 1.     
            spwpoint = self._spw_powerLawMod(model, spectrum, reffreq, spwbeam)    
        elif spectrum_type == 'powerDust':
            if parttouv_check: spwbeam = self._spw_SimpleSource(model) 
            else:              spwbeam = 1.     
            spwpoint = self._spw_powerLawDust(model, spectrum, reffreq, spwbeam)    
        elif spectrum_type == 'tSZ':
            spwpoint = self._spw_tSZ(model)  

        return spwpoint

    
    ################################################################################
        
    def partouv(self, spwpoints, model, spectrum, model_type, spectrum_type):
        
        spwpoint = self._spw_correct(model, spectrum, spectrum_type)
        model_type    = model_type.split('_')[0]  
        
        if model_type == 'pointSource':
            spwpoints += self._uvpoint(dx  = (np.deg2rad(model['RA']-self.pb_head['CRVAL1']))*np.cos(np.deg2rad(self.pb_head['CRVAL2'])), 
                                       dy  = np.deg2rad(model['Dec']-self.pb_head['CRVAL2']),
                                       u   = self.uvdata[0], 
                                       v   = self.uvdata[1],
                                       off = model['Offset'],
                                       amp = spwpoint)

            
        elif model_type == 'gaussSource':
            spwpoints += self._uvgauss( dx    = np.deg2rad(model['RA']-self.pb_head['CRVAL1'])*np.cos(np.deg2rad(self.pb_head['CRVAL2'])), 
                                        dy    = np.deg2rad(model['Dec']-self.pb_head['CRVAL2']),
                                        u     = self.uvdata[0], 
                                        v     = self.uvdata[1],
                                        off   = model['Offset'],
                                        amp   = spwpoint,
                                        scale = model['Major'],
                                        e     = model['e'],
                                        theta = model['Angle'])
        
        return spwpoints
    
    ################################################################################
    
    def _sztouv(self, imdata, model, spectrum, spectrum_type, Tout, osz, iscompton):
        
        freqs = np.unique(self.uvdata[2])
        if not iscompton: 
            imfreq = [freqs[0], freqs[-1]]
            imdelt = [np.abs(self.pb_head['CDELT1']),np.abs(self.pb_head['CDELT2'])]
            imconv = comptonCorrect(np.array([computeFlatCompton(imfreq,imdelt,order) for order in range(1+osz)]),Tout)
            imdata = imdata/imconv

        uvimage = np.ascontiguousarray(np.flip(np.multiply(imdata,self.pb_im),axis=0))
        uvmodel = sampleImage(uvimage, self.pb_head['CDELT2']*np.pi/180., self.uvdata[0], self.uvdata[1])

        uvconv  = comptonCorrect(np.array([comptonToJyPix(self.uvdata[2],
                                                          self.pb_head['CDELT1'],
                                                          self.pb_head['CDELT2'])*comptonRelativ(self.uvdata[2],order) for order in range(1+osz)]
                                         ),
                                 Tout)
        
        uvmodel = np.multiply(uvconv, uvmodel) 

        return uvmodel
    
    def imtouv(self, imdata, model, spectrum, model_type, spectrum_type, Tout=0.00, osz=4, iscompton=True): #QuickFix hardcoded iscompt

        model_type    = model_type.split('_')[0]  
        if COMPONENTS[model_type]['SZ'] == True: 
            return self._sztouv(imdata, model, spectrum, spectrum_type, Tout, osz, iscompton)
        else:
            print("Nee nee")
            pass
            