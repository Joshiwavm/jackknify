import scipy
from scipy.optimize import curve_fit
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse  
plt.style.use('./plots/scripts/thesis.mplstyle')
matplotlib.rc('image', origin='lower')

from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
from astropy import units as u
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.wcs import WCS

import casatools
from casatasks import *

import tqdm
####### Own Tools #######
from src.MsManager import *
from src.Settings import *
from src.UnitTransform import *

Tcmb    = 2.7255
mec2    = ((const.m_e*const.c*const.c).to(u.keV)).value
clight  = const.c.value
kboltz  = const.k_B.value
hplanck = const.h.value

class Test_SZonly:
    def __init__(self, vis, typ, model_name, seed = 195758392, spws = [['0','1','2','3']], fields = ['0'], band = 'band6', test = False):
        
        self.test       = test 
        self.seed       = seed
        self.vis        = vis
        self.vis_jacked = 'Model_'+ (vis.split('/')[-1]).split('.ms')[0] + typ + '_Jacked_seed'+str(self.seed)+'.ms'
                            
        self.model_name = model_name
        self.typ        = typ
        self.reader     = MsReader(self.vis, self.typ, spws = spws, fields = fields, band = band)
        
        
        #--- moved below from saver to here
        self.reader.ms_copydir = './output/add_ons/'
        self.reader.ms_modelfile = self.reader.ms_copydir + self.reader.ms_file.split('/')[-1]
        
        self.vis_jacked =  self.reader.ms_copydir + self.vis_jacked  
                
        if typ =='com07m':
            self.imsize = 256
            self.imcell = '1.50arcsec'
        elif typ =='com12m':
            self.imsize = 256
            self.imcell = '0.15arcsec'
            
    def _loader(self):
        uvdist, UVreal, UVimag, UVwgts, _, _ = self.reader.uvdata_loader()

        self.UVreal = UVreal
        self.UVimag = UVimag
        self.uvdist = uvdist
        self.UVwgts = UVwgts
        
    def _saver(self):
            
        if self.test is not True:
            self.reader.model_to_ms(self.UVreal_jacked + self.UVimag_jacked*1j,
                                    1., 
                                    'replace', 
                                    (self.vis.split('/')[-1]).split('.ms')[0] + self.typ + '_Jacked_seed'+str(self.seed),
                                    self.wgt)
        else:
            self.reader.model_to_ms(self.UVreal + self.UVimag*1j,
                                    1., 
                                    'replace', 
                                    (self.vis.split('/')[-1]).split('.ms')[0] + self.typ + '_Jacked_seed'+str(self.seed),
                                    self.wgt)

    def _stw_manual(self):
                
        bin_frac = 0.001
        bins = np.logspace(np.log10(np.nanmin(self.uvdist)), np.log10(np.nanmax(self.uvdist)), int(bin_frac*len(self.uvdist)))
        self.UVreal_binned = np.zeros(len(bins)-1)
        self.UVimag_binned = np.zeros(len(bins)-1)

        for i in tqdm.tqdm(range(len(bins)-1)):
            maskb =(self.uvdist >= bins[i]) & (self.uvdist < bins[i+1])
            self.UVreal_binned[i] = np.nansum(self.UVreal_jacked[maskb])
            self.UVimag_binned[i] = np.nansum(self.UVimag[maskb])
            
        self.wgt = 1/(np.nanstd(self.UVreal_binned))**2/ bin_frac
    
    def Jack_it(self):
        indexing = np.ones_like(self.UVreal)
        indexing[:len(indexing)//2] = 0

        np.random.seed(self.seed)
        np.random.shuffle(indexing)

        UVreal_jacked = np.copy(self.UVreal)
        UVreal_jacked[indexing.astype(bool)] *= -1.
        self.UVreal_jacked = UVreal_jacked
        
        UVimag_jacked = np.copy(self.UVimag)
        UVimag_jacked[indexing.astype(bool)] *= -1. 
        self.UVimag_jacked = UVimag_jacked        

        self.wgt = None
        # self._stw_manual()
        
    #######################
    ###     Imaging     ###
    #######################
    
    def _image(self): #normal
        
        #HD1 stuff
        # --------------------
        tclean(       vis     =   self.vis_jacked.replace('Model_', 'Model_0_'), 
                  imagename   =       self.vis_jacked.replace('.ms','_cube.im'),  
                  niter       =                                               0,
                  # spw         =                                            '25',
                  imsize      =                                     self.imsize, 
                  cell        =                                     self.imcell, 
                  gridder     =                                      'standard', 
                  weighting   =                                       'natural', 
                  specmode    =                                          'cube') 
        
        exportfits(self.vis_jacked.replace('.ms','_cube.im.image'), 
                   self.vis_jacked.replace('.ms','_cube.im.fits'),
                   overwrite=True)  

        image_name = self.vis_jacked.replace('.ms','_cube.im.fits')

        os.system('rm -rf {0}.pb'.format(self.vis_jacked.replace('.ms','_cube.im')))
        os.system('rm -rf {0}.psf'.format(self.vis_jacked.replace('.ms','_cube.im')))
        os.system('rm -rf {0}.model'.format(self.vis_jacked.replace('.ms','_cube.im')))
        os.system('rm -rf {0}.sumwt'.format(self.vis_jacked.replace('.ms','_cube.im')))
        os.system('rm -rf {0}.weight'.format(self.vis_jacked.replace('.ms','_cube.im')))
        os.system('rm -rf {0}.residual'.format(self.vis_jacked.replace('.ms','_cube.im')))
        os.system('rm -rf {0}.image'.format(self.vis_jacked.replace('.ms','_cube.im')))

        # Standard Stuff
        # ----------------------
        
#         tclean( vis         =                      self.vis_jacked, 
#                 imagename   = self.vis_jacked.replace('.ms','.im'),  
#                 niter       =                                    0, 
#                 imsize      =                          self.imsize, 
#                 cell        =                          self.imcell, 
#                 gridder     =                           'standard', 
#                 weighting   =                            'natural', 
#                 specmode    =                               'mfs')     
        
#         exportfits(imagename = self.vis_jacked.replace('.ms','.im.image'), 
#                    fitsimage = self.vis_jacked.replace('.ms','.im.fits'), 
#                    overwrite = True)
        
#         image_name = self.vis_jacked.replace('.ms','.im.fits')

#         os.system('rm -rf {0}.pb'.format(self.vis_jacked.replace('.ms','.im')))
#         os.system('rm -rf {0}.psf'.format(self.vis_jacked.replace('.ms','.im')))
#         os.system('rm -rf {0}.model'.format(self.vis_jacked.replace('.ms','.im')))
#         os.system('rm -rf {0}.sumwt'.format(self.vis_jacked.replace('.ms','.im')))
#         os.system('rm -rf {0}.weight'.format(self.vis_jacked.replace('.ms','.im')))
#         os.system('rm -rf {0}.residual'.format(self.vis_jacked.replace('.ms','.im')))

        noise, header = fits.getdata(image_name, header=True)
        noise = noise[0,20] *1e3

        self.noise  = noise
        self.header = header
        return noise
    

    def _getmodel(self):
        
        SZ       = fits.getdata(self.model_name)[0,0]*1e3
        Combined = self.noise + SZ
        
        self.SZ       = SZ
        return SZ, Combined
    
    def _visualize(self, im, vis_noise = False, vis_model = False, vis_SZ = False):
            
        FOV     = [int(self.imsize/4 *1.6),int(self.imsize - self.imsize/4*1.6)]

        fig, ax = plt.subplots(1,1)
        im0 = ax.imshow(im[FOV[0]:FOV[1],FOV[0]:FOV[1]])
        cbar = plt.colorbar(im0, ax = ax)
        cbar.set_label(r'$S_{\nu}$ [mJy/Beam]')
        ax.text(0.05*(FOV[1]-FOV[0]),0.90*(FOV[1]-FOV[0]), 
                '{} = {:.4f} [mJy/Beam]'.format(r"$\sigma_{S_\nu}$", np.nanstd(self.noise)), 
                color='white')
        
        if self.typ == 'com12m':
            x = np.arange(0.05*(FOV[1]-FOV[0]), 0.05*(FOV[1]-FOV[0])+20, 1)
        elif self.typ == 'com07m':
            x = np.arange(0.05*(FOV[1]-FOV[0]), 0.05*(FOV[1]-FOV[0])+2, 1)
            
        y = np.copy(x)
        y[:] = 0.05*(FOV[1]-FOV[0])

        ax.plot(x, y, '-', linewidth=3, color='white')
        ax.text(0.07*(FOV[1]-FOV[0]), 0.10*(FOV[1]-FOV[0]), r'3"', 
                verticalalignment='top', color='white', size=11)

        ax.set_xticks([])
        ax.set_yticks([])
        
        if vis_noise:
            cont_SZ = np.nanstd(self.noise)*np.array([-7,-6, -5, -4, -3, -2, 2, 3, 4])
            ax.contour(im[FOV[0]:FOV[1],FOV[0]:FOV[1]], cont_SZ)
            plt.savefig('./plots/add_ons/Jacked_knifed_obs_'+self.typ+'.pdf', dpi = 300)

        if vis_model:
            self.modelcont = np.nanstd(im)*np.array([-7,-6, -5, -4, -3, -2, 2, 3, 4])
            ax.contour(im[FOV[0]:FOV[1],FOV[0]:FOV[1]], self.modelcont)
            plt.savefig('./plots/add_ons/Veszee_model_'+self.typ+'.pdf', dpi = 300)
        
        if vis_SZ:  
            cont_SZ = np.nanstd(self.noise)*np.array([-7,-6, -5, -4, -3, -2, 2, 3, 4])
            ax.contour(self.SZ[FOV[0]:FOV[1],FOV[0]:FOV[1]], self.modelcont, colors = 'red', alpha = 0.5)
            ax.contour(im[FOV[0]:FOV[1],FOV[0]:FOV[1]], cont_SZ)
            plt.savefig('./plots/add_ons/Combined_model_jack_knifed_'+self.typ+'.pdf', dpi = 300)
            
        plt.close()
        
        
    def run(self):
        print('.. Loading in MS')
        self._loader()
        print('.. Jack Knife it')
        self.Jack_it()
        print('.. Saving to MS')
        self._saver()

        print('.. Get Jack Knifed image')
        Noise = self._image()
        print('.. Get SZ and Combine')
        SZ, Combined = self._getmodel()
        
        print()
        print('.. Visualize')
        self._visualize(Noise,    vis_noise = True,  vis_model = False, vis_SZ = False)
        self._visualize(SZ,       vis_noise = False, vis_model = True,  vis_SZ = False)
        self._visualize(Combined, vis_noise = False, vis_model = False, vis_SZ = True)
        
        os.system('rm -vf *.last')
        os.system('rm -vf *.log')
        
        return Noise, SZ, Combined