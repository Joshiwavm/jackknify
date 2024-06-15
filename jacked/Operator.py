import os
import tqdm
import numpy as np

from casatasks import tclean, exportfits
from ImSettings import * 

####### Own Tools #######
from MsManager import *

class Jack:
    def __init__(self, 
                 vis, 
                 spws, 
                 fields, 
                 band,
                 array,
                 samples, 
                 test = False,
                 update_weights = False,
                 ):
        
        # initialize variables
        # --------------------
        
        self.test           = test 
        self.update_weights = update_weights
        self.seed           = 42
        
        self.fields         = fields
        self.spws           = spws
        self.N              = samples
        self.vis_input      = vis
        self.imcase         = Imparams(config=array, band=band)
        
    @property
    def manager(self):
        return MSmanager(self.vis_input, 
                         self.vis_output,
                         spws = self.spws, 
                         fields = self.fields, 
                         band = self.band, 
                         array = self.array)

    @property
    def vis_jacked(self):    
        return self.manager.ms_copydir + self.vis_output  

    @property
    def vis_output(self):
        return (self.vis_input.split('/')[-1]).split('.ms')[0] + '_Jacked_seed'+str(self.seed)+'.ms' 
        
    def _loader(self):
        uvdist, UVreal, UVimag, UVwgts, _, _ = self.manager.uvdata_loader()

        self.UVreal = UVreal
        self.UVimag = UVimag
        self.uvdist = uvdist
        self.UVwgts = UVwgts
        
    def _saver(self):
            
        if self.test is not True:
            self.manager.model_to_ms(self.UVreal_jacked + self.UVimag_jacked*1j,
                                    self.wgt)
        else:
            self.manager.model_to_ms(self.UVreal + self.UVimag*1j,
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

        self.UVreal_jacked = np.copy(self.UVreal)
        self.UVreal_jacked[indexing.astype(bool)] *= -1.
        
        self.UVimag_jacked = np.copy(self.UVimag)
        self.UVimag_jacked[indexing.astype(bool)] *= -1. 

        if not self.update_weights: self.wgt = None
        else: self._stw_manual()
    
    def _image(self):

        spw = ''
        for i,sp in enumerate(np.array(self.manager.spws).flatten()):
            if not i ==0: spw += ','+sp
            else: spw += sp
        
        tclean( vis         =                      self.vis_jacked, 
                imagename   = self.vis_jacked.replace('.ms','.im'),  
                niter       =                                    0, 
                spw         =                                  spw,
                imsize      =                self.imcase.imsize//2, 
                cell        =   str(self.imcase.cellsize)+'arcsec', 
                gridder     =                           'standard', 
                weighting   =                            'natural', 
                specmode    =                               'cube')     
        
        exportfits(imagename = self.vis_jacked.replace('.ms','.im.image'), 
                   fitsimage = self.vis_jacked.replace('.ms','.im.fits'), 
                   overwrite = True)
        
        os.system('rm -rf {0}.pb'.format(self.vis_jacked.replace('.ms','.im')))
        os.system('rm -rf {0}.psf'.format(self.vis_jacked.replace('.ms','.im')))
        os.system('rm -rf {0}.model'.format(self.vis_jacked.replace('.ms','.im')))
        os.system('rm -rf {0}.sumwt'.format(self.vis_jacked.replace('.ms','.im')))
        os.system('rm -rf {0}.weight'.format(self.vis_jacked.replace('.ms','.im')))
        os.system('rm -rf {0}.residual'.format(self.vis_jacked.replace('.ms','.im')))
        os.system('rm -rf {0}.ms'.format(self.vis_jacked))
        
    def run(self):

        for i in range(0,self.N,1):
            
            self.seed += i

            print('.. Loading in MS')
            self._loader()
            print('.. Jack Knife it')
            self.Jack_it()
            print('.. Saving to MS')
            self._saver()

            print('.. Image')
            self._image()
                    
            os.system('rm -vf *.last')
            os.system('rm -vf *.log')
