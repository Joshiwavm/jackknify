from casatasks import tclean, exportfits
import os
import tqdm
from astropy.io import fits
from IPython.display import clear_output

from .ImSettings import * 
from .MsManager import MSmanager
from .Plot import SLP, IMAGE

from . import JAXknife

class Jack:
    def __init__(self, 
                 fname: str, 
                 outdir: str,
                 spws: list, 
                 fields: list, 
                 band: str,
                 array: str,
                 test: bool = False,
                 update_weights: bool = False,
                 ):
        
        # initialize variables
        # --------------------
        self.test           = test 
        self.update_weights = update_weights
        self.band           = band
        self.array          = array
        self.seed           = 42
        self.fields         = fields
        self.spws           = spws
        self.vis_input      = fname
        self.outdir         = outdir

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
                
    @property
    def manager(self):
        return MSmanager(filename_in = self.vis_input, 
                         filename_out = self.vis_output,
                         outdir = self.outdir+'ms_files/',
                         spws = self.spws, 
                         fields = self.fields, 
                        )
    
    @property
    def imcase(self):
        return Imparams(config=self.array, band=self.band)

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
    
    def _jack_it(self):
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
    
    def clean(self,
               vis: str,
               reffreq:str = ''
               ):

        self.outfile = '{0}{1}'.format(self.outdir,'cubes/')

        if not os.path.exists(self.outfile):
            os.makedirs(self.outfile)
        
        self.outfile = self.outfile + '{0}'.format(vis.split('/')[-1])
        self.outfile = self.outfile.replace('.ms','.im')

        spw = ''
        for i,sp in enumerate(np.array(self.manager.spws).flatten()):
            if not i ==0: spw += ','+sp
            else: spw += sp

        tclean( vis         =                                  vis, 
                imagename   =                         self.outfile,  
                niter       =                                    0, 
                spw         =                                  spw,
                imsize      =                   self.imcase.imsize, 
                cell        =   str(self.imcase.cellsize)+'arcsec', 
                gridder     =                           'standard', 
                weighting   =                            'natural', 
                reffreq     =                              reffreq,
                specmode    =                               'cube')     
        
        exportfits(imagename = self.outfile + '.image', 
                   fitsimage = self.outfile + '.fits', 
                   overwrite = True)
        
        os.system('rm -rf {0}.pb'.format(self.outfile))
        os.system('rm -rf {0}.psf'.format(self.outfile))
        os.system('rm -rf {0}.model'.format(self.outfile))
        os.system('rm -rf {0}.weight'.format(self.outfile))
        os.system('rm -rf {0}.residual'.format(self.outfile))
        os.system('rm -rf {0}.image'.format(self.outfile))
        os.system('rm -rf {0}.sumwt'.format(self.outfile))
    
    def run(self, 
            samples: int = 1,
            seed: int = 42
            ):
        
        for i in tqdm.tqdm(range(0,samples,1)):

            self.seed = i + seed

            print('.. Loading in MS')
            self._loader()
            print('.. Jack Knife it')
            self._jack_it()
            print('.. Saving to MS')
            self._saver()

            print('.. Image')
            self.clean(self.vis_jacked)
                    
            os.system('rm -vf *.last')
            os.system('rm -vf *.log')

            if i != samples-1:
                clear_output(True)

    def gofast(self,
               samples: int = 1,
               seed: int = 42):
        print('.. Loading in MS')
        _, UVreal, UVimag, UVwgts, uw, vw = self.manager.uvdata_loader()
        
        rng = np.random.default_rng(seed=seed)

        cdelt = np.deg2rad(self.imcase.cellsize/3.60E+03)
        xw = -2.00*np.pi*vw*cdelt
        yw =  2.00*np.pi*uw*cdelt

        print('.. Running JAXknife')
        for i in tqdm.tqdm(range(0,samples,1)):
            img = JAXknife.run(xw,yw,UVreal+1j*UVimag,UVwgts,self.imcase.imsize,rng)

            os.system(f'mkdir -p {self.outdir}jaxknife/')
            outpath = self.vis_jacked.split('/')[-1]
            outpath = f'{self.outdir}jaxknife/{outpath}_jax_samp_{i:05d}'

            np.savez_compressed(outpath,img)

            if i != samples-1:
                clear_output(True)

    def plot_map(self, 
                 savedir:str = './',
                 moment: str = 'continuum',
                 channels:int = None,
                 center_coord: tuple = (None, None), 
                 box_size_arcsec: float = None,
                 idx:int = 0,
                 ):
        
        if not os.path.exists(savedir):
            os.makedirs(savedir)
            
        if not hasattr(self, 'outfile'):  
            raise RuntimeError("Please provide a fits file in self.outfile, or simply run clean.")
        else:
            fname = self.outfile + '.fits'
        
        obj = IMAGE(fname = fname,
                    moment = moment, 
                    channels = channels, 
                    center_coord = center_coord, 
                    idx = idx,
                    box_size_arcsec = box_size_arcsec)
        
        obj.plot(savedir = savedir)

    def plot_slp(self,
                 savedir:str = './',
                 idx:int = 0,
                 size:float = 1.,
                 savenp:bool=False,
                ):

        if not os.path.exists(savedir):
            os.makedirs(savedir)

        if not hasattr(self, 'outfile'):
            raise RuntimeError("Please provide a fits file in self.outfile, or simply run clean.")
        else:
            fname = self.outfile + '.fits'

        obj = SLP(fname, size, idx, amount = 2000, visualize=False)
        obj.plot(savedir = savedir, savenp = savenp)