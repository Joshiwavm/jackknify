from src.Utillities import *

from astropy.constants import c
from astropy.io import fits
import sys
import shutil
import os

import casatools
from casatasks import *

class MsReader:
    def __init__(self, filename, obs_type, spws = [['0','1','2','3']], fields = ['0'], band = 'band3'): #quickfix hardcoded spws and fields etc.
      
        self.ms_file  = filename
        self.obs_type = obs_type
        self.band     = band
        self.fields   = fields
        self.spws     = spws
        self._set_obs()
        
        self.binvis     = './output/{1}/output_{0}_{1}.ms.field-fid.spw-sid'.format(band,self.obs_type)
        self.ms_copydir = './output/{1}/model/'.format(band,self.obs_type)
        self.timebin    = '0s'
        
        self.ms_modelfile = self.ms_copydir + self.ms_file.split('/')[-1]

    def _set_obs(self):
        if self.obs_type=='com07m':
            self.imsize = 256
            self.imcell = '1.50arcsec'
        elif self.obs_type=='com12m':
            self.imsize = 1024
            self.imcell = '0.15arcsec'
        elif self.obs_type=='ext12m':
            self.imsize = 512
            self.imcell = '0.10arcsec'

        else:
            raise print("Wrong array element")
        
    def prep_pb(self):
        
        ms=casatools.ms()
        self.splitms_file = []
        
        print('Pre-Processing Primairy Beam')
        for f, field in enumerate(self.fields):
            print('- Processing field {0}'.format(field))
            for s, spw in enumerate(self.spws[f]):
                print('-- Spectral window {0}'.format(spw))

                outvis = self.binvis.replace('-fid','-{0}'.format(field))
                outvis = outvis.replace('-sid','-{0}'.format(spw))

                tclean(vis         =                  self.ms_file,
                       imagename   =   outvis.replace('.ms','.im'),
                       datacolumn  =                        'data',
                       field       =                         field,
                       spw         =                           spw,
                       niter       =                             0,
                       pblimit     =                           0.0,
                       imsize      =                   self.imsize,
                       cell        =                   self.imcell,
                       gridder     =                    'standard',
                       weighting   =                     'natural',
                       specmode    =                         'mfs',
                       parallel    =                         False)
                
                exportfits('{0}.pb'.format(outvis.replace('.ms','.im')),
                           '{0}.pbeam.fits'.format(outvis.replace('.ms','.im')),
                           overwrite=True)
            
                os.system('rm -rf {0}.pb'.format(outvis.replace('.ms','.im')))
                os.system('rm -rf {0}.psf'.format(outvis.replace('.ms','.im')))

                os.system('rm -rf {0}.image'.format(outvis.replace('.ms','.im')))
                os.system('rm -rf {0}.model'.format(outvis.replace('.ms','.im')))
                os.system('rm -rf {0}.sumwt'.format(outvis.replace('.ms','.im')))
                os.system('rm -rf {0}.weight'.format(outvis.replace('.ms','.im')))
                os.system('rm -rf {0}.residual'.format(outvis.replace('.ms','.im')))

            
                hdu = fits.open('{0}.pbeam.fits'.format(outvis.replace('.ms','.im')))[0]
                hdu.data[np.isnan(hdu.data)] = 0.0

                data = np.copy(hdu.data[0,0])
                diff0a = np.ones(np.shape(data))
                diff0b = np.ones(np.shape(data))
                diff1a = np.ones(np.shape(data))
                diff1b = np.ones(np.shape(data))

                diff0a[1:,   :  ] = np.diff(data,axis=0); diff0a[int(diff0a.shape[0]/2)+1:,:] *= -1.0 
                diff0b[ :-1, :  ] = np.diff(data,axis=0); diff0b[int(diff0b.shape[0]/2):,:]   *= -1.0 
                diff1a[ :  ,1:  ] = np.diff(data,axis=1); diff1a[:,int(diff1a.shape[1]/2)+1:] *= -1.0 
                diff1b[ :  , :-1] = np.diff(data,axis=1); diff1b[:,int(diff1b.shape[1]/2):]   *= -1.0 

                mask0 = np.logical_or(diff0a<0,diff0b<0)
                mask1 = np.logical_or(diff1a<0,diff1b<0)
                mask  = np.logical_or(mask0,mask1)
                mask  = np.logical_or(mask,data==0)
                data[mask] = np.nan

                hdu.data[0,0] = np.copy(data)
                hdu.writeto('{0}.pbeam.fits'.format(outvis.replace('.ms','.im')),overwrite=True)

        os.system('rm -vf *.last')
        os.system('rm -vf *.log')
        print()

    def uvdata_loader(self):

        UVreal = np.empty(0)
        UVimag = np.empty(0)
        uvdist = np.empty(0)
        uvwghts = np.empty(0)
        us = np.empty(0)
        vs = np.empty(0)

        
        for f, field in enumerate(self.fields):
            print('- Processing field {0}'.format(field))

            for s, spw in enumerate(self.spws[f]):
                print('-- Spectral window {0}'.format(spw))
                
                ms=casatools.ms()
                ms.open(self.ms_file)
                ms.selectinit(reset=True)
                ms.selectinit(datadescid=int(spw))
                ms.select({'field_id': int(field)})
                
                rec = ms.getdata(['u','v','data','weight'])
                uvreal = ((rec['data'][0][:].real+rec['data'][1][:].real)/2.0)
                uvimag = ((rec['data'][0][:].imag+rec['data'][1][:].imag)/2.0)
                uvwght = 4.0/(1.0/rec['weight'][0]+1.0/rec['weight'][1])
                
                u = rec['u']
                v = rec['v']
                freqs  = ms.range('chan_freq')['chan_freq'][:,0]                
                
                ms.close()
            
                uwave  = (u.reshape(-1,1)*freqs/const.c.value)
                vwave  = (v.reshape(-1,1)*freqs/const.c.value)
                
                uwave  = np.swapaxes(uwave, 0, 1)
                vwave  = np.swapaxes(vwave, 0, 1)
                shapes  = np.ones_like(uwave)

                uwave  = uwave.flatten()
                vwave  = vwave.flatten()
                                
                uvwght = (shapes*uvwght.reshape(1,-1)).flatten()
              
                uvdist  = np.append(uvdist, (uwave**2 + vwave**2)**0.5*1e-3)
                UVreal  = np.append(UVreal,  uvreal)
                UVimag  = np.append(UVimag,  uvimag)
                uvwghts = np.append(uvwghts, uvwght)
                us      = np.append(us, u) 
                vs      = np.append(vs, v)

        return uvdist, UVreal, UVimag, uvwghts, us, vs
    
    def uvloader(self, spw, field):

        ms=casatools.ms()
        ms.open(self.ms_file)
        ms.selectinit(reset=True)
        ms.selectinit(datadescid=int(spw))
        ms.select({'field_id': int(field)})

        u = np.copy(ms.getdata(['u'])['u'])
        v = np.copy(ms.getdata(['v'])['v'])
        freqs  = ms.range('chan_freq')['chan_freq'][:,0]
        ms.close()

        uwave  = (freqs * u.reshape(-1,1)/const.c.value).flatten()
        vwave  = (freqs * v.reshape(-1,1)/const.c.value).flatten()
        uvfreq = (np.ones((len(u), len(freqs)))*freqs).flatten() 
            
        uvdata = np.array([uwave, vwave, uvfreq])
        return uvdata        
    
    def model_to_ms(self, model, scale, todo, savename, sigma = None, iters = 0, mock = True):
        
        try: shutil.copytree(self.ms_file, self.ms_modelfile)
        except:pass
    
        ms = casatools.ms()
        ms.open(self.ms_modelfile,nomodify=False)
        
        index = 0
        for f, field in enumerate(self.fields):
            print('- Processing field {0}'.format(field))
            for s, spw in enumerate(self.spws[f]):
                print('-- Spectral window {0}'.format(spw))

                ms.selectinit(datadescid=int(spw))
                ms.select({'field_id': int(field)})

                rec = ms.getdata(['data', 'weight'])
                freqs = ms.range('chan_freq')['chan_freq'][:,0]
                       
                if not mock:
                    uvreal = (model[index:index+len(rec['data'][0][0]) * len(freqs)].real).reshape(len(rec['data'][0][0]), len(freqs)) 
                    uvimag = (model[index:index+len(rec['data'][0][0]) * len(freqs)].imag).reshape(len(rec['data'][0][0]), len(freqs)) 
                    uvreal = np.swapaxes(uvreal, 0,1) #/ 1e4 #Quicckkkfix 
                    uvimag = np.swapaxes(uvimag, 0,1) #/ 1e4 #Quicckkkfix
                else:
                    uvreal = (model[index:index+len(rec['data'][0][0]) * len(freqs)].real).reshape(len(freqs), len(rec['data'][0][0])) 
                    uvimag = (model[index:index+len(rec['data'][0][0]) * len(freqs)].imag).reshape(len(freqs), len(rec['data'][0][0])) 

                if sigma != None:
                    self.wgts  = rec['weight']
                    self.wgts  = np.zeros_like(self.wgts) + 1/sigma**2
                    if   (todo== 'replace'): 
                        rec['data'] = scale*np.array([(uvreal+1.0j*uvimag),(uvreal+1.0j*uvimag)])
                        rec['weight'] = self.wgts
                    elif (todo=='subtract'): 
                        rec['data'] = np.copy(rec['data']) - scale*np.array([(uvreal+1.0j*uvimag),(uvreal+1.0j*uvimag)]) 
                        rec['weight'] = self.wgts
                else:
                    if   (todo== 'replace'): 
                        rec['data'] = scale*np.array([(uvreal+1.0j*uvimag),(uvreal+1.0j*uvimag)])

                    elif (todo=='subtract'): 
                        rec['data'] = np.copy(rec['data']) - scale*np.array([(uvreal+1.0j*uvimag),(uvreal+1.0j*uvimag)]) 
                
                ms.putdata(rec)
                ms.reset()
                
                index += len(rec['data'][0][0]) * len(freqs)
        ms.close()
        
        if  (todo== 'replace'): 
            os.rename(self.ms_modelfile, self.ms_copydir + 'Model_'+ str(iters)+'_' + savename + '.ms')
        elif (todo=='subtract'): 
            os.rename(self.ms_modelfile, self.ms_copydir + 'Model-residue_'+ str(iters)+'_' + savename + '.ms')        