import shutil
import numpy as np
import os
import astropy.constants as const
import casatools

# The class for handling ms files
class MSmanager:
    def __init__(   self, 
                    filename_in,
                    filename_out,
                    outdir,
                    spws, 
                    fields,
                ):
      
        self.ms_file  = filename_in
        self.ms_save  = filename_out
        self.spws     = spws
        self.fields   = fields

        self.ms_copydir = outdir
        self.ms_modelfile = self.ms_copydir + self.ms_file.split('/')[-1]
        
        if not os.path.exists(self.ms_copydir):
            os.makedirs(self.ms_copydir)
        
    def uvdata_loader(self):

        UVreal = np.empty(0)
        UVimag = np.empty(0)
        uvdist = np.empty(0)
        uvwghts = np.empty(0)
        us = np.empty(0)
        vs = np.empty(0)

        
        for f, field in enumerate(self.fields):
            for s, spw in enumerate(self.spws[f]):                
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
    
    def model_to_ms(self, model, sigma = None):
        
        try: shutil.copytree(self.ms_file, self.ms_modelfile)
        except:pass
    
        ms = casatools.ms()
        ms.open(self.ms_modelfile,nomodify=False)
        
        index = 0
        for f, field in enumerate(self.fields):
            for s, spw in enumerate(self.spws[f]):

                ms.selectinit(datadescid=int(spw))
                ms.select({'field_id': int(field)})

                rec = ms.getdata(['data', 'weight'])
                freqs = ms.range('chan_freq')['chan_freq'][:,0]

                uvreal = (model[index:index+len(rec['data'][0][0]) * len(freqs)].real).reshape(len(freqs), len(rec['data'][0][0])) 
                uvimag = (model[index:index+len(rec['data'][0][0]) * len(freqs)].imag).reshape(len(freqs), len(rec['data'][0][0])) 
                
                rec['data'] = np.array([(uvreal+1.0j*uvimag),(uvreal+1.0j*uvimag)])

                if sigma != None:
                    self.wgts  = rec['weight']
                    self.wgts  = np.zeros_like(self.wgts) + 1/sigma**2
                    rec['weight'] = self.wgts
                
                ms.putdata(rec)
                ms.reset()
                
                index += len(rec['data'][0][0]) * len(freqs)
        ms.close()

        os.rename(self.ms_modelfile, self.ms_copydir + self.ms_save.split('/')[-1])
