import numpy as np
from interferopy.cube import Cube
import os

direc = '/Users/jvanmarr/Documents/GitHub/Jack-knife/output/cubes/'
outdir = '/Users/jvanmarr/Documents/GitHub/Jack-knife/output/findclumps/'
files = os.listdir(direc)

fits_files = [file for file in files if file.endswith('.fits')]

for fname in fits_files:
    print(fname)
    print()
    
    outrun = outdir + fname.replace('.im.fits', '/')
    if not os.path.exists(outrun):
        os.mkdir(outrun)
    
    cube = Cube(direc + fname)
    _ = cube.findclumps_full(output_file=outrun + 'findlcumps',
                        kernels=np.arange(5, 22, 2),
                        SNR_min=0,                                               
                        delta_offset_arcsec=0.2, 
                        delta_freq=0.1,

                        run_search=True,
                        run_crop=True,
                        run_fidelity=False, #True

                        verbose=True,
                        sextractor_param_file='default.sex',
                        ncores=4)
