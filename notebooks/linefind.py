import numpy as np
from interferopy.cube import Cube
import os

gamma = 3

direc = '/Users/jvanmarr/Documents/GitHub/Jack-knife/output/cubes/'
outdir = '/Users/jvanmarr/Documents/GitHub/Jack-knife/output/findclumps/'
files = os.listdir(direc)

# Filter files with .fits extension
fits_files = [file for file in files if file.endswith('.fits')]


for fname in fits_files:
    print(fname)
    print()
    
    outrun = outdir + fname.replace('.im.fits', '/')
    if not os.path.exists(outrun):
        os.mkdir(outrun)
    
    cube = Cube(direc + fname)
    _ = cube.findclumps_full(output_file=outrun + 'findlcumps',
                        kernels=np.arange(1, 12, 2),
                        SNR_min=3,                                               
                        delta_offset_arcsec=1,
                        delta_freq=0.4,

                        run_search=True,
                        run_crop=True,#True
                        run_fidelity=True, #True

                        fidelity_bins=np.arange(0, 10, 0.2),
                        min_SN_fit=3.0,
                        fidelity_threshold=0.5,

                        verbose=True,
                        sextractor_param_file='default.sex',
                        ncores=4)

#     if not os.path.exists(outrun+'findclumpoutput/'):
#         os.mkdir(outrun+'findclumpoutput/')
    
#     np.save(file = outrun+'findclumpoutput/catP_catN_CandP_CandN.npy',
#             arr  = np.array([catP, catN, candP, candN]))