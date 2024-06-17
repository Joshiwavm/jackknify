import numpy as np
from interferopy.cube import Cube

fname = '/Users/jvanmarr/Documents/GitHub/Jack-knife/data/Glass-z13_target_concat_tbin30s_cwidth38MHz_60spw.im_test.fits'
cube = Cube(fname)

catP, catN, candP, candN = cube.findclumps_full(output_file='../output/findclumps/findclumps_',
                                                kernels=np.arange(3, 16, 2),
                                                run_search=True,
                                                run_crop=True,
                                                SNR_min=3,
                                                delta_offset_arcsec=1,
                                                delta_freq=0.4,
                                                run_fidelity=True,
                                                fidelity_bins=np.arange(0, 10, 0.2),
                                                min_SN_fit=3.0, fidelity_threshold=0.5,
                                                verbose=True,
                                                sextractor_param_file = 'default.sex',
                                                ncores=4)