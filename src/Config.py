# file for hard coded stuff like imsize, cell size per band and configuration
import numpy as np

BINARY_VALUES = np.array([128, 256, 512, 1024, 2048])

def FOV(freq, D):
    return 2.44 * np.rad2deg(299792458/freq/D) * 3600 #arcsec

def imsize(PB, cell_size):
    ideal = PB/cell_size
    bit_value = BINARY_VALUES[np.argmin(abs(BINARY_VALUES - ideal))]
    return bit_value
    

RESOLUTIONS = {
    ('7m', 'Band1'): ...,
    ('7m', 'Band3'): 12.5,
    ('7m', 'Band4'): 8.35,
    ('7m', 'Band5'): 6.77,
    ('7m', 'Band6'): 5.45,
    ('7m', 'Band7'): 3.63,
    ('7m', 'Band8'): 2.72,
    ('7m', 'Band9'): 1.93,
    ('7m', 'Band10'):1.44,

    ('C1', 'Band1'): ...,
    ('C1', 'Band3'): 3.38,
    ('C1', 'Band4'): 2.25,
    ('C1', 'Band5'): 1.83,
    ('C1', 'Band6'): 1.47,
    ('C1', 'Band7'): 0.98,
    ('C1', 'Band8'): 0.74,
    ('C1', 'Band9'): 0.52,
    ('C1', 'Band10'):0.39,

    ('C2', 'Band1'): ...,
    ('C2', 'Band3'): 2.30,
    ('C2', 'Band4'): 1.53,
    ('C2', 'Band5'): 1.24,
    ('C2', 'Band6'): 1.00,
    ('C2', 'Band7'): 0.67,
    ('C2', 'Band8'): 0.50,
    ('C2', 'Band9'): 0.35,
    ('C2', 'Band10'):0.26,

    ('C3', 'Band1'): ...,
    ('C3', 'Band3'): 1.42,
    ('C3', 'Band4'): 0.94,
    ('C3', 'Band5'): 0.77,
    ('C3', 'Band6'): 0.62,
    ('C3', 'Band7'): 0.41,
    ('C3', 'Band8'): 0.31,
    ('C3', 'Band9'): 0.22,
    ('C3', 'Band10'):0.16,

    ('C4', 'Band1'): ...,
    ('C4', 'Band3'): 0.92,
    ('C4', 'Band4'): 0.61,
    ('C4', 'Band5'): 0.50,
    ('C4', 'Band6'): 0.40,
    ('C4', 'Band7'): 0.27,
    ('C4', 'Band8'): 0.20,
    ('C4', 'Band9'): 0.14,
    ('C4', 'Band10'):0.11,

    ('C5', 'Band1'): ...,
    ('C5', 'Band3'): value,
    ('C5', 'Band4'): value,
    ('C5', 'Band5'): value,
    ('C5', 'Band6'): value,
    ('C5', 'Band7'): value,
    ('C5', 'Band8'): value,
    ('C5', 'Band9'): value,
    ('C5', 'Band10'): value,

    ('C6', 'Band1'): ...,
    ('C6', 'Band3'): value,
    ('C6', 'Band4'): value,
    ('C6', 'Band5'): value,
    ('C6', 'Band6'): value,
    ('C6', 'Band7'): value,
    ('C6', 'Band8'): value,
    ('C6', 'Band9'): value,
    ('C6', 'Band10'): value,

    ('C7', 'Band1'): ...,
    ('C7', 'Band3'): value,
    ('C7', 'Band4'): value,
    ('C7', 'Band5'): value,
    ('C7', 'Band6'): value,
    ('C7', 'Band7'): value,
    ('C7', 'Band8'): value,
    ('C7', 'Band9'): value,
    ('C7', 'Band10'): value,

    ('C8', 'Band1'): ...,
    ('C8', 'Band3'): value,
    ('C8', 'Band4'): value,
    ('C8', 'Band5'): value,
    ('C8', 'Band6'): value,
    ('C8', 'Band7'): value,
    ('C8', 'Band8'): value,
    ('C8', 'Band9'): value,
    ('C8', 'Band10'): value,

    ('C9', 'Band1'): ...,
    ('C9', 'Band3'): value,
    ('C9', 'Band4'): value,
    ('C9', 'Band5'): value,
    ('C9', 'Band6'): value,
    ('C9', 'Band7'): value,
    ('C9', 'Band8'): value,
    ('C9', 'Band9'): value,
    ('C9', 'Band10'): value,

    ('C10', 'Band1'): ...,
    ('C10', 'Band3'): value,
    ('C10', 'Band4'): value,
    ('C10', 'Band5'): value,
    ('C10', 'Band6'): value,
    ('C10', 'Band7'): value,
    ('C10', 'Band8'): value,
    ('C10', 'Band9'): value,
    ('C10', 'Band10'): value
}


    
    
    
    
    
    
    
    
    
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
                      
        if typ =='com07m':
            self.imsize = 256
            self.imcell = '1.50arcsec'
        elif typ =='com12m':
            self.imsize = 256
            self.imcell = '0.15arcsec'
            
        else:
            raise print("Wrong array element")
