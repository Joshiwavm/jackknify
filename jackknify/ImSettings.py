import numpy as np

BINARY_VALUES = np.array([128, 256, 512, 1024, 2048])
RESOLUTIONS   = { #arcsec
    ('7m', 'Band1'): 31.5,
    ('7m', 'Band2'): ...,
    ('7m', 'Band3'): 12.5,
    ('7m', 'Band4'): 8.35,
    ('7m', 'Band5'): 6.77,
    ('7m', 'Band6'): 5.45,
    ('7m', 'Band7'): 3.63,
    ('7m', 'Band8'): 2.72,
    ('7m', 'Band9'): 1.93,
    ('7m', 'Band10'):1.44,

    ('C1', 'Band1'): 8.45,
    ('C1', 'Band2'): ...,
    ('C1', 'Band3'): 3.38,
    ('C1', 'Band4'): 2.25,
    ('C1', 'Band5'): 1.83,
    ('C1', 'Band6'): 1.47,
    ('C1', 'Band7'): 0.98,
    ('C1', 'Band8'): 0.74,
    ('C1', 'Band9'): 0.52,
    ('C1', 'Band10'):0.39,

    ('C2', 'Band1'): 5.75,
    ('C2', 'Band2'): ...,
    ('C2', 'Band3'): 2.30,
    ('C2', 'Band4'): 1.53,
    ('C2', 'Band5'): 1.24,
    ('C2', 'Band6'): 1.00,
    ('C2', 'Band7'): 0.67,
    ('C2', 'Band8'): 0.50,
    ('C2', 'Band9'): 0.35,
    ('C2', 'Band10'):0.26,

    ('C3', 'Band1'): 3.55,
    ('C3', 'Band2'): ...,
    ('C3', 'Band3'): 1.42,
    ('C3', 'Band4'): 0.94,
    ('C3', 'Band5'): 0.77,
    ('C3', 'Band6'): 0.62,
    ('C3', 'Band7'): 0.41,
    ('C3', 'Band8'): 0.31,
    ('C3', 'Band9'): 0.22,
    ('C3', 'Band10'):0.16,

    ('C4', 'Band1'): 2.30,
    ('C4', 'Band2'): ...,
    ('C4', 'Band3'): 0.92,
    ('C4', 'Band4'): 0.61,
    ('C4', 'Band5'): 0.50,
    ('C4', 'Band6'): 0.40,
    ('C4', 'Band7'): 0.27,
    ('C4', 'Band8'): 0.20,
    ('C4', 'Band9'): 0.14,
    ('C4', 'Band10'):0.11,

    ('C5', 'Band1'): 1.27,
    ('C5', 'Band2'): ...,
    ('C5', 'Band3'): 0.55,
    ('C5', 'Band4'): 0.36,
    ('C5', 'Band5'): 0.30,
    ('C5', 'Band6'): 0.24,
    ('C5', 'Band7'): 0.16,
    ('C5', 'Band8'): 0.12,
    ('C5', 'Band9'): 0.084,
    ('C5', 'Band10'):0.063,

    ('C6', 'Band1'): 0.78,
    ('C6', 'Band2'): ...,
    ('C6', 'Band3'): 0.31,
    ('C6', 'Band4'): 0.20,
    ('C6', 'Band5'): 0.17,
    ('C6', 'Band6'): 0.13,
    ('C6', 'Band7'): 0.089,
    ('C6', 'Band8'): 0.067,
    ('C6', 'Band9'): 0.047,
    ('C6', 'Band10'):0.035,

    ('C7', 'Band1'): 0.53,
    ('C7', 'Band2'): ...,
    ('C7', 'Band3'): 0.21,
    ('C7', 'Band4'): 0.14,
    ('C7', 'Band5'): 0.11,
    ('C7', 'Band6'): 0.092,
    ('C7', 'Band7'): 0.061,
    ('C7', 'Band8'): 0.046,
    ('C7', 'Band9'): 0.033,
    ('C7', 'Band10'):0.024,

    ('C8', 'Band1'): 0.240,
    ('C8', 'Band2'): ...,
    ('C8', 'Band3'): 0.096,
    ('C8', 'Band4'): 0.064,
    ('C8', 'Band5'): 0.052,
    ('C8', 'Band6'): 0.042,
    ('C8', 'Band7'): 0.028,
    ('C8', 'Band8'): 0.021,
    ('C8', 'Band9'): 0.015,
    ('C8', 'Band10'):0.011,

    ('C9', 'Band1'): 0.143,
    ('C9', 'Band2'): ...,
    ('C9', 'Band3'): 0.057,
    ('C9', 'Band4'): 0.038,
    ('C9', 'Band5'): 0.031,
    ('C9', 'Band6'): 0.025,
    ('C9', 'Band7'): 0.017,
    ('C9', 'Band8'): 0.012,
    ('C9', 'Band9'): 0.0088,
    ('C9', 'Band10'):...,

    ('C10', 'Band1'): 0.105,
    ('C10', 'Band2'): ...,
    ('C10', 'Band3'): 0.042,
    ('C10', 'Band4'): 0.028,
    ('C10', 'Band5'): 0.023,
    ('C10', 'Band6'): 0.018,
    ('C10', 'Band7'): 0.012,
    ('C10', 'Band8'): 0.0091,
    ('C10', 'Band9'): ...,
    ('C10', 'Band10'):...
}

FREQUENCIES = { #GHz
    ('Band1'):  40,
    ('Band2'): ...,
    ('Band3'): 100,
    ('Band4'): 150,
    ('Band5'): 185,
    ('Band6'): 230,
    ('Band7'): 345,
    ('Band8'): 460,
    ('Band9'): 650,
    ('Band10'):870 
}

class Imparams:
    def __init__(   self,
                    config,
                    band 
                ):
        
        self.config = config
        self.band   = band
        self.pb     = self._find_pb()
        self.freq   = FREQUENCIES[band]*1e9
        
        self.cellsize = self.find_cellsize()
        self.imsize   = self.find_imsize(self.find_FOV())//2

    def _find_pb(self):
        if self.config == '7m': return 7
        else: return 12

    def find_FOV(self):
        return 2.44 * np.rad2deg(299792458/self.freq/self.pb) * 3600  #arcsec

    def find_imsize(self, pb):
        ideal = pb/self.cellsize
        bit_value = BINARY_VALUES[np.argmin(abs(BINARY_VALUES - ideal))]
        if bit_value >= ideal:
            return bit_value
        elif (np.argmin(abs(BINARY_VALUES - ideal)) + 1) < 5: 
            return BINARY_VALUES[np.argmin(abs(BINARY_VALUES - ideal)) + 1]
        else: return BINARY_VALUES[-1]

    def _round_to_1(self, x):
        return round(x, -int(np.floor(np.log10(abs(x)))))

    def find_cellsize(self):
        res = RESOLUTIONS[self.config, self.band]/2./10.
        return self._round_to_1(res)
