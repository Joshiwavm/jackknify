from src.Models import *

COMPONENTS = {
    "pointSource": {
        "id": 1,
        "variables": [ 
            'RA','Dec','Amplitude','Offset'
        ],
        "make_image"  :   False,
        "SZ"          :   False, 
        "spectrum"    :   False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','Amplitude (Jy)','Offset'
        ]
    },
   "gaussSource": {
        "id": 2,
        "variables": [ 
            'RA','Dec','Amplitude','Major','e','Angle','Offset'
        ],
        "make_image"  :   False,
        "SZ"          :   False, 
        "spectrum"    :   False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','Amplitude (Jy)','Major (deg)','e','Angle (deg)','Offset'
        ]
    },
    "gaussSurface": {
        "id": 3,
        "variables": [ 
            'RA', 'Dec', 'Amplitude', 'Major', 'e', 'Angle', 'Offset', 'Temperature'
        ],
        "make_image"  :    True,
        "SZ"          :   False, 
        "spectrum"    :   False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','Amplitude (Jy)','Major (deg)','e','Angle (deg)','Offset','Temperature (keV)'
        ]
    },
    "betaPressure": {
        "id": 4,
        "function": betaProfile,
        "variables": [ 
           'RA', 'Dec', 'Amplitude', 'Major', 'e', 'Angle', 'Offset', 'Temperature', 'Beta', 'z'
        ],
        "make_image"  :   True,
        "SZ"          :   True, 
        "spectrum"    :  False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','Amplitude (keV/cm3)','Major (deg)','e','Angle (deg)','Offset','Temperature (keV)','Beta','z'
        ]
    },
    "gnfwPressure": {
        "id": 5,
        "function": gnfwProfile,
        "variables": [ 
            'RA', 'Dec', 'Amplitude', 'Major', 'e', 'Angle', 'Offset', 'Temperature', 'Alpha', 'Beta', 'Gamma', 'z'
        ],
        "make_image"  :   True,
        "SZ"          :   True, 
        "spectrum"    :  False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','Amplitude (keV/cm3)','Major (deg)','e','Angle (deg)','Offset','Temperature (keV)','Beta', 'Gamma', 'z'
        ]
    },
    "A10Pressure": {
        "id": 6,
        "function": a10Profile,
        "variables": [ 
            'RA', 'Dec', 'log10', 'c500', 'e', 'Angle', 'Offset', 'Temperature', 'Alpha', 'Beta', 'Gamma', 'P0', 'Alpha_p', 'z', 'bias'
        ],
        "make_image"  :   True,
        "SZ"          :   True, 
        "spectrum"    :  False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','log10(M500/Msun)','c500','e','Angle (deg)','Offset','Temperature (keV)','Alpha','Beta','Gamma','P0','Alpha_p','z','bias'
        ]
    },
    "A10PressureLS": {
        "id": 6,
        "function": a10Profile,
        "variables": [ 
            'RA', 'Dec', 'log10', 'c500', 'e', 'Angle', 'Offset', 'Temperature', 'Alpha', 'Beta', 'Gamma', 'P0', 'Alpha_p', 'z', 'bias', 'depth'
        ],
        "make_image"  :   True,
        "SZ"          :   True, 
        "spectrum"    :  False, 
        "var_to_print":[
            'RA (deg)','Dec (deg)','log10(M500/Msun)','c500','e','Angle (deg)','Offset','Temperature (keV)','Alpha','Beta','Gamma','P0','Alpha_p','z','bias', 'depth'
        ]
    },
    "powerLaw": {
        "id": 7,
        "variables": [ 
            "SpecIndex"
        ],
        "make_image"  :  False,
        "SZ"          :  False, 
        "spectrum"    :   True, 
        "var_to_print":[
            'SpecIndex'
        ]
    },
    "powerLawMod": {
        "id": 8,
        "variables": [ 
            "SpecIndex", "SpecCurv"
        ],
        "make_image"  :  False,
        "SZ"          :  False, 
        "spectrum"    :   True, 
        "var_to_print":[
            'SpecIndex','SpecCurv'
        ]
    },
    "powerDust": {
        "id": 9,
        "variables": [ 
            "SpecIndex", "SpecCurv", "Temp", "beta", "z", "kappa0", "nu0"
        ],
        "make_image"  :  False,
        "SZ"          :  False, 
        "spectrum"    :   True, 
        "var_to_print":[
            'SpecIndex','log(Mass/M_sun)','Temperature (keV)','Beta','z','kappa0','nu0'
        ]
    },
    "tSZ": {
        'id':  10, 
        "variables": [ 
        ],
        "make_image"  :  False,
        "SZ"          :  False, 
        "spectrum"    :   True, 
        "var_to_print":[

        ]
    },
    "Scaling": {
        'id':  11, 
        "variables": [ 
            'amp'
        ],
        "make_image"  :  False,
        "SZ"          :  False, 
        "spectrum"    :   True, 
        "var_to_print":[
            'alpha'
        ]
    }
}
