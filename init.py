from src import jacked

#which observations do you want to corrupt
VIS        = './PATH/TO/MSfile' 
FIELDS     = ['field1', 'field2']
SPWS       = [['spw1','spw2'], ['spw1', 'spw2']]
BAND       = ['BandX']
ARRAY      = ['Your configuration'] # could be any of the following: 7m, C1, C2, CX...

#how many new samples do you need
N          = ...

#initilaize
tool = jacked.Jack(vis     = VIS, 
                   spws    = SPWS, 
                   fields  = FIELDS, 
                   band    = BAND,
                   array   = ARRAY,
                   samples = N
                   )

#run
tool.run()