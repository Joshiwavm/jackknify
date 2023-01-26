from src import jacked

#which observations do you want to corrupt --> could do this in CLE
VIS        = './data/Workable/HD1_Band6_tbin30s_cwidth30kms.ms' 
FIELDS     = ['0']
SPWS       = [['25','27','29','31','57','58','59','60','85','86','87','88','113','114','115','116']]
BAND       = ['6']
ARRAY      = ['com12m']

#how many new samples do you need
N          = 100

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