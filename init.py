from Jack_it.jacked import *

vis        = './data/Workable/HD1_Band6_tbin30s_cwidth30kms.ms' 
fields     = ['0']
spws       = [['25','27','29','31','57','58','59','60','85','86','87','88','113','114','115','116']]
typ        = 'com12m'
N          = 20

tool = Jack(vis, typ, seed = 195758392, spws = spws, fields = fields, test = False)
tool.run()