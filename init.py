from Jack_it import *

vis        = './data/Workable/HD1_Band6_tbin30s_cwidth30kms.ms' 
spws       = [['25','27','29','31','57','58','59','60','85','86','87','88','113','114','115','116']]
fields     = ['0']

for i in range(0,20):
    typ        = 'com12m'
    model_name = ''

    tool = Test_SZonly(vis, typ, model_name, seed = 195758392+i, spws = spws, fields = fields, test = False)
    print('..loading')
    tool._loader()
    print("..Jack_knifing")
    tool.Jack_it()
    print('..saving')
    tool._saver()
    Img = tool._image()
