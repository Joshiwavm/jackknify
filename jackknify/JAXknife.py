import jax; jax.config.update('jax_enable_x64',True)
import jax.numpy as jp

import jax_finufft

import numpy as np

def run(xw,yw,dw,w,csize,rng):
    flips = np.ones(xw.shape[0])
    flips[:xw.shape[0]//2] = -1.00
    rng.shuffle(flips)
    
    return image(xw,yw,dw*flips,w,csize)

def image(x,y,c,w,csize):
    x = jp.append(x,-x)
    y = jp.append(y,-y)
    w = jp.append(w, w)
    c = jp.append(c,c.conj())
    return jax_finufft.nufft1((csize,csize),c*w/np.sum(w),x,y).real