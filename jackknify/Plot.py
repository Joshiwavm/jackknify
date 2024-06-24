import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
from tqdm import tqdm
from astropy.io import fits
from astropy.wcs import WCS
from astropy.constants import c
from astropy import units as u

from IPython.display import clear_output

def circle_mask(im, xc, yc, rcirc):
        ny, nx = im.shape
        y,x = np.mgrid[0:nx,0:ny]
        r = np.sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc))
        return ( (r < rcirc))

def ellipsoid_mask(im, xc, yc, a, b, Theta):
    ny, nx = im.shape
    y,x = np.mgrid[0:nx,0:ny]
    e = (((x-xc) * np.cos(Theta) + (y - yc)*np.sin(Theta))**2)/a**2 +  (((x-xc) * np.sin(Theta) - (y - yc)*np.cos(Theta))**2)/b**2
    return (e<=1)

class SLP():
    
    def __init__(self, 
                 fname: str,
                 size: float, 
                 idx:int = 0,
                 amount: int = 200, 
                 visualize: bool = False,
                 center_coord: tuple = (None, None)
                 ):

        # Load in Cube
        # ------------
        self.fname = fname
        self.hdu = fits.open(self.fname)
        self.cube   = self.hdu[idx].data[0] if self.hdu[idx].data.ndim > 3 else self.hdu[idx].data[None]
        self.header = self.hdu[idx].header
        self.pixel_size = self.header['CDELT1']
        self.center_coord = center_coord
        
        # Get Beam
        try:
            BMAJ = np.average(self.hdu[1].data['BMAJ'])/3600  # arcseconds
            BMIN = np.average(self.hdu[1].data['BMIN'])/3600  # arcseconds
            BPA = np.average(self.hdu[1].data['BPA'])  # Position angle in degrees
        except:
            BMAJ = self.header['BMAJ']   # arcseconds
            BMIN = self.header['BMIN']   # arcseconds
            BPA = self.header['BPA']  # Position angle in degrees

        # convert
        BeamArea     = np.pi *BMIN*BMAJ/self.pixel_size**2/(4*np.log(2))   # pixels per beam 
        self.cube    = self.cube/BeamArea*1e3                               # mJy/pixel

        # make mask
        wcs = WCS(self.header, naxis = 2)
        if self.center_coord[0] == None:
            center_x, center_y = (self.cube.shape[-2]//2, self.cube.shape[-1]//2)
        else:
            center_x, center_y = wcs.world_to_pixel_values(self.center_coord[0], self.center_coord[1])
        
        self.mask = ellipsoid_mask(self.cube[0], 
                                   center_x, 
                                   center_y, 
                                   BMIN/self.pixel_size * size, 
                                   BMAJ/self.pixel_size * size, 
                                   np.deg2rad(BPA)
                                   )
        

        self.visualize = visualize
        self.amount = amount
        
        xx, yy         = np.meshgrid(np.arange(0,self.cube.shape[-2], 1), np.arange(0, self.cube.shape[-1],1))
        self.CoM_mask  =  (np.mean(xx[self.mask]), np.mean(yy[self.mask]))
        self.rr        = ((xx-self.CoM_mask[0])**2 + (yy-self.CoM_mask[1])**2)**0.5
        self.mask_size = int(np.sum(self.mask)**0.5+0.1*self.cube.shape[-2]) # --> quick fix

        self.run()

    def _position(self):

        r =  np.sqrt(np.random.uniform(0.1,1, size = self.amount))
        theta = np.random.uniform(0,1, size = self.amount) * 2 * np.pi

        x = self.CoM_mask[0] + self.cube.shape[1]/3 *r* np.cos(theta) #hard coded
        y = self.CoM_mask[1] + self.cube.shape[1]/3 *r* np.sin(theta) #hard coded
        
        pos = np.array([x,y], dtype = np.int).T
        return np.vstack(([int(self.CoM_mask[0]), int(self.CoM_mask[1])], pos))

    def _make_mask(self, pos):
        new_masks = np.zeros((len(pos[1:]),self.im.shape[-2], self.im.shape[-1]))
        for idx, xy in enumerate(pos[1:]):
            new_masks[idx, 
                      xy[0]-self.mask_size//2:xy[0]+self.mask_size//2, 
                      xy[1]-self.mask_size//2:xy[1]+self.mask_size//2] = self.mask[pos[0][0]-self.mask_size//2:pos[0][0]+self.mask_size//2, 
                                                                                   pos[0][1]-self.mask_size//2:pos[0][1]+self.mask_size//2]
            
        return new_masks.astype(bool)  
    
    def _estimates(self, new_masks):
        estimates = []
        for m in new_masks:
            estimates.append(np.sum(self.im[m]))
        return np.nanstd(estimates), np.nanmean(estimates)
    
    def _visualize(self, masks):
        for idx, m in enumerate(masks):
            clear_output(True)
            plt.imshow(m, origin='lower')
            plt.show()
            if idx>30: break
    
    def run(self):
        bootstrap_std = []
        bootstrap_means = []
        for i in tqdm(range(len(self.cube))):
            self.im = self.cube[i]
            pos     = self._position()
            masks   = self._make_mask(pos)            

            std, mean = self._estimates(masks)
            bootstrap_std.append(std)
            bootstrap_means.append(mean)
            
        if self.visualize: self._visualize(masks)

        self.stds  = np.array(bootstrap_std)
        self.means = np.array(bootstrap_means)
    
    def plot(self,savedir,savenp):

        Dfreq    = self.header['CDELT3']
        restfreq = self.header['CRVAL3']

        xaxis      = (np.arange(len(self.cube))*Dfreq + restfreq)/1e9    
        slp        = np.nansum(self.cube[:, self.mask], axis = 1)

        if savenp:
            np.save(file = savedir + self.fname.split('/')[-1].replace('.fits', '_slp.npy'), arr =[slp, self.stds])

        # Visualize SLP
        # -----------------
        fig, ax = plt.subplots(1,1, sharex=True, figsize=(7,4))#, gridspec_kw={'height_ratios': [2, 1]})

        ax.step(xaxis, slp, label = 'Central Beam', c = 'C0')
        ax.step(xaxis, self.stds, label = r'1$\sigma$', c = 'C1', ls =':')
        ax.step(xaxis, -1*self.stds, c = 'C1', ls =':')
        ax.axhline(0, c='gray', ls='--')
        ax.axvline(restfreq/1e9, c = 'k')
        ax.set_xlim(xmin = xaxis[0], xmax = xaxis[-1])
        ax.set_ylabel('SNR')
        ax.legend(loc = 1, frameon=False)

        plt.tight_layout()
        plt.savefig(savedir + self.fname.split('/')[-1].replace('.fits', '_slp.pdf'), dpi = 300)
        plt.show()

class IMAGE():
    def __init__(self,
                 fname:str,
                 moment:str = 'continuum',
                 idx:int = 0,
                 channels:int = None,
                 center_coord: tuple = (None, None), 
                 box_size_arcsec: float = None
                ):
        
        self.fname    = fname
        self.hdu      = fits.open(fname)
        self.header   = self.hdu[idx].header
        self.cube     = self.hdu[idx].data[0] if self.hdu[idx].data.ndim > 3 else self.hdu[idx].data[None]
        self.channels = channels
        self.moment   = moment

        self.center_coord = center_coord
        self.box_size_arcsec = box_size_arcsec
        self._make_image()

    def _make_image(self,):
        data = self.cube.copy()
        if self.moment == 'continuum':
            self.img = np.nanmean(data, axis = 0)*1e3
        elif self.moment == 'moment-0':
            dv = self.header['CDELT3']/self.header['CRVAL3'] * c.to(u.km/u.s).value
            self.img = np.nansum(data[self.channels[0]:self.channels[1]]*dv, axis=0)*1e3
        else:
            raise ValueError(f"Moment does not match available inputs: continuum, or moment-0")

    def plot(self,
              savedir:str = './',
              ):
        
        # Get the WCS information from the header
        wcs = WCS(self.header, naxis = 2)
        
        try:
            beam_major = np.average(self.hdu[1].data['BMAJ'])  # arcseconds
            beam_minor = np.average(self.hdu[1].data['BMIN'])  # arcseconds
            beam_pa = np.average(self.hdu[1].data['BPA'])  # Position angle in degrees
        except:
            beam_major = self.header['BMAJ'] * 3600  # arcseconds
            beam_minor = self.header['BMIN'] * 3600  # arcseconds
            beam_pa = self.header['BPA']  # Position angle in degrees
            
        if self.center_coord[0] == None:
            center_x, center_y = (self.img.shape[-2]//2, self.img.shape[-1]//2)
        else:
            center_x, center_y = wcs.world_to_pixel_values(self.center_coord[0], self.center_coord[1])
        
        # Convert box size from arcseconds to pixels
        pix_scale_x = np.abs(wcs.wcs.cdelt[0]) * 3600  # arcsec per pixel in x direction
        pix_scale_y = np.abs(wcs.wcs.cdelt[1]) * 3600  # arcsec per pixel in y direction
        
        if self.box_size_arcsec == None:
            box_size_x_pix = self.img.shape[0]
            box_size_y_pix = self.img.shape[1]
        else:
            box_size_x_pix = self.box_size_arcsec / pix_scale_x
            box_size_y_pix = self.box_size_arcsec / pix_scale_y

        
        # Determine the bounding box in pixel coordinates
        x_min = int(center_x - box_size_x_pix / 2)
        x_max = int(center_x + box_size_x_pix / 2)
        y_min = int(center_y - box_size_y_pix / 2)
        y_max = int(center_y + box_size_y_pix / 2)
        
        # Estimate the standard deviation outside the bounding box
        mask = np.ones(self.img.shape, dtype=bool)
        mask[y_min:y_max, x_min:x_max] = False
        
        if np.sum(mask) == 0:
            std_dev = np.nanstd(self.img)
        else:
            std_dev = np.nanstd(self.img[mask])

        # Crop the image to the bounding box
        cropped_img = self.img[y_min:y_max, x_min:x_max]

        # Adjust WCS for the cropped image
        cropped_wcs = wcs[y_min:y_max, x_min:x_max]

        # Plot the cropped image with WCS projection
        fig, ax = plt.subplots(subplot_kw={'projection': cropped_wcs}, figsize=(5,4))
        im = ax.imshow(cropped_img, origin='lower', cmap='RdBu_r', vmin = np.nanmin(self.img), vmax = -1*np.nanmin(self.img))
        ax.set_xlabel('RA (J2000)', fontsize = 12)
        ax.set_ylabel('DEC (J2000)', fontsize = 12)
        
        # Plot contours based on the estimated standard deviation
        levels = [-4*std_dev, -3*std_dev, -2*std_dev, -1*std_dev, 1*std_dev, 2*std_dev, 3*std_dev, 4*std_dev]  # Example contour levels
        ax.contour(cropped_img, levels=levels, colors='k', alpha=0.7)
        
        # Add a colorbar
        cbar = plt.colorbar(im, ax=ax, orientation='vertical')
        if self.moment == 'continuum':
            cbar.set_label(r'$S_{\nu}$ [mJy/Beam]', fontsize = 12)
        elif self.moment == 'moment-0':
            cbar.set_label(r'$S_{\nu}$ [mJy/Beam km/s]', fontsize = 12)

        # Plot the beam size in the lower right corner
        beam_major_pix = beam_major / pix_scale_x
        beam_minor_pix = beam_minor / pix_scale_y    
        beam = Ellipse((0.1*box_size_x_pix, 0.1*box_size_y_pix), width=beam_minor_pix, height=beam_major_pix,
                    angle=beam_pa, color='gray', fill='//')
        ax.add_patch(beam)
        
        plt.tight_layout()
        plt.savefig(savedir + self.fname.split('/')[-1].replace('.fits', '_'+self.moment+'.pdf'), dpi = 300)
        plt.show()
