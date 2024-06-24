.. jackknify documentation master file, created by
   sphinx-quickstart on Tue Jun 18 20:11:29 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to jackknify's documentation!
=====================================

`jackknify` is a Python library used to create noise realizations of radio/mm interferometry data.

The package was developed to make realistic noise maps, which are usefull for feasibillity studies, as to better understand the false likelihood of detecting weak spectral lines. 

The package has been tested and used for ALMA data but the technique is broadly applicable to any interferometric data s.

To get started, check out the :doc:`install` page.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   tutorials

Please report any issues on the development page: https://github.com/Joshiwavm/jackknify


Acknowledgment
===============

If you make use of ``jackknify`` in your work, please cite it as: van Marrewijk (2024), using the following BibTeX entry::
   @software{joshiwa_van_marrewijk_2024_12516585,
     author       = {Joshiwa van Marrewijk and
                     Luca Di Mascolo},
     title        = {Joshiwavm/jackknify: Release Jackknify 0.1.1},
     month        = jun,
     year         = 2024,
     publisher    = {Zenodo},
     version      = {first\_release},
     doi          = {10.5281/zenodo.12516585},
     url          = {https://doi.org/10.5281/zenodo.12516585}
   }


