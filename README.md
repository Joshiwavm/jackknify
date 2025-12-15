Jack-knife
==========
[![DOI](https://zenodo.org/badge/593247898.svg)](https://zenodo.org/doi/10.5281/zenodo.12516584)
[![Docs](https://img.shields.io/badge/docs-v1.0.0-2ea44f)](https://joshiwavm.github.io/jackknify/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


``jackknify``is a Python-based package that jackknifes ALMA visibilities to create noise realizations from the observations. 

Methodology
==========

Jackknifing is a simple but effective tool to characterize the underlying noise distribution of any type of data set. This tool specifically is implemented for interferometric data. ``jackknify`` splits half the visibilities randomly in two subsets, then multiplies one half with -1 so that when the data is binned, any signal present in the data is averaged out. This creates observation-specific noise realization of the data, which can be used to for instance, sample the likelihood a false detection. 

The full methodology can be found [here](https://ui.adsabs.harvard.edu/abs/2025A%26A...695A.204V/abstract) and is also used in [this work](https://arxiv.org/abs/2210.03754). 

Installation
============

``jackknify`` itself can be installed through

    pip install jackknify
    
or alternatively

    python -m pip install git+https://github.com/Joshiwavm/jackknify

or from the source

    git clone https://github.com/Joshiwavm/jackknify
    cd jackknify
    pip install -e .


## Dependancies

``jackknify`` uses ``casatask`` and ``casatools`` to interface with CASA measurements. ``casatask`` and ``casatools`` requires ``casadata`` to load. Sadly, this is a  ~350 MB sized file making the installment a bit slow. Further, when performing line searches, we make use of the package ``interferopy``, which is a Python-based package for common tasks used in the observational radio/mm interferometry data analysis.

## Trouble shooting casatask installation (if needed) 

If you want to run `jackknify` on a Mac with an Apple Silicon chip, run it in a Rosetta terminal. To open a Rosetta session in your terminal, run:
    
    /usr/bin/arch -x86_64 /bin/zsh --login

Further, `casadata` will download and store examples sets into the folder:  ~/.casa/data. However, it might not have permission from the local machine to do so. If such an error comes up. Just run:

    mkdir ~/.casa/data

To make the folder. That should solve most problems. 


Documentation
============

For your convenience, there are notebooks on how to run and use ``jackknify`` for line inference. You can find them in the docs/notebooks folder. Also, check out the documentation [here](https://joshiwavm.github.io/jackknify/).
