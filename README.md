Jack-knife
==========
[![DOI](https://zenodo.org/badge/593247898.svg)](https://zenodo.org/doi/10.5281/zenodo.12516584)

``jackknify``is a Python-based package that jackknifes ALMA visibilities to create noise realizations from the observations. 
``jackknify``is a Python-based package that jackknifes ALMA visibilities to create noise realizations from the observations. 

Methodology
==========

Jackknifing is a simplistic but effective tool to quantify the underlying noise distribution of any type of data set. This tool specifically is implemented for interferometric data. ``jackknify`` separates half the visibilities randomly in two subsets, then multiplies one half with -1 so that when the data is binned or imaged, any signal present in the data is averaged out. This creates observation-specific noise realization of the data, which can be used to for instance, sample the likelihood a false detection. 

The full methodology can be found [here](https://arxiv.org/abs/2210.03754) and in an upcoming paper, which is still in preparation. 

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

``jackknify`` uses ``casatask`` and ``casatools`` to interface with CASA measurements. ``casatask`` and ``casatools`` requires ``casadata`` to load. Sadly, this is a  ~350 MB sized file. Further, when performing line searches, we make use of the package ``interferopy``, which is a Python-based package for common tasks used in the observational radio/mm interferometry data analysis.

## Mac 

If you want to run `jackknify` on a Mac with an Apple Silicon chip, run it in a Rosetta terminal. To open a Rosetta session in your terminal, run:
    
    /usr/bin/arch -x86_64 /bin/zsh --login


Documentation
============

For your convenience, there are notebooks on how to run ``jackknify``. You can find them in the docs/notebooks folder. Also, check out the documentation [here](https://joshiwavm.github.io/jackknify/).
