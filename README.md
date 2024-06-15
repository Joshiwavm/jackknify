Jack-knife
==========

``jacked``, is a python-based package that jackknifes ALMA visibillites to create noise realizations of the observations. 

Methodology
==========

Jackknifing is a simplistic but effictive tool to quantifiy the underlying noise distribution of any type of data set. This tool specifically is implemented for interferometric data. ``jacked`` seperates half the visibillities randomly in two subsets, then multiplies one half with -1, so that when the data is binned or imaged any signal present in the data is averaged out. This creates observation-specific noise realization of the data, which can be used to for instance, sample the likelihood a false detection. 

The full methodology can be found [here](https://arxiv.org/abs/2210.03754)

Installation
============

``jacked`` uses ``casatask`` and ``casatools`` to interface with CASA measurements. ``casatask`` and ``casatools`` requires `python=3.6`. This requirement limits the compatabillity of jacked sadly enough with other packages such as [interferopy](https://interferopy.readthedocs.io/en/latest/index.html). Therefore, we recommand (even though annoying) to install ``jacked`` in a seperate environment from the one used to do the data analysis on. 

``jacked`` itselves can be installed through::

    pip install jacked

or alternatively through::

python -m pip install git+https://github.com/Joshiwavm/jacked

or 

    git clone https://github.com/Joshiwavm/jacked
    cd jacked
    pip install -e .

Documentation
============

tutorials. 
