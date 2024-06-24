Jack-knife
==========

``jacked``is a Python-based package that jackknifes ALMA visibilites to create noise realizations from the observations. 

Methodology
==========

Jackknifing is a simplistic but effective tool to quantify the underlying noise distribution of any type of data set. This tool specifically is implemented for interferometric data. ``jacked`` separates half the visibilities randomly in two subsets, then multiplies one half with -1 so that when the data is binned or imaged, any signal present in the data is averaged out. This creates observation-specific noise realization of the data, which can be used to for instance, sample the likelihood a false detection. 

The full methodology can be found [here](https://arxiv.org/abs/2210.03754) and in an upcoming paper which is still in preparations. 

Installation
============

``jacked`` itself can be installed through

    pip install jacked
    
or alternatively

    python -m pip install git+https://github.com/Joshiwavm/jacked

or from the source

    git clone https://github.com/Joshiwavm/jacked
    cd jacked
    pip install -e .


## Dependancies

``jacked`` uses ``casatask`` and ``casatools`` to interface with CASA measurements. ``casatask`` and ``casatools`` requires ``casadata`` to even laod. Sadly, this is a \approx 350 MB sized file. Further, when performing line searches, we make use of the package ``interferopy``, which is a Python-based package for common tasks used in the observational radio/mm interferometry data analysis.

## Mac 

If you want to run `jacked` on a Mac that has an Apple Silicon chip, run it in a Rosetta terminal. To do this, just run:
    
    /usr/bin/arch -x86_64 /bin/zsh --login

Documentation
============

For your convenience, there are notebooks on how to run ``jacked``. You can find them in the docs/notebooks folder. Also, check out the documentation [here](...).

Acknowledgment
============

If you make use of jacked in your work, please cite it as **, using the following BibTeX entry:

```

```
