Jack-knife
==========

``jacked``, is a python-based package that jackknifes ALMA visibillites to create noise realizations of the observations.  Tutorials for installation and
usage can be found in the Notebooks folder. 

Methodology
==========

Jackknifing is a simplistic but effictive tool to quantifiy the underlying noise distribution of any type of data set. This tool specifically is implemented for interferometric data. ``jacked`` seperates half the visibillities randomly in two subsets, then multiplies one half with -1, so that when the data is binned or imaged any signal present in the data is averaged out. This creates observation-specific noise realization of the data, which can be used to for instance, sample the likelihood a false detection. 

The full methodology can be found 
`here <https://arxiv.org/abs/2210.03754>`_.

