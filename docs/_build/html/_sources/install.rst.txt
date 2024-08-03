Installation
============

*jackknify* can be installed with ``pip``:

.. code-block:: bash

   pip install jackknify

Alternatively, install using GitHub

.. code-block:: bash

    python -m pip install git+https://github.com/Joshiwavm/jackknify

Or directly from the source

.. code-block:: bash

        git clone https://github.com/Joshiwavm/jackknify
        cd jackknify
        pip install -e .

*jackknify* uses ``casatasks`` and ``casatools``. This might cause a problem, so be wairy. 

Running *jackknify* on Apple Silicon chips
-------------

Due to the limited python verssion support of ``casatasks`` and ``casatools``, the only way for unning *jackknify* on an Apple Silicon computer is to work within a Rosetta terminal. To do so, open a standard terminal tab and type

.. code-block:: bash

   env /usr/bin/arch -x86_64 /bin/zsh --login

You can then proceed with the standard procecure for installing your favorite ``python`` or environment manager, and, finally, *jackknify*. To our knowledge, the most up-to-date version of ``python`` that can support the ``casatasks`` and ``casatools`` libraries is ``v3.8``.
