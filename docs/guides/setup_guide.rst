Installation and setup
======================

The package's code was written in Python 3.5.x and is available on GitHub in its entirety, but relies on several dependencies. In case you are using a distribution like anaconda or similar, most should already be present on your machine.

Requirements
------------

After making sure your Python environment is up to date, for the data analysis part, install any the following packages not already present:

* NumPy
* SciPy
* matplotlib
* numba (version 0.35 or above)

Additionally, if you will be using this package to extract particle positions from your recordings, the data extraction functionality relies on the packages

* trackpy
* pims
* ctypes
* pyqt5
* pyqtgraph
* pandas

for full functionality, although it can be run without *pims* if you are using common formats.

Installation
------------

You can install the latest version of the package either directly from GitHub,

.. code-block:: bash

    pip install --upgrade https://github.com/fmf-pkp-2019/tweezer/tarball/master

for which Git is not required, or by first copying the repository,

.. code-block:: bash

    git clone https://github.com/fmf-pkp-2019/tweezer.git
    sudo python setup.py install

Verification
------------

Once ``Optical tweezer tools`` is installed, we suggest running the provided unit tests to make sure the code functions as intended on your machine. The file *run_tests.py* in the *tweezer/test* directory tests most functions; if all tests pass without failures or errors, the package has been (most likely) properly installed.
