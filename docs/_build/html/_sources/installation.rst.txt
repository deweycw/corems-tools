Installation
============

CoreMSTools can be installed directly from TestPyPi or from the source code. In both cases, we recommend installing in a virtual Python environment (unless running through a container).

.. note::

    Much of the functionality provided by these tools requires an installation of CoreMS **in the same Python environment**. If you are using a containerized version of CoreMS, you should install and run CoreMSTools in the CoreMS container. 

.. raw:: html

   <hr>

Install from TestPyPi    
---------------------

If you don't plan to make changes to source code and are simply interested in using coremstools, we recommend installing through PyPi. 

.. code-block::

    python -m pip install --index-url https://test.pypi.org/simple/ --no-deps coremstools

.. raw:: html

   <hr>

Install from source
-------------------------------

If you would like to make changes to the source code while also using the package, you can install coremstools by building it from source. Note that you can install an editable version of the package by including ``-e`` with ``pip install``.

.. code-block::

    git clone https://github.com/deweycw/corems-tools.git corems-tools
    cd corems-tools/src
    python -m pip install ./ # include '-e' flag during pip install for editable install