Pycork Pre-Compiled Python Library for the Cork Library
=============================================================================

.. image:: https://github.com/drlukeparry/pycork/actions/workflows/pythonpublish.yml/badge.svg
    :target: https://github.com/drlukeparry/pycork/actions
.. image:: https://badge.fury.io/py/pycork.svg
    :target: https://badge.fury.io/py/pycork
.. image:: https://static.pepy.tech/personalized-badge/pycork?period=total&units=international_system&left_color=black&right_color=orange&left_text=Downloads
 :target: https://pepy.tech/project/pycork


Pycork is a Python library offering the functionality of the Cork boolean CSG library in a compiler friendly form suitable across all platforms. The library includes the dependencies for the Multi-Precision Integer and Rationals (MPIR) 3.0 built-in used by the Cork libary. The package aims to provide a simpler route for compiling the package for individuals and in addition, generating python bindings for use across other projects. Refactoring has been authored to tidy up the existing codebase so that it can be built across multiple platforms, in particular Windows, using the CMake build-system. At this stage, no further optimisations or improvements will be made specifically to the cork library, inclusive of its algorithms.

The python bindings are simple and offer access to the core functionality offered by the Cork library to perform boolean operations on watertight meshes. Additionally, it removes the awkward step of generatig .off files that are used in the command-line interface of the Cork library. The user may pass triangular meshes (vertices, tri-faces indices) as numpy arrays to each function.

For further information, see the latest `release notes <https://github.com/drlukeparry/pycork/blob/master/CHANGELOG.md>`_.

Installation
*************

Installation is currently supported on Windows. No special requiremnets are necessary for using pycork, except having the numpy library available. It is recommend to also install the `trimesh <https://github.com/mikedh/trimesh>`_ library to provide an interface to processing meshes as input for pycork.

.. code:: bash

    conda install -c numpy
    pip install trimesh

Installation of pycork can then be performed using pre-built python packages using the PyPi repository.

.. code:: bash

    pip install pycork

Alternatively, pycork may be compiled directly from source. Currently the prerequisites are the a compliant c++ build environment, include CMake build system. Currently the package has been tested on Windows 10, using VS2019.

.. code:: bash

    git clone https://github.com/drlukeparry/pycork.git && cd ./pycork
    git submodule update --init --recursive

    python setup.py install

Usage
******

The Cork CSG library, by design, has a simple interface for is functionality. Further detailed description of the function is therefore not necessary.

.. code:: python

    import numpy as np
    import trimesh

    import pycork

    # Note any manifold, watertight mesh can be used in conjuction with the Trimesh library
    meshA = trimesh.load_mesh('meshA.off')
    meshB = trimesh.load_mesh('meshB.off')

    # Extra list of vertices and triangular faces from the meshes
    vertsA = meshA.vertices
    trisA = meshA.faces

    vertsB = meshB.vertices
    trisB = meshB.faces

    pycork.isSolid(vertsA, trisA)
    pycork.isSolid(vertsB, trisB)

    #Perform the boolean opertions directly with Cork library
    vertsC, trisC = pycork.union(vertsA, trisA,
                                 vertsB, trisB)

    vertsD, trisD = pycork.intersection(vertsA, trisA,
                                        vertsB, trisB)


    meshC = trimesh.Trimesh(vertices=vertsC, faces=trisC, process=True)


