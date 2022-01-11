Pycork Pre-Compiled Python Library for the Cork Library
=============================================================================

.. image:: https://github.com/drlukeparry/pycork/actions/workflows/pythonpublish.yml/badge.svg
    :target: https://github.com/drlukeparry/pycork/actions
.. image:: https://badge.fury.io/py/pycork.svg
    :target: https://badge.fury.io/py/pycork
.. image:: https://static.pepy.tech/personalized-badge/pycork?period=total&units=international_system&left_color=black&right_color=orange&left_text=Downloads
 :target: https://pepy.tech/project/pycork


Pycork is a Python library offering the functionality of the Cork boolean CSG library in a compiler friendly form suitable across all platforms. The library includes the dependencies for the Multi-Precision Integer and Rationals (MPIR) used by the Cork libary to provide a simpler route for compiling the package individually and also python bindings. Refactoring has been done to tidy up the existing code base that can be built using CMake and on Windows.

The python bindings are simple and offer access to  the core functionality offered by the Cork library. Additionally, it removes the awkard requirement to generate .off files used in the command-line option for Cork library. This is simply done by passing the triangular meshes (vertices, tri faces) as numpy arrays as arugments for each function.

For further information, see the latest `release notes <https://github.com/drlukeparry/pycork/blob/master/CHANGELOG.md>`_.

Installation
*************

Installation is currently supported on Windows. No special requiremets are necessary for using pycork, except having the numpy library available. It is recommends to additionally the `trimesh <https://github.com/mikedh/trimesh>`_ library to provide an interface to processing meshes as input for pycork.

.. code:: bash

    conda install -c numpy
    pip install trimesh

Installation of pycork can then be performed using pre-built python packages using the PyPi repository.

.. code:: bash

    pip install pycork

Alternatively, pycork may be compiled directly from source. Currently the prerequisites are the a compliant c++ build environment, include CMake build system. Currently tested on Windows 10, using VS19.0

.. code:: bash

    git clone https://github.com/drlukeparry/pyslm.git && cd ./pyslm
    git submodule update --init --recursive

    python setup.py install

Usage
******

The Cork CSG library is simple in nature and there are few functions that require any extra explicit description.

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


