# -*- coding: utf-8 -*-

import unittest
import platform
import tempfile

import numpy as np
import trimesh

import context

import pycork

class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    tolerance = 1e-6

    def test_isSolid(self):

        
        # Note any manifold, watertight mesh can be used in conjuction with the Trimesh library
        meshA = trimesh.load_mesh('../models/meshA.off')
        meshB = trimesh.load_mesh('../models/meshB.off')

        # Extra list of vertices and triangular faces from the meshes
        vertsA = meshA.vertices
        trisA = meshA.faces

        vertsB = meshB.vertices
        trisB = meshB.faces

        pycork.isSolid(vertsA, trisA)
        pycork.isSolid(vertsB, trisB)

    
     def test_union(self):

        
        # Note any manifold, watertight mesh can be used in conjuction with the Trimesh library
        meshA = trimesh.load_mesh('../models/meshA.off')
        meshB = trimesh.load_mesh('../models/meshB.off')

        # Extra list of vertices and triangular faces from the meshes
        vertsA = meshA.vertices
        trisA = meshA.faces

        vertsB = meshB.vertices
        trisB = meshB.faces

        #Perform the boolean opertions directly with Cork library
        vertsC, trisC = pycork.union(vertsA, trisA,
                                     vertsB, trisB)

        meshC = trimesh.Trimesh(vertices=vertsC, faces=trisC, process=True)
        

    def test_difference(self):

        
        # Note any manifold, watertight mesh can be used in conjuction with the Trimesh library
        meshA = trimesh.load_mesh('../models/meshA.off')
        meshB = trimesh.load_mesh('../models/meshB.off')

        # Extra list of vertices and triangular faces from the meshes
        vertsA = meshA.vertices
        trisA = meshA.faces

        vertsB = meshB.vertices
        trisB = meshB.faces


        vertsD, trisD = pycork.difference(vertsA, trisA,
                                        vertsB, trisB)


        meshC = trimesh.Trimesh(vertices=vertsC, faces=trisC, process=True)
    
    
    def test_intersection(self):

        myPart = pyslm.Part('partName')

        # Create test geometry cube
        dims = np.array((2,3,4))
        cube = trimesh.creation.box(dims)

        myPart.setGeometry(cube)

        # Test the geometry loading and the volume
        assert abs(myPart.volume - np.prod(dims)) < self.tolerance

        # Test the rotation and transformation
        myPart.rotation = [0,90,0]
        assert abs(myPart.extents[2] - dims[0]) < self.tolerance
        

if __name__ == '__main__':
    unittest.main()
