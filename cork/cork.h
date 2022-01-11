// +-------------------------------------------------------------------------
// | cork.h
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2013
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the Cork library.
// |
// |    Cork is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    Cork is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with Cork.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#ifndef CORK_CORK_H_HEADER_HAS_BEEN_INCLUDED
#define CORK_CORK_H_HEADER_HAS_BEEN_INCLUDED

#include "CORK_Export.h"

#include <cork/mesh/mesh.h>
#include <cork/rawmesh/rawMesh.h>

#include <Eigen/Dense>

#ifndef uint
typedef unsigned int uint;
#endif

// if a mesh is taken as input, the client must manage the memory
// if a mesh is given as output, please use the provided
// function to free the allocated memory.
struct CorkTriMesh
{
    uint    n_triangles;
    uint    n_vertices;
    uint    *triangles;
    float   *vertices;
};

struct CorkTriangle;

struct CORK_EXPORT CorkVertex : public MinimalVertexData,
                                public RemeshVertexData,
                                public IsctVertexData,
                                public BoolVertexData

{
    void merge(const CorkVertex &v0, const CorkVertex &v1);
    void interpolate(const CorkVertex &v0, const CorkVertex &v1);
    void isct(IsctVertEdgeTriInput<CorkVertex,CorkTriangle> input);
    void isct(IsctVertTriTriTriInput<CorkVertex,CorkTriangle> input);
    void isctInterpolate(const CorkVertex &v0, const CorkVertex &v1);
};


struct CORK_EXPORT CorkTriangle :
        public MinimalTriangleData,
        public RemeshTriangleData,
        public IsctTriangleData,
        public BoolTriangleData
{
    void merge(const CorkTriangle &, const CorkTriangle &) {}
    static void split(CorkTriangle &, CorkTriangle &,
                      const CorkTriangle &) {}
    void move(const CorkTriangle &) {}
    void subdivide(SubdivideTriInput<CorkVertex,CorkTriangle> input)
    {
        bool_alg_data = input.pt->bool_alg_data;
    }
};



typedef Mesh<CorkVertex, CorkTriangle> CorkMesh;

CORK_EXPORT void eigenToCorkMesh(const Eigen::Matrix<double, Eigen::Dynamic, 3> &verts,
                                 const Eigen::Matrix<uint64_t, Eigen::Dynamic, 3> &tris,
                                 CorkMesh *triMesh);

CORK_EXPORT void corkMesh2Eigen(const CorkMesh &mesh,
                                Eigen::Matrix<double, Eigen::Dynamic, 3> &verts,
                                Eigen::Matrix<uint64_t, Eigen::Dynamic, 3> &tris);

CORK_EXPORT void freeCorkTriMesh(CorkTriMesh *mesh);

// the inputs to Boolean operations must be "solid":
//  -   closed (aka. watertight; see comment at bottom)
//  -   non-self-intersecting
// additionally, inputs should use a counter-clockwise convention
// for triangle facing.  If the triangles are presented in clockwise
// orientation, the object is interpreted as its unbounded complement

// This function will test whether or not a mesh is solid
CORK_EXPORT bool isSolid(CorkTriMesh mesh);

// Boolean operations follow
// result = A U B
CORK_EXPORT void computeUnion(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

// result = A - B
CORK_EXPORT void computeDifference(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

// result = A ^ B
CORK_EXPORT void computeIntersection(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

// result = A XOR B
CORK_EXPORT void computeSymmetricDifference(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

// Not a Boolean operation, but related:
//  No portion of either surface is deleted.  However, the
//  curve of intersection between the two surfaces is made explicit,
//  such that the two surfaces are now connected.
CORK_EXPORT void resolveIntersections(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

#endif
