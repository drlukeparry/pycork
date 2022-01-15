// +-------------------------------------------------------------------------
// | cork.cpp
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
#include "cork.h"

#include <cork/mesh/mesh.h>

typedef Mesh<CorkVertex, CorkTriangle> CorkMesh;
typedef RawMesh<CorkVertex, CorkTriangle> RawCorkMesh;


void CorkVertex::merge(const CorkVertex &v0, const CorkVertex &v1) {
    double                              a0 = 0.5;
    if(v0.manifold && !v1.manifold)     a0 = 0.0;
    if(!v0.manifold && v1.manifold)     a0 = 1.0;

    // interpoalte between vertices
    pos = a0 * v0.pos + (1.0 - a0) * v1.pos;
}

void CorkVertex::interpolate(const CorkVertex &v0, const CorkVertex &v1) {
    double a0 = 0.5;
    double a1 = 0.5;
    pos = a0 * v0.pos + a1 * v1.pos;
}

void CorkVertex::isct(IsctVertEdgeTriInput<CorkVertex,CorkTriangle> input)
{
    Vec2d a_e = Vec2d(1,1)/2.0;
    Vec3d a_t = Vec3d(1,1,1)/3.0;

    a_e /= 2.0;
    a_t /= 2.0;
}
void CorkVertex::isct(IsctVertTriTriTriInput<CorkVertex,CorkTriangle> input)
{
    Vec3d       a[3];
    for(uint k=0; k<3; k++) {
        a[k]  = Vec3d(1,1,1)/3.0;
        a[k] /= 3.0;
    }
    for(uint i=0; i<3; i++) {
      for(uint j=0; j<3; j++) {
    }}
}
void CorkVertex::isctInterpolate(const CorkVertex &v0, const CorkVertex &v1) {
    double a0 = len(v1.pos - pos);
    double a1 = len(v0.pos - pos);

    if(a0 + a1 == 0.0)
        a0 = a1 = 0.5; // safety
    double sum = a0+a1;
    a0 /= sum;
    a1 /= sum;
}

void freeCorkTriMesh(CorkTriMesh *mesh)
{
    delete[] mesh->triangles;
    delete[] mesh->vertices;
    mesh->n_triangles = 0;
    mesh->n_vertices = 0;
}


//using RawCorkMesh = RawMesh<CorkVertex, CorkTriangle>;
//using CorkMesh = Mesh<CorkVertex, CorkTriangle>;

void eigenToCorkMesh(const Eigen::Matrix<double, Eigen::Dynamic, 3> &verts,
                     const Eigen::Matrix<uint64_t, Eigen::Dynamic, 3> &tris,
                     CorkMesh *mesh )
{
    RawCorkMesh raw;

    raw.vertices.resize(verts.rows());
    raw.triangles.resize(tris.rows());

    if(verts.rows() == 0 || tris.rows() == 0) {
        CORK_ERROR("empty mesh input to Cork routine.");
        *mesh = CorkMesh(raw);
        return;
    }

    /*
     * Due to Row-Major ordering in Numpy Arrays, it is more efficient to copy each row at a time
     */

    // X Coordinatees
    for(uint64_t i=0; i< tris.rows(); i++)
        raw.triangles[i].a = tris(i,0);

    for(uint64_t i=0; i< tris.rows(); i++)
        raw.triangles[i].b = tris(i,1);

    for(uint64_t i=0; i< tris.rows(); i++)
        raw.triangles[i].c = tris(i,2);

    uint64_t max_ref_idx = tris.maxCoeff();

    if(max_ref_idx > verts.rows()) {

        CORK_ERROR("mesh input to Cork routine has an out of range reference "
              "to a vertex.");

        raw.vertices.clear();
        raw.triangles.clear();

        *mesh = CorkMesh(raw);
        return;
    }

    for(uint64_t i=0; i<verts.rows(); i++)
        raw.vertices[i].pos.x = verts(i, 0);

    for(uint64_t i=0; i<verts.rows(); i++)
        raw.vertices[i].pos.y = verts(i, 1);

    for(uint64_t i=0; i<verts.rows(); i++)
        raw.vertices[i].pos.z = verts(i, 2);

    *mesh = CorkMesh(raw);

}

void corkMesh2Eigen(const CorkMesh &mesh,
                    Eigen::Matrix<double, Eigen::Dynamic, 3> &verts,
                    Eigen::Matrix<uint64_t, Eigen::Dynamic, 3> &tris)
{

    const RawCorkMesh raw = mesh.raw();

    verts.resize(raw.vertices.size(), Eigen::NoChange);
    tris.resize(raw.triangles.size(), Eigen::NoChange);

    // Note: Unfortunatly this is a relatively expensive copy due to change in Row/Column Ordering

    for(uint64_t i=0; i < tris.rows(); i++) {
        tris(i, 0) = raw.triangles[i].a;
        tris(i, 1) = raw.triangles[i].b;
        tris(i, 2) = raw.triangles[i].c;
    }

    for(uint64_t i=0; i< verts.rows(); i++) {
        verts(i,0) = raw.vertices[i].pos.x;
        verts(i,1) = raw.vertices[i].pos.y;
        verts(i,2) = raw.vertices[i].pos.z;
    }
}


void corkTriMesh2CorkMesh(CorkTriMesh in, CorkMesh *mesh_out) {
    RawCorkMesh raw;
    raw.vertices.resize(in.n_vertices);
    raw.triangles.resize(in.n_triangles);
    if(in.n_vertices == 0 || in.n_triangles == 0) {
        CORK_ERROR("empty mesh input to Cork routine.");
        *mesh_out = CorkMesh(raw);
        return;
    }
    
    uint max_ref_idx = 0;
    for(uint i=0; i<in.n_triangles; i++) {
        raw.triangles[i].a = in.triangles[3*i+0];
        raw.triangles[i].b = in.triangles[3*i+1];
        raw.triangles[i].c = in.triangles[3*i+2];
        max_ref_idx = std::max(
                        std::max(max_ref_idx,
                                 in.triangles[3*i+0]),
                        std::max(in.triangles[3*i+1],
                                 in.triangles[3*i+2])
                      );
    }
    if(max_ref_idx > in.n_vertices) {
        CORK_ERROR("mesh input to Cork routine has an out of range reference "
              "to a vertex.");
        raw.vertices.clear();
        raw.triangles.clear();
        *mesh_out = CorkMesh(raw);
        return;
    }
    
    for(uint i=0; i<in.n_vertices; i++) {
        raw.vertices[i].pos.x = in.vertices[3*i+0];
        raw.vertices[i].pos.y = in.vertices[3*i+1];
        raw.vertices[i].pos.z = in.vertices[3*i+2];
    }
    
    *mesh_out = CorkMesh(raw);
}

void corkMesh2CorkTriMesh(CorkMesh *mesh_in, CorkTriMesh *out)
{

    RawCorkMesh raw = mesh_in->raw();
    
    out->n_triangles = raw.triangles.size();
    out->n_vertices  = raw.vertices.size();
    
    out->triangles = new uint[(out->n_triangles) * 3];
    out->vertices  = new float[(out->n_vertices) * 3];
    
    for(uint i=0; i<out->n_triangles; i++) {
        (out->triangles)[3*i+0] = raw.triangles[i].a;
        (out->triangles)[3*i+1] = raw.triangles[i].b;
        (out->triangles)[3*i+2] = raw.triangles[i].c;
    }
    
    for(uint i=0; i<out->n_vertices; i++) {
        (out->vertices)[3*i+0] = raw.vertices[i].pos.x;
        (out->vertices)[3*i+1] = raw.vertices[i].pos.y;
        (out->vertices)[3*i+2] = raw.vertices[i].pos.z;
    }
}


bool isSolid(CorkTriMesh cmesh)
{
    CorkMesh mesh;
    corkTriMesh2CorkMesh(cmesh, &mesh);
    
    bool solid = true;
    
    if(mesh.isSelfIntersecting()) {
        CORK_ERROR("isSolid() was given a self-intersecting mesh");
        solid = false;
    }
    
    if(!mesh.isClosed()) {
        CORK_ERROR("isSolid() was given a non-closed mesh");
        solid = false;
    }
    
    return solid;
}

void computeUnion(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out) {

    CorkMesh cmIn0, cmIn1;

    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);
    
    cmIn0.boolUnion(cmIn1);
    
    corkMesh2CorkTriMesh(&cmIn0, out);
}

void computeDifference(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out) {

    CorkMesh cmIn0, cmIn1;
    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);
    
    cmIn0.boolDiff(cmIn1);
    
    corkMesh2CorkTriMesh(&cmIn0, out);
}

void computeIntersection(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out) {

    CorkMesh cmIn0, cmIn1;
    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);
    
    cmIn0.boolIsct(cmIn1);
    
    corkMesh2CorkTriMesh(&cmIn0, out);
}

void computeSymmetricDifference(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out){

    CorkMesh cmIn0, cmIn1;

    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);
    
    cmIn0.boolXor(cmIn1);
    
    corkMesh2CorkTriMesh(&cmIn0, out);
}

void resolveIntersections(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out) {

    CorkMesh cmIn0, cmIn1;

    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);
    
    cmIn0.disjointUnion(cmIn1);
    cmIn0.resolveIntersections();
    
    corkMesh2CorkTriMesh(&cmIn0, out);
}

