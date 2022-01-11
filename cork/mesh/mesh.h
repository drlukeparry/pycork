// +-------------------------------------------------------------------------
// | mesh.h
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
#ifndef CORK_MESH_H_HEADER_HAS_BEEN_INCLUDED
#define CORK_MESH_H_HEADER_HAS_BEEN_INCLUDED

// SIMPLE USAGE:
//  In order to get standard template inclusion behavior/usage
//  just include this file.  Then the entire template code
//  will be included in the compilation unit.
// ADVANCED USAGE:
//  Only include "mesh.decl.h" where-ever you would normally include a
//  header file.  This will avoid including the implementation code in
//  the current compilation unit.
//  Then, create a seperate cpp file which includes "mesh.h" and
//  explicitly instantiates the template with the desired template
//  parameters.
//  By following this scheme, you can prevent re-compiling the entire
//  template implementation in every usage compilation unit and every
//  time those compilation units are recompiled during development.

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#include <set>
#include <sstream>

#include <map>

#include <cork/accel/aabvh.h>

#include <cork/isct/empty3d.h>
#include <cork/isct/quantization.h>
#include <cork/isct/unsafeRayTriIsct.h>
#include <cork/isct/triangle.h>

#include <cork/math/bbox.h>
#include <cork/math/vec.h>
#include <cork/math/ray.h>

#include <cork/rawmesh/rawMesh.h>

#include <cork/util/iterPool.h>
#include <cork/util/memPool.h>
#include <cork/util/prelude.h>
#include <cork/util/shortVec.h>

#include <cork/util/unionFind.h>

#define REAL double

extern "C" {
    #include <cork/isct/triangle.h>
}


struct BoolVertexData {
};

struct BoolTriangleData {
    byte bool_alg_data; // internal use by algorithm
    // please copy value when the triangle is subdivided
};

template<class VertData, class TriData>
struct IsctVertEdgeTriInput
{
    VertData*   e[2];
    VertData*   t[3];
};

template<class VertData, class TriData>
struct IsctVertTriTriTriInput
{
    VertData*   t[3][3];
};

template<class VertData, class TriData>
struct SubdivideTriInput
{
    TriData*    pt;
    VertData*   pv[3];
    VertData*   v[3];
};

// in order to perform intersections, VertData and TriData must support
struct IsctVertexData {
    // specify how to compute new data for vertices formed by intersections
    /*
    // vertices on edge and triangle forming an intersection...
    void isct(IsctVertEdgeTriInput input);
    void isct(IsctVertTriTriTriInput input);
    void isctInterpolate(const VertData &v0, const VertData &v1);
    */
};

struct IsctTriangleData {
    // specify how to compute new data for a triangle in the event
    // that it is merged with another triangle (merge)
    // split into two triangles (split)
    // or that the triangle is moved (move)
    /*
    void subdivide(SubdivideTriInput input);
    */
};

// in order to perform remeshing, VertData and TriData must support
struct RemeshVertexData {
    bool manifold; // whether this point is manifold.
    // useful for modifying interpolation behavior
    // specify how to compute new data for a vertex in the event of
    // either an edge collapse (via merge) or edge split (via interpolate)
    /*
    void merge(const VertData &v0, const VertData &v1);
    void interpolate(const VertData &v0, const VertData &v1);
    */
};

struct RemeshTriangleData {
    // specify how to compute new data for a triangle in the event
    // that it is merged with another triangle (merge)
    // split into two triangles (split)
    // or that the triangle is moved (move)
    /*
    void merge(const TriData &t0, const TriData &t1);
    static void split(TriData &t0, TriData &t1, const TriData &t_orig);
    void move(const TriData &t_old);
    */
};

struct RemeshOptions
{
    double maxEdgeLength;
    double minEdgeLength;
    double minAngle;
    double maxAngle;

    RemeshOptions() :
            maxEdgeLength(1.0),
            minEdgeLength(0.3),
            minAngle(5.0),
            maxAngle(170.0)
    {}
};

// only for internal use, please do not use as client
struct TopoVert;
struct TopoEdge;
struct TopoTri;

typedef TopoVert* Vptr;
typedef TopoEdge* Eptr;
typedef TopoTri*  Tptr;
//using Vptr = TopoVert*;
//using Eptr = TopoEdge*;
//using Tptr = TopoTri*;
// end internal items

template<class VertData, class TriData>
class Mesh
{
public:
    Mesh();
    Mesh(Mesh &&src);
    Mesh(const RawMesh<VertData,TriData> &raw);
    virtual ~Mesh();

    void operator=(Mesh &&src);

    // validity check:
    //  - all numbers are well-defined and finite
    //  - all triangle vertex indices are in the right range
    bool valid() const;

    RawMesh<VertData,TriData> raw() const;

    inline int numVerts() const { return verts.size(); }
    inline int numTris() const { return tris.size(); }

    inline void for_verts(std::function<void(VertData &)> func);
    inline void for_tris(std::function<void(TriData &, VertData &, VertData &, VertData &)> func);
    inline void for_edges(std::function<void(VertData &, VertData &)> start,
                          std::function<void(TriData &t,VertData &, VertData &, VertData &)> each_tri);

    // form the disjoint union of two meshes
    void disjointUnion(const Mesh &cp);

    struct Isct {
        Ray3d   ray;
        bool    exists;

        uint    tri_id;
        Vec3d   isct;
        Vec3d   bary;
    };

    Isct pick(Ray3d ray);
    inline void accessIsct(const Isct &isct,
                           std::function<void(TriData &,
                                              VertData &, VertData &, VertData &)> func);

    // checks if the mesh is closed
    bool isClosed();

public: // REMESHING module
    // REQUIRES:
    //  - MinimalData
    //  - RemeshData
    void remesh();
    RemeshOptions remesh_options;

public: // ISCT (intersections) module
    void resolveIntersections(); // makes all intersections explicit
    bool isSelfIntersecting(); // is the mesh self-intersecting?

    void testingComputeStaticIsctPoints(std::vector<Vec3d> *points);
    void testingComputeStaticIsct(std::vector<Vec3d> *points,
                                  std::vector< std::pair<Vec3d,Vec3d> > *edges);

public: // BOOLean operation module
    // all of the form
    //      this = this OP rhs
    void boolUnion(Mesh &rhs);
    void boolDiff(Mesh &rhs);
    void boolIsct(Mesh &rhs);
    void boolXor(Mesh &rhs);

private:    // Internal Formats
    struct Tri {
        TriData data;
        union {
            struct {
                uint a, b, c; // vertex ids
            };
            uint v[3];
        };

        inline Tri() {}
    };

    inline void merge_tris(uint tid_result, uint tid0, uint tid1);
    inline void split_tris(uint t0ref, uint t1ref, uint t_orig_ref);
    inline void move_tri(Tri &t_new, Tri &t_old);
    inline void subdivide_tri(uint t_piece_ref, uint t_parent_ref);

private:    // DATA
    std::vector<Tri>        tris;
    std::vector<VertData>   verts;

private:    // caches
    struct NeighborEntry {
        uint vid;
        ShortVec<uint, 2> tids;
        inline NeighborEntry() {}
        inline NeighborEntry(uint vid_) : vid(vid_) {}
    };

    struct NeighborCache {
        std::vector< ShortVec<NeighborEntry, 8> > skeleton;
        inline NeighborEntry& operator()(uint i, uint j) {
            uint N = skeleton[i].size();
            for(uint k = 0; k < N; k++) {
                if(skeleton[i][k].vid == j)
                    return skeleton[i][k];
            }
            skeleton[i].push_back(NeighborEntry(j));
            return skeleton[i][N];
        }
    };

    NeighborCache createNeighborCache();

    // parallel to vertex array
    std::vector<uint> getComponentIds();

    // like the neighbor cache, but more customizable
    template<class Edata>
    struct EGraphEntry {
        uint                vid;
        ShortVec<uint, 2>   tids;
        Edata               data;
        inline EGraphEntry() {}
        inline EGraphEntry(uint vid_) : vid(vid_) {}
    };
    template<class Edata>
    struct EGraphCache {
        std::vector< ShortVec<EGraphEntry<Edata>, 8> > skeleton;
        inline EGraphEntry<Edata> & operator()(uint i, uint j) {
            uint N = skeleton[i].size();
            for(uint k = 0; k < N; k++) {
                if(skeleton[i][k].vid == j)
                    return skeleton[i][k];
            }
            skeleton[i].push_back(EGraphEntry<Edata>(j));
            return skeleton[i][N];
        }
        inline void for_each(std::function<void(
                uint i, uint j, EGraphEntry<Edata> &entry
        )> action
        ) {
            for(uint i=0; i<skeleton.size(); i++) {
                for(auto &entry : skeleton[i]) {
                    action(i, entry.vid, entry);
                }
            }
        }
    };
    template<class Edata>
    EGraphCache<Edata> createEGraphCache();


private:    // TopoCache Support
    struct TopoCache;
private:    // Isct Support
    class  IsctProblem; // implements intersection functionality
    class TriangleProblem; // support type for IsctProblem
    typedef TriangleProblem* Tprob;
    //using Tprob = TriangleProblem*;
private:    // Bool Support
    class BoolProblem;

private:    // Remeshing Support
    struct RemeshScratchpad;

    Eptr allocateRemeshEdge(RemeshScratchpad &);
    void deallocateRemeshEdge(RemeshScratchpad &, Eptr);

    void edgeSplit(RemeshScratchpad &,
                   Eptr e_split);
    void edgeCollapse(RemeshScratchpad &,
                      Eptr e_collapse,
                      bool collapsing_tetrahedra_disappear);

    // Need edge scoring routines...
    void scoreAndEnqueue(std::set< std::pair<double, Eptr> > &queue, Eptr edge);
    void dequeue(std::set< std::pair<double, Eptr> > &queue, Eptr edge);
    double computeEdgeScore(Eptr edge);

    // support functions
    void populateTriFromTopoTri(Tptr t);
    // calls the first function once, then the second once for each triangle
    inline void edgeNeighborhood(
            Eptr edge,
            std::function<void(VertData &v0, VertData &v1)> once,
            std::function<void(VertData &v0, VertData &v1,
                               VertData &vopp, TriData &t)> each_tri
    );
};

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::for_verts(std::function<void(VertData &v)> func) {
    for(auto &v : verts)
        func(v);
}

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::for_tris(std::function<void(TriData &, VertData &, VertData &, VertData &)> func) {
    for(auto &tri : tris) {
        auto &a = verts[tri.a];
        auto &b = verts[tri.b];
        auto &c = verts[tri.c];
        func(tri.data, a, b, c);
    }
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::for_edges(std::function<void(VertData &, VertData &)> start,
                                       std::function<void(TriData &t, VertData &, VertData &, VertData &)> each_tri
) {
    NeighborCache cache = createNeighborCache();
    for(uint i=0; i<cache.skeleton.size(); i++) {
        for(auto &entry : cache.skeleton[i]) {
            uint j = entry.vid;
            start(verts[i], verts[j]);
            for(uint tid : entry.tids) {
                Tri &tri = tris[tid];
                each_tri(tri.data, verts[tri.a], verts[tri.b], verts[tri.c]);
            }
        }
    }
}

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::accessIsct(const Isct &isct,
                                               std::function<void(TriData &, VertData &, VertData &, VertData &)> func) {
    Tri &tri = tris[isct.tri_id];
    auto &a = verts[tri.a];
    auto &b = verts[tri.b];
    auto &c = verts[tri.c];
    func(tri.data, a, b, c);
}


/*
 * Implementation of mesh
 */
// constructors
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh() {}
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh(Mesh &&cp)
        : tris(cp.tris), verts(cp.verts)
{}
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh(const RawMesh<VertData,TriData> &raw) :
        tris(raw.triangles.size()), verts(raw.vertices)
{
    // fill out the triangles
    for(uint i=0; i<raw.triangles.size(); i++) {
        tris[i].data = raw.triangles[i];
        tris[i].a = raw.triangles[i].a;
        tris[i].b = raw.triangles[i].b;
        tris[i].c = raw.triangles[i].c;
    }
}
template<class VertData, class TriData>
Mesh<VertData,TriData>::~Mesh()
{

}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::operator=(Mesh &&src)
{
    tris = src.tris;
    verts = src.verts;
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::valid() const
{
    for(uint i=0; i<verts.size(); i++) {
        if(!std::isfinite(verts[i].pos.x) ||
           !std::isfinite(verts[i].pos.y) ||
           !std::isfinite(verts[i].pos.z)) {
            std::ostringstream message;
            message << "vertex #" << i << " has non-finite coordinates: "
                    << verts[i].pos;
            CORK_ERROR(message.str());
            return false;
        }
    }

    for(uint i=0; i<tris.size(); i++) {
        if(tris[i].a >= verts.size() ||
           tris[i].b >= verts.size() ||
           tris[i].c >= verts.size()) {
            std::ostringstream message;
            message << "triangle #" << i << " should have indices in "
                    << "the range 0 to " << (verts.size()-1)
                    << ", but it has invalid indices: "
                    << tris[i].a << ", " << tris[i].b << ", " << tris[i].c;
            CORK_ERROR(message.str());
            return false;
        }
    }

    return true;
}

template<class VertData, class TriData>
RawMesh<VertData,TriData> Mesh<VertData,TriData>::raw() const
{
    RawMesh<VertData,TriData> result;
    result.vertices = verts;
    result.triangles.resize(tris.size());
    for(uint i=0; i<tris.size(); i++) {
        result.triangles[i]   = tris[i].data;
        result.triangles[i].a = tris[i].a;
        result.triangles[i].b = tris[i].b;
        result.triangles[i].c = tris[i].c;
    }
    return result;
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::disjointUnion(const Mesh &cp)
{
    uint oldVsize = verts.size();
    uint oldTsize = tris.size();
    uint cpVsize  = cp.verts.size();
    uint cpTsize  = cp.tris.size();
    uint newVsize = oldVsize + cpVsize;
    uint newTsize = oldTsize + cpTsize;

    std::vector<int> v_remap(cpVsize); // oh this is obvious...
    verts.resize(newVsize);
    tris.resize(newTsize);

    for(uint i=0; i<cpVsize; i++)
        verts[oldVsize + i] = cp.verts[i];

    for(uint i=0; i<cpTsize; i++) {
        auto &tri = tris[oldTsize + i];
        tri = cp.tris[i];
        tri.a += oldVsize;
        tri.b += oldVsize;
        tri.c += oldVsize;
    }
}

// Picking.
// Dumb Implementation just passes over all triangles w/o any precomputed
// acceleration structure
template<class VertData, class TriData>
typename Mesh<VertData,TriData>::Isct
Mesh<VertData,TriData>::pick(Ray3d ray)
{
    Isct result;
    result.ray = ray;
    result.exists = false;

    double mint = DBL_MAX;

    // pass all triangles over ray
    for(uint i=0; i<tris.size(); i++) {
        const Tri  &tri = tris[i];

        uint   a = tri.a;
        uint   b = tri.b;
        uint   c = tri.c;
        Vec3d va = verts[a].pos;
        Vec3d vb = verts[b].pos;
        Vec3d vc = verts[c].pos;
        // normalize vertex order (to prevent leaks)
        if(a > b) { std::swap(a, b); std::swap(va, vb); }
        if(b > c) { std::swap(b, c); std::swap(vb, vc); }
        if(a > b) { std::swap(a, b); std::swap(va, vb); }

        double t;
        Vec3d  bary;
        if(isct_ray_triangle(ray, va, vb, vc, &t, &bary)) {
            if(t > 0 && t < mint) {
                result.exists = true;
                mint = t;
                result.tri_id = i;
                result.isct = ray.p + t * ray.r;
                result.bary = bary;
            }
        }
    }

    return result;
}




template<class VertData, class TriData>
bool Mesh<VertData,TriData>::isClosed()
{
    EGraphCache<int> chains = createEGraphCache<int>();
    chains.for_each([&](uint i, uint j, EGraphEntry<int> &entry) {
        entry.data = 0;
    });
    // count up how many times each edge is encountered in one
    // orientation vs. the other
    for(Tri &tri : tris) {
        chains(tri.a, tri.b).data ++;
        chains(tri.b, tri.a).data --;

        chains(tri.b, tri.c).data ++;
        chains(tri.c, tri.b).data --;

        chains(tri.c, tri.a).data ++;
        chains(tri.a, tri.c).data --;
    }
    // now go through and see if any of these are non-zero
    bool closed = true;
    chains.for_each([&](uint i, uint j, EGraphEntry<int> &entry) {
        if(entry.data != 0)
            closed = false;
    });
    return closed;
}




static inline
bool contains(const ShortVec<uint, 8> &list, uint item)
{
    for(uint k : list)
        if(k == item)
            return true;
    return false;
}

template<class VertData, class TriData>
typename Mesh<VertData,TriData>::NeighborCache
Mesh<VertData,TriData>::createNeighborCache()
{
    NeighborCache result;
    result.skeleton.resize(verts.size());

    for(uint tid = 0; tid < tris.size(); tid++) {
        const Tri &tri = tris[tid];

        result(tri.a, tri.b).tids.push_back(tid);
        result(tri.b, tri.a).tids.push_back(tid);

        result(tri.a, tri.c).tids.push_back(tid);
        result(tri.c, tri.a).tids.push_back(tid);

        result(tri.b, tri.c).tids.push_back(tid);
        result(tri.c, tri.b).tids.push_back(tid);
    }

    return result;
}

// This function signature is an amazing disaster...
#ifdef _WIN32
template<class VertData, class TriData>
template<class Edata>
typename Mesh<VertData,TriData>::EGraphCache<Edata>
#else
template<class VertData, class TriData>
template<class Edata>
typename Mesh<VertData,TriData>::template EGraphCache<Edata>
#endif
Mesh<VertData,TriData>::createEGraphCache()
{
    EGraphCache<Edata> result;
    result.skeleton.resize(verts.size());

    for(uint tid = 0; tid < tris.size(); tid++) {
        const Tri &tri = tris[tid];

        result(tri.a, tri.b).tids.push_back(tid);
        result(tri.b, tri.a).tids.push_back(tid);

        result(tri.a, tri.c).tids.push_back(tid);
        result(tri.c, tri.a).tids.push_back(tid);

        result(tri.b, tri.c).tids.push_back(tid);
        result(tri.c, tri.b).tids.push_back(tid);
    }

    return result;
}


template<class VertData, class TriData>
std::vector<uint> Mesh<VertData,TriData>::getComponentIds()
{
    UnionFind uf(verts.size());
    for(const Tri &tri : tris) {
        uf.unionIds(tri.a, tri.b);
        uf.unionIds(tri.a, tri.c);
    }

    return uf.dump();
}

#include <cork/util/iterPool.h>

/*
 *  Allows for topological algorithms to manipulate
 *  a more familiar pointer data structure based on a simplicial complex.
 *  This structure can be regenerated from the more basic
 *  vertex/triangle arrays using
 *      createTopoCache()
 *  Once manipulations have been satisfactorily performed,
 *  the underlying vertex/triangle arrays can be cleaned up for
 *  further use by topologically insensitive algorithms by
 *      commitTopoCache()
 */

#define INVALID_ID uint(-1)

struct TopoVert {
    uint                    ref;        // index to actual data
    void*                   data;       // algorithm specific handle

    ShortVec<Tptr, 8>       tris;       // triangles this vertex is incident on
    ShortVec<Eptr, 8>       edges;      // edges this vertex is incident on
};

struct TopoEdge {
    void*                   data;       // algorithm specific handle

    Vptr                    verts[2];   // endpoint vertices
    ShortVec<Tptr, 2>       tris;       // incident triangles
};

struct TopoTri {
    uint                    ref;        // index to actual data
    void*                   data;       // algorithm specific handle

    Vptr                    verts[3];   // vertices of this triangle
    Eptr                    edges[3];   // edges of this triangle
                                        // opposite to the given vertex
};


template<class VertData, class TriData>
struct Mesh<VertData, TriData>::TopoCache {
    IterPool<TopoVert>    verts;
    IterPool<TopoEdge>    edges;
    IterPool<TopoTri>     tris;

    Mesh *mesh;
    TopoCache(Mesh *owner);
    virtual ~TopoCache() {}

    // until commit() is called, the Mesh::verts and Mesh::tris
    // arrays will still contain garbage entries
    void commit();

    bool isValid();
    void print();

    // helpers to create bits and pieces
    inline Vptr newVert();
    inline Eptr newEdge();
    inline Tptr newTri();

    // helpers to release bits and pieces
    inline void freeVert(Vptr);
    inline void freeEdge(Eptr);
    inline void freeTri(Tptr);

    // helper to delete geometry in a structured way
    inline void deleteTri(Tptr);

    // helper to flip triangle orientation
    inline void flipTri(Tptr);

private:
    void init();
};


template<class VertData, class TriData> inline
Vptr Mesh<VertData, TriData>::TopoCache::newVert()
{
    uint        ref         = mesh->verts.size();
                mesh->verts.push_back(VertData());
    Vptr        v           = verts.alloc(); // cache.verts
                v->ref      = ref;
                return v;
}
template<class VertData, class TriData> inline
Eptr Mesh<VertData, TriData>::TopoCache::newEdge()
{
    Eptr        e           = edges.alloc(); // cache.edges
                return e;
}
template<class VertData, class TriData> inline
Tptr Mesh<VertData, TriData>::TopoCache::newTri()
{
    uint        ref         = mesh->tris.size();
                mesh->tris.push_back(Tri());
    Tptr        t           = tris.alloc(); // cache.tris
                t->ref      = ref;
                return t;
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::freeVert(Vptr v)
{
    verts.free(v);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::freeEdge(Eptr e)
{
    edges.free(e);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::freeTri(Tptr t)
{
    tris.free(t);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::deleteTri(Tptr tri)
{
    // first, unhook the triangle from its faces
    for(uint k=0; k<3; k++) {
        Vptr            v                   = tri->verts[k];
                        v->tris.erase(tri);
        Eptr            e                   = tri->edges[k];
                        e->tris.erase(tri);
    }
    // now, let's check for any edges which no longer border triangles
    for(uint k=0; k<3; k++) {
        Eptr            e                   = tri->edges[k];
        if(e->tris.size() == 0) {
            // delete edge
            // unhook from vertices
            Vptr        v0                  = e->verts[0];
                        v0->edges.erase(e);
            Vptr        v1                  = e->verts[1];
                        v1->edges.erase(e);
            freeEdge(e);
        }
    }
    // now, let's check for any vertices which no longer border triangles
    for(uint k=0; k<3; k++) {

        Vptr v = tri->verts[k];

        if(v->tris.size() == 0) {
            freeVert(v);
        }

    }

    // finally, release the triangle
    freeTri(tri);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::flipTri(Tptr t)
{
    std::swap(t->verts[0], t->verts[1]);
    std::swap(t->edges[0], t->edges[1]);
    std::swap(mesh->tris[t->ref].v[0], mesh->tris[t->ref].v[1]);
}

template<class VertData, class TriData>
Mesh<VertData, TriData>::TopoCache::TopoCache(Mesh *owner) : mesh(owner)
{
    init();
}


// support structure for cache construction
struct TopoEdgePrototype {
    uint vid;
    ShortVec<Tptr, 2> tris;
    TopoEdgePrototype() {}
    TopoEdgePrototype(uint v) : vid(v) {}
};

inline TopoEdgePrototype& getTopoEdgePrototype(uint a, uint b,
                                               std::vector< ShortVec<TopoEdgePrototype, 8> > &prototypes) {
    uint N = prototypes[a].size();

    for(uint i=0; i<N; i++) {
        if(prototypes[a][i].vid == b)
            return prototypes[a][i];
    }
    prototypes[a].push_back(TopoEdgePrototype(b));
    return prototypes[a][N];
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::TopoCache::init()
{
    // first lay out vertices
    std::vector<Vptr> temp_verts(mesh->verts.size()); // need temp. reference
    for(uint i=0; i<mesh->verts.size(); i++) {
        Vptr vert = verts.alloc(); // cache.verts.alloc()
        vert->ref = i;
        temp_verts[i] = vert;
    }

    // We need to still do the following
    //  * Generate TopoTris
    //  * Generate TopoEdges
    // ---- Hook up references between
    //  * Triangles and Vertices
    //  * Triangles and Edges
    //  * Vertices and Edges

    // We handle two of these items in a pass over the triangles,
    //  * Generate TopoTris
    //  * Hook up Triangles and Vertices
    // building a structure to handle the edges as we go:
    std::vector< ShortVec<TopoEdgePrototype, 8> > edgeacc(mesh->verts.size());
    for(uint i=0; i<mesh->tris.size(); i++) {
        Tptr tri = tris.alloc(); // cache.tris.alloc()
        tri->ref = i;
        const Tri &ref_tri = mesh->tris[i];

        // triangles <--> verts
        uint vids[3];
        for(uint k=0; k<3; k++) {
            uint vid = vids[k] = ref_tri.v[k];
            tri->verts[k] = temp_verts[vid];
            temp_verts[vid]->tris.push_back(tri);
        }
        // then, put these in arbitrary but globally consistent order
        if(vids[0] > vids[1])   std::swap(vids[0], vids[1]);
        if(vids[1] > vids[2])   std::swap(vids[1], vids[2]);
        if(vids[0] > vids[1])   std::swap(vids[0], vids[1]);
        // and accrue in structure
        getTopoEdgePrototype(vids[0], vids[1], edgeacc).tris.push_back(tri);
        getTopoEdgePrototype(vids[0], vids[2], edgeacc).tris.push_back(tri);
        getTopoEdgePrototype(vids[1], vids[2], edgeacc).tris.push_back(tri);
    }

    // Now, we can unpack the edge accumulation to
    //  * Generate TopoEdges
    //  * Hook up Triangles and Edges
    //  * Hook up Vertices and Edges
    for(uint vid0=0; vid0 < edgeacc.size(); vid0++) {
      for(TopoEdgePrototype &proto : edgeacc[vid0]) {
        uint vid1 = proto.vid;
        Vptr v0 = temp_verts[vid0];
        Vptr v1 = temp_verts[vid1];

        Eptr edge = edges.alloc(); // cache.edges.alloc()
        // edges <--> verts
        edge->verts[0] = v0;
        v0->edges.push_back(edge);
        edge->verts[1] = v1;
        v1->edges.push_back(edge);
        // edges <--> tris
        for(Tptr tri : proto.tris) {
            edge->tris.push_back(tri);
            for(uint k=0; k<3; k++) {
                if(v0 != tri->verts[k] && v1 != tri->verts[k]) {
                    tri->edges[k] = edge;
                    break;
                }
            }
        }
    }}

    //ENSURE(isValid());
    //print();
}




template<class VertData, class TriData>
void Mesh<VertData, TriData>::TopoCache::commit()
{
    //ENSURE(isValid());

    // record which vertices are live
    std::vector<bool> live_verts(mesh->verts.size(), false);
    verts.for_each([&](Vptr vert) { // cache.verts
        live_verts[vert->ref] = true;
    });

    // record which triangles are live, and record connectivity
    std::vector<bool> live_tris(mesh->tris.size(), false);
    tris.for_each([&](Tptr tri) { // cache.tris
        live_tris[tri->ref] = true;
        for(uint k=0; k<3; k++)
            mesh->tris[tri->ref].v[k] = tri->verts[k]->ref;
    });

    // compact the vertices and build a remapping function
    std::vector<uint> vmap(mesh->verts.size());
    uint write = 0;
    for(uint read = 0; read < mesh->verts.size(); read++) {
        if(live_verts[read]) {
            vmap[read] = write;
            mesh->verts[write] = mesh->verts[read];
            write++;
        } else {
            vmap[read] = INVALID_ID;
        }
    }
    mesh->verts.resize(write);

    // rewrite the vertex reference ids
    verts.for_each([&](Vptr vert) { // cache.verts
        vert->ref = vmap[vert->ref];
    });

    std::vector<uint> tmap(mesh->tris.size());
    write = 0;
    for(uint read = 0; read < mesh->tris.size(); read++) {
        if(live_tris[read]) {
            tmap[read] = write;
            mesh->tris[write] = mesh->tris[read];
            for(uint k=0; k<3; k++)
                mesh->tris[write].v[k] = vmap[mesh->tris[write].v[k]];
            write++;
        } else {
            tmap[read] = INVALID_ID;
        }
    }
    mesh->tris.resize(write);

    // rewrite the triangle reference ids
    tris.for_each([&](Tptr tri) { // cache.tris
        tri->ref = tmap[tri->ref];
    });
}



// support functions for validity check
template<class T, class Container> inline
bool count(const Container &contain, const T &val) {
    uint c=0;
    for(const T &t : contain)
        if(t == val)    c++;
    return c;
}
template<class T> inline
bool count2(const T arr[], const T &val) {
    return ((arr[0] == val)? 1 : 0) + ((arr[1] == val)? 1 : 0);
}
template<class T> inline
bool count3(const T arr[], const T &val) {
    return ((arr[0] == val)? 1 : 0) + ((arr[1] == val)? 1 : 0)
                                    + ((arr[2] == val)? 1 : 0);
}

template<class VertData, class TriData>
bool Mesh<VertData, TriData>::TopoCache::isValid()
{
    //print();
    std::set<Vptr> vaddr;
    std::set<Eptr> eaddr;
    std::set<Tptr> taddr;
    verts.for_each([&vaddr](Vptr v) { vaddr.insert(v); });
    edges.for_each([&eaddr](Eptr e) { eaddr.insert(e); });
    tris.for_each( [&taddr](Tptr t) { taddr.insert(t); });

    // check verts
    verts.for_each([&](Vptr v) {
        ENSURE(v->ref < mesh->verts.size());
        // make sure each edge pointer goes somewhere and that
        // the pointed-to site also points back correctly
        for(Eptr e : v->edges) {
            ENSURE(eaddr.count(e) > 0); // pointer is good
            ENSURE(count2(e->verts, v) == 1); // back-pointer is good
        }
        for(Tptr t : v->tris) {
            ENSURE(taddr.count(t) > 0);
            ENSURE(count3(t->verts, v) == 1);
        }
    });

    // check edges
    edges.for_each([&](Eptr e) {
        // check for non-degeneracy
        ENSURE(e->verts[0] != e->verts[1]);
        for(uint k=0; k<2; k++) {
            Vptr v = e->verts[k];
            ENSURE(vaddr.count(v) > 0);
            ENSURE(count(v->edges, e) == 1);
        }
        for(Tptr t : e->tris) {
            ENSURE(taddr.count(t) > 0);
            ENSURE(count3(t->edges, e) == 1);
        }
    });

    // check triangles
    tris.for_each([&](Tptr t) {
        // check for non-degeneracy
        ENSURE(t->verts[0] != t->verts[1] && t->verts[1] != t->verts[2]
                                          && t->verts[0] != t->verts[2]);
        for(uint k=0; k<3; k++) {
            Vptr v = t->verts[k];
            ENSURE(vaddr.count(v) > 0);
            ENSURE(count(v->tris, t) == 1);

            Eptr e = t->edges[k];
            ENSURE(eaddr.count(e) == 1);
            ENSURE(count(e->tris, t) == 1);

            // also need to ensure that the edges are opposite the
            // vertices as expected
            Vptr v0 = e->verts[0];
            Vptr v1 = e->verts[1];
            ENSURE((v0 == t->verts[(k+1)%3] && v1 == t->verts[(k+2)%3])
                || (v0 == t->verts[(k+2)%3] && v1 == t->verts[(k+1)%3]));
        }
    });

    return true;
}


CORK_EXPORT std::ostream& operator<<(std::ostream &out, const TopoVert& vert);
CORK_EXPORT std::ostream& operator<<(std::ostream &out, const TopoEdge& edge);
CORK_EXPORT std::ostream& operator<<(std::ostream &out, const TopoTri& tri);

template<class VertData, class TriData>
void Mesh<VertData, TriData>::TopoCache::print()
{
    using std::cout;
    using std::endl;

    cout << "dumping remeshing cache for debug..." << endl;
    cout << "TRIS" << endl;
    int tri_count = 0;
    tris.for_each([&](Tptr t) {
        cout << " " << t << ": " << *t << endl;
        tri_count++;
    });
    cout << "There were " << tri_count << " TRIS" << endl;
    cout << "EDGES" << endl;
    int edge_count = 0;
    edges.for_each([&](Eptr e) {
        cout << " " << e << ": " << endl;
        cout << "  v " << e->verts[0] << "; "
                       << e->verts[1] << endl;
        cout << "  t (" << e->tris.size() << ")" << endl;
        for(Tptr t : e->tris)
        cout << "    " << t << endl;
        edge_count++;
    });
    cout << "There were " << edge_count << " EDGES" << endl;
    cout << "VERTS" << endl;
    int vert_count = 0;
    verts.for_each([&](Vptr v) {
        cout << " " << v << ": ref(" << v->ref << ")" << endl;
        cout << "  e (" << v->edges.size() << ")" << endl;
        for(Eptr e : v->edges)
        cout << "    " << e << endl;
        cout << "  t (" << v->tris.size() << ")" << endl;
        for(Tptr t : v->tris)
        cout << "    " << t << endl;
        vert_count++;
    });
    cout << "There were " << vert_count << " VERTS" << endl;
}


// Implementation

#include "mesh.remesh.tpp"
#include "mesh.isct.tpp"
#include "mesh.bool.tpp"


#endif

