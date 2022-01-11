#include "mesh.h"

std::ostream& operator<<(std::ostream &out, const TopoVert& vert)
{
    out << "ref(" << vert.ref << ") "
        << "e(" << vert.edges.size() << "):";
    for(Eptr e : vert.edges)
        out << e << ";";
    out << " "
        << "t(" << vert.tris.size() << "):";
    for(Tptr t : vert.tris)
        out << t << ";";
    return out;
}

std::ostream& operator<<(std::ostream &out, const TopoEdge& edge)
{
    out << "v(2):" << edge.verts[0] << "(" << edge.verts[0]->ref << ");"
                   << edge.verts[1] << "(" << edge.verts[1]->ref << ");";
    out << " "
        << "t(" << edge.tris.size() << "):";
    for(Tptr t : edge.tris)
        out << t << ";";
    return out;
}

std::ostream& operator<<(std::ostream &out, const TopoTri& tri)
{
    out << "ref(" << tri.ref << ") ";
    out << "v(3):" << tri.verts[0] << "(" << tri.verts[0]->ref << ");"
                   << tri.verts[1] << "(" << tri.verts[1]->ref << ");"
                   << tri.verts[2] << "(" << tri.verts[2]->ref << ");";
    out << " ";
    out << "e(3):" << tri.edges[0] << ";"
                   << tri.edges[1] << ";"
                   << tri.edges[2] << ";";
    return out;
}
