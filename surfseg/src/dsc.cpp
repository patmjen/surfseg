#include "dsc.h"

float DSC::tetQuality(TetKey tk) const
{
    // Computes volume length ratio
    constexpr float oneDivSix = 1.0f / 6.0f;
    constexpr float sqrtTwo = static_cast<float>(1.4142135623730950488016887242097);
    auto vs = getTetVerts(tk);
    Vec3f p4 = vertices[vs[0]].pos;
    Vec3f p1 = vertices[vs[1]].pos;
    Vec3f p2 = vertices[vs[2]].pos;
    Vec3f p3 = vertices[vs[3]].pos;    
    float vol = fabsf(dot(p1 - p4, cross(p2 - p4, p3 - p4)))*oneDivSix;
    float lRms = 0.0f;
    for (EdgeKey ek : getTetEdges(tk)) {
        Edge e = edges[ek];
        lRms += sqr_length(vertices[e.v1].pos - vertices[e.v0].pos);
    }
    lRms = sqrt(lRms*oneDivSix);
    return 6.0f*sqrtTwo*vol / (lRms*lRms*lRms);
}

std::vector<float> DSC::allTetQuality() const
{
    std::vector<float> out;
    out.reserve(tets.size());
    for (const Tetrahedron& t : tets) {
        out.push_back(tetQuality(t.self));
    }
    return out;
}

void DSC::smartLaplacianSmooth(bool interfaceFlagsSet)
{
    auto tetQ = allTetQuality();
    smartLaplacianSmooth(tetQ, interfaceFlagsSet);
}

void DSC::smartLaplacianSmooth(std::vector<float>& tetQualities, bool interfaceFlagsSet)
{
    // Smart Laplacian smoothing of non-interface vertices
    if (!interfaceFlagsSet) {
        setInterfaceFlags();
    }
    for (Vertex& v : vertices) {
        if (v.interf) {
            // Skip all interface edges
            //continue;
        }
        std::vector<TetKey> coTets;
        std::vector<EdgeKey> coEdges;
        std::tie(coEdges, std::ignore, coTets) = getVertCoboundary(v.self);
        float oldMinQual = INFINITY;
        for (TetKey tk : coTets) {
            oldMinQual = (tetQualities[tk] < oldMinQual) ? tetQualities[tk] : oldMinQual;
        }
        Vec3f newPos(0.0f);
        for (EdgeKey ek : coEdges) {
            const Edge& e = edges[ek];
            newPos += (edges[ek].v0 == v.self) ? vertices[e.v1].pos : vertices[e.v0].pos;
        }
        newPos /= coEdges.size();
        Vec3f oldPos = v.pos;
        v.pos = newPos;
        float newMinQual = INFINITY;
        for (TetKey tk : coTets) {
            float q = tetQuality(tk);
            newMinQual = (q < newMinQual) ? q : newMinQual;
        }
        if (newMinQual > oldMinQual) {
            // Keep new position and update tet. qualities
            for (TetKey tk : coTets) {
                tetQualities[tk] = tetQuality(tk);
            }
        } else {
            // Revert position
            v.pos = oldPos;
        }
    }
}
