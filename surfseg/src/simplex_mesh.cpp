#include "simplex_mesh.h"
#include <algorithm>
#include <fstream>
#include <queue>

#include <GL/glew.h>

#include "util.h"

SimplexMesh::SimplexMesh(const SimplexMesh& other) :
    vertices(other.vertices),
    edges(other.edges),
    triangles(other.triangles),
    tets(other.tets)
{}

SimplexMesh::SimplexMesh(SimplexMesh&& other) :
    vertices(std::move(other.vertices)),
    edges(std::move(other.edges)),
    triangles(std::move(other.triangles)),
    tets(std::move(other.tets))
{}

SimplexMesh& SimplexMesh::operator=(const SimplexMesh& other)
{
    if (this != &other) {
        vertices = other.vertices;
        edges = other.edges;
        triangles = other.triangles;
        tets = other.tets;
    }
    return *this;
}

SimplexMesh& SimplexMesh::operator=(SimplexMesh&& other)
{
    if (this != &other) {
        vertices = std::move(other.vertices);
        edges = std::move(other.edges);
        triangles = std::move(other.triangles);
        tets = std::move(other.tets);
    }
    return *this;
}

void SimplexMesh::clear()
{
    vertices.clear();
    edges.clear();
    triangles.clear();
    tets.clear();
}

void SimplexMesh::build(const std::vector<Vec3f>& vertexPositions,
    const std::vector<std::array<int, 4>>& tetIdxs, const std::vector<int> *tetValues)
{
    clear();
    std::unordered_map<std::pair<int, int>, EdgeKey> edgeMap;
    std::unordered_map<std::array<int, 3>, TriKey> triMap;
    std::unordered_map<std::array<int, 4>, TetKey> tetMap;

    // Add vertices
    for (const Vec3f& p : vertexPositions) {
        Vertex v;
        v.pos = p;
        v.self = vertices.size();
        vertices.push_back(v);
    }

    int i = 0;
    for (const auto& idxs : tetIdxs) {
        // Add edges
        addEdge_(idxs[0], idxs[1], edgeMap);
        addEdge_(idxs[0], idxs[2], edgeMap);
        addEdge_(idxs[0], idxs[3], edgeMap);
        addEdge_(idxs[1], idxs[2], edgeMap);
        addEdge_(idxs[1], idxs[3], edgeMap);
        addEdge_(idxs[2], idxs[3], edgeMap);

        // Add triangles
        addTriangle_(idxs[0], idxs[1], idxs[2], triMap, edgeMap);
        addTriangle_(idxs[0], idxs[1], idxs[3], triMap, edgeMap);
        addTriangle_(idxs[0], idxs[2], idxs[3], triMap, edgeMap);
        addTriangle_(idxs[1], idxs[2], idxs[3], triMap, edgeMap);

        // Add tetrahedron
        TetKey tk = addTetrahedron_(idxs[0], idxs[1], idxs[2], idxs[3], tetMap, triMap);
        if (tetValues != nullptr) {
            tets[tk].value = (*tetValues)[i];
        }
        ++i;
    }
}

SimplexMesh::EdgeKey SimplexMesh::addEdge_(int i1, int i2,
    std::unordered_map<std::pair<int, int>, EdgeKey>& edgeMap)
{
    // Ensure vertices are sorted
    int a = (i1 < i2) ? i1 : i2;
    int b = (i1 < i2) ? i2 : i1;

    std::pair<int, int> mapKey = std::make_pair(a, b);
    auto it = edgeMap.find(mapKey);
    EdgeKey ek;
    if (it == edgeMap.end()) {
        // Edge has not been added previosly so add it now
        ek = edges.size();
        Edge e;
        e.self = ek;
        e.v0 = a;
        e.v1 = b;

        // Register new edge with its vertices - it doesn't matter that we overwrite every time
        vertices[a].edge = ek;
        vertices[b].edge = ek;

        edgeMap[mapKey] = ek;
        edges.push_back(e);
    } else {
        ek = it->second;
    }
    return ek;
}

SimplexMesh::TriKey SimplexMesh::addTriangle_(int i1, int i2, int i3,
    std::unordered_map<std::array<int, 3>, TriKey>& triMap,
    const std::unordered_map<std::pair<int, int>, EdgeKey>& edgeMap)
{
    std::array<int, 3> mapKey = { i1, i2, i3 };
    std::sort(mapKey.begin(), mapKey.end());

    auto it = triMap.find(mapKey);
    TriKey tk;
    if (it == triMap.end()) {
        tk = triangles.size();
        Triangle t;
        t.self = tk;

        // We know all edges for this triangle have been added and mapKey is already sorted
        t.edges = {
            edgeMap.at(std::make_pair(mapKey[0], mapKey[1])),
            edgeMap.at(std::make_pair(mapKey[0], mapKey[2])),
            edgeMap.at(std::make_pair(mapKey[1], mapKey[2]))
        };

        // Register new triangle with its edges - it doesn't matter that we overwrite every time
        edges[t.edges[0]].triangle = tk;
        edges[t.edges[1]].triangle = tk;
        edges[t.edges[2]].triangle = tk;

        triMap[mapKey] = tk;
        triangles.push_back(t);
    } else {
        tk = it->second;
    }
    return tk;
}

SimplexMesh::TetKey SimplexMesh::addTetrahedron_(int i1, int i2, int i3, int i4,
    std::unordered_map<std::array<int, 4>, TetKey>& tetMap,
    const std::unordered_map<std::array<int, 3>, TriKey>& triMap)
{
    std::array<int, 4> mapKey = { i1, i2, i3, i4 };
    std::sort(mapKey.begin(), mapKey.end());

    // Small helper function
    const auto setTriangleTet = [&](TriKey triKey, TetKey tetKey) {
        Triangle& tri = triangles[triKey];
        if (tri.t0 == INVALID_TETRAHEDRON) {
            tri.t0 = tetKey;
        } else {
            tri.t1 = tetKey;
        }
    };

    auto it = tetMap.find(mapKey);
    TetKey tk;
    if (it == tetMap.end()) {
        tk = tets.size();
        Tetrahedron t;
        t.self = tk;

        // We know all edges for this tet. have been added and mapKey is already sorted
        t.tris = {
            triMap.at({ mapKey[0], mapKey[1], mapKey[2] }),
            triMap.at({ mapKey[0], mapKey[1], mapKey[3] }),
            triMap.at({ mapKey[0], mapKey[2], mapKey[3] }),
            triMap.at({ mapKey[1], mapKey[2], mapKey[3] })
        };

        // Register new tet. with its triangles - it doesn't matter that we overwrite every time
        setTriangleTet(t.tris[0], tk);
        setTriangleTet(t.tris[1], tk);
        setTriangleTet(t.tris[2], tk);
        setTriangleTet(t.tris[3], tk);

        tetMap[mapKey] = tk;
        tets.push_back(t);
    } else {
        tk = it->second;
    }
    return tk;
}

float SimplexMesh::computeBoundingSphere() const noexcept
{
    Vec3f tmp;
    return computeBoundingSphere(tmp);
}

float SimplexMesh::computeBoundingSphere(Vec3f& center) const noexcept
{
    // Compute centroid
    center = Vec3f(0);
    for (const Vertex& v : vertices) {
        center += v.pos;
    }
    center /= vertices.size();

    // Compute max. distance to centroid
    float r = 0;
    for (const Vertex& v : vertices) {
        r = std::max(r, length(v.pos - center));
    }
    return r;
}

bool SimplexMesh::isTriInterface(TriKey tk, int insideVal) const
{
    return isTriInterface(triangles[tk], insideVal);
}

bool SimplexMesh::isTriInterface(const Triangle& t, int insideVal) const
{
    // TODO: Make this more robust
    int val1 = t.t0 != INVALID_TETRAHEDRON ? tets[t.t0].value : 0;
    int val2 = t.t1 != INVALID_TETRAHEDRON ? tets[t.t1].value : val1;
    return (val1 == insideVal || val2 == insideVal) && val1 != val2;
}

void SimplexMesh::setInterfaceFlags()
{
    for (Triangle& t : triangles) {
        bool interf = isTriInterface(t);
        if (interf) {
            t.interf = interf;
            for (EdgeKey ek : t.edges) {
                Edge& e = edges[ek];
                e.interf = interf;
                vertices[e.v0].interf = interf;
                vertices[e.v1].interf = interf;
            }
        }
    }
}

void SimplexMesh::clearFlags()
{
    std::for_each(vertices.begin(), vertices.end(), [](auto& x) { x.flag = false; });
    std::for_each(edges.begin(), edges.end(), [](auto& x) { x.flag = false; });
    std::for_each(triangles.begin(), triangles.end(), [](auto& x) { x.flag = false; });
    std::for_each(tets.begin(), tets.end(), [](auto& x) { x.flag = false; });
}

void SimplexMesh::saveToDsc(const std::string& fname) const
{
    std::ofstream file(fname);
    assert(file); // Ensure file was opened

    // Write vertices
    for (const Vertex& v : vertices) {
        file << v.pos[0] << ' ' << v.pos[1] << ' ' << v.pos[2] << '\n';
    }

    // Write tetrahedrons
    for (const Tetrahedron& t : tets) {
        file << 't';
        for (VertKey vk : getTetVerts(t.self)) {
            file << ' ' << vk;
        }
        file << ' ' << t.value << '\n';
    }
}

void SimplexMesh::loadFromDsc(const std::string& fname)
{
    std::ifstream file(fname);
    assert(file); // Ensure file was opened
    std::vector<Vec3f> vertexPositions;
    std::vector<std::array<int, 4>> tetIdxs;
    std::vector<int> tetValues;

    while (!file.eof()) {
        char c = '\0';
        file >> c;
        if (c == 'v') {
            float x, y, z;
            file >> x;
            file >> y;
            file >> z;
            vertexPositions.push_back(Vec3f(x, y, z));
        } else if (c == 't') {
            int i1, i2, i3, i4, value;
            file >> i1;
            file >> i2;
            file >> i3;
            file >> i4;
            file >> value;
            tetIdxs.push_back({ i1, i2, i3, i4 });
            tetValues.push_back(value);
        }
    }
    build(vertexPositions, tetIdxs, &tetValues);
}

void SimplexMesh::renderInterface(int insideVal) const
{
    glBegin(GL_TRIANGLES);
    for (const Triangle& t : triangles) {
        if (isTriInterface(t, insideVal)) {
            std::array<VertKey, 3> vs = getTriVerts(t.self);

            const Vec3f p0 = vertices[vs[0]].pos, p1 = vertices[vs[1]].pos, p2 = vertices[vs[2]].pos;

            glNormal3fv(normalize(cross(p1 - p0, p2 - p0)).get());
            glVertex3fv(p0.get());
            glVertex3fv(p1.get());
            glVertex3fv(p2.get());
        }
    }
    glEnd();
}

void SimplexMesh::renderWireframe() const
{
    // Save current polygon mode
    GLint oldMode[2];
    glGetIntegerv(GL_POLYGON_MODE, oldMode);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glBegin(GL_TRIANGLES);
    for (const Triangle& t : triangles) {
        glNormal3fv(Vec3f(0).get()); // Hack to get black lines
        for (VertKey vk : getTriVerts(t.self)) {
            glVertex3fv(vertices[vk].pos.get());
        }
    }
    glEnd();

    // Restore old polygon mode
    glPolygonMode(GL_FRONT, oldMode[0]);
    glPolygonMode(GL_BACK, oldMode[1]);
}

std::array<SimplexMesh::VertKey, 3> SimplexMesh::getTriVerts(TriKey tk) const
{
    const Triangle& t = triangles[tk];
    VertKey v0 = edges[t.edges[0]].v0;
    VertKey v1 = edges[t.edges[0]].v1;
    VertKey v2 = edges[t.edges[1]].v0;
    v2 = (v2 == v1 || v2 == v0) ? edges[t.edges[1]].v1 : v2;
    return { v0, v1, v2 };
}

std::array<SimplexMesh::VertKey, 4> SimplexMesh::getTetVerts(TetKey tk) const
{
    const Tetrahedron& t = tets[tk];
    std::array<VertKey, 4> verts;

    // Add all vertices for the first triangle
    int i = 0;
    for (VertKey vk : getTriVerts(t.tris[0])) {
        verts[i] = vk;
        ++i;
    }

    // Find final vertex in second triangle
    for (VertKey vk : getTriVerts(t.tris[1])) {
        if (vk != verts[0] && vk != verts[1] && vk != verts[2]) {
            // Vertex not yet added - add it now and stop looking
            verts[3] = vk;
            break;
        }
    }

    return verts;
}

std::array<SimplexMesh::EdgeKey, 6> SimplexMesh::getTetEdges(TetKey tk) const
{
    std::array<EdgeKey, 6> edges;
    // Add all vertices for the first triangle
    int i = 0;
    for (EdgeKey ek : triangles[tets[tk].tris[0]].edges) {
        edges[i] = ek;
        ++i;
    }
    // Add remaining vertices from next two triangles
    for (EdgeKey ek : triangles[tets[tk].tris[1]].edges) {
        if (edges[0] != ek && edges[1] != ek && edges[2] != ek) {
            edges[i] = ek;
            ++i;
        }
    }
    for (EdgeKey ek : triangles[tets[tk].tris[2]].edges) {
        if (edges[0] != ek && edges[1] != ek && edges[2] != ek && edges[3] != ek && edges[4] != ek) {
            edges[i] = ek;
            break;
        }
    }
    return edges;
}

std::tuple<std::vector<SimplexMesh::EdgeKey>,
    std::vector<SimplexMesh::TriKey>,
    std::vector<SimplexMesh::TetKey>> SimplexMesh::getVertCoboundary(VertKey vk)
{
    std::vector<EdgeKey> edgeOut;
    std::vector<TriKey> triOut;
    std::vector<TetKey> tetOut;

    std::queue<TetKey> toProcess;
    toProcess.push(triangles[edges[vertices[vk].edge].triangle].t0);
    while (!toProcess.empty()) {
        Tetrahedron& tet = tets[toProcess.front()];
        toProcess.pop();
        tet.flag = true;
        tetOut.push_back(tet.self);

        // Add needed triangles and edges to coboundary
        for (TriKey tk : tet.tris) {
            Triangle& t = triangles[tk];
            if (!t.flag) {
                bool incident = false;
                for (EdgeKey ek : t.edges) {
                    Edge& e = edges[ek];
                    if (e.v0 == vk || e.v1 == vk) {
                        // Vertex vk is a face of this edge, so both the edge and triangle
                        // are in the coboundary
                        incident = true;
                        if (!e.flag) {
                            edgeOut.push_back(ek);
                            e.flag = true;
                        }
                    }
                }
                if (incident) {
                    triOut.push_back(tk);
                    t.flag = true;
                }
            }
        }

        // Find next tets to process, if any
        for (TetKey tk : getAdjTets(tet.self)) {
            if (tk != INVALID_TETRAHEDRON && !tets[tk].flag) {
                auto vs = getTetVerts(tk);
                if (vs[0] == vk || vs[1] == vk || vs[2] == vk || vs[3] == vk) {
                    toProcess.push(tk);
                }
            }
        }
    }
    // Reset touched flags
    for (EdgeKey ek : edgeOut) {
        edges[ek].flag = false;
    }
    for (TriKey tk : triOut) {
        triangles[tk].flag = false;
    }
    for (TetKey tk : tetOut) {
        tets[tk].flag = false;
    }

    return std::make_tuple(edgeOut, triOut, tetOut);
}

std::tuple<std::vector<SimplexMesh::VertKey>,
    std::vector<SimplexMesh::EdgeKey>,
    std::vector<SimplexMesh::TriKey>> SimplexMesh::getVertLink(VertKey vk)
{
    // TODO: This function is similar to getVertCoboundary - find way to unify them
    std::vector<VertKey> vertOut;
    std::vector<EdgeKey> edgeOut;
    std::vector<TriKey> triOut;
    std::vector<TetKey> processedTets;

    std::queue<TetKey> toProcess;
    toProcess.push(triangles[edges[vertices[vk].edge].triangle].t0);
    while (!toProcess.empty()) {
        Tetrahedron& tet = tets[toProcess.front()];
        toProcess.pop();
        tet.flag = true;
        processedTets.push_back(tet.self);

        // Find the triangle face which is not incident with vk and add it- and its
        // sub-simplicies to the link
        for (TriKey tk : tet.tris) {
            if (tk == INVALID_TRIANGLE) {
                continue;
            }
            auto vs = getTriVerts(tk);
            if (vs[0] != vk && vs[1] != vk && vs[2] != vk) {
                // Triangle is not incident with vk so add it to the link
                Triangle& tri = triangles[tk];
                triOut.push_back(tk);
                // We don't need to set the triangle flag, since there is only one 
                // relevant triangle per tet

                // Add triangle edges and vertices to link
                for (EdgeKey ek : tri.edges) {
                    Edge& e = edges[ek];
                    if (!e.flag) {
                        e.flag = true;
                        edgeOut.push_back(ek);
                    }
                }
                for (VertKey vk : vs) {
                    if (!vertices[vk].flag) {
                        vertices[vk].flag = true;
                        vertOut.push_back(vk);
                    }
                }
            }
        }

        // Find next tets to process, if any
        for (TetKey tk : getAdjTets(tet.self)) {
            if (tk != INVALID_TETRAHEDRON && !tets[tk].flag) {
                auto vs = getTetVerts(tk);
                if (vs[0] == vk || vs[1] == vk || vs[2] == vk || vs[3] == vk) {
                    toProcess.push(tk);
                }
            }
        }
    }

    // Reset touched flags
    for (VertKey vk : vertOut) {
        vertices[vk].flag = false;
    }
    for (EdgeKey ek : edgeOut) {
        edges[ek].flag = false;
    }
    for (TriKey tk : triOut) {
        triangles[tk].flag = false;
    }

    return std::make_tuple(vertOut, edgeOut, triOut);
}

std::array<SimplexMesh::TetKey, 4> SimplexMesh::getAdjTets(TetKey tk) const
{
    std::array<TetKey, 4> out;
    int i = 0;
    for (TriKey triKey : tets[tk].tris) {
        if (triangles[triKey].t0 == tk) {
            out[i] = triangles[triKey].t1;
        } else {
            out[i] = triangles[triKey].t0;
        }
        ++i;
    }
    return out;
}

bool operator==(const SimplexMesh& lhs, const SimplexMesh& rhs)
{
    return lhs.vertices == rhs.vertices && lhs.edges == rhs.edges &&
        lhs.triangles == rhs.triangles && lhs.tets == rhs.tets;
}

bool operator!=(const SimplexMesh& lhs, const SimplexMesh& rhs)
{
    return !(lhs == rhs);
}

bool operator==(const SimplexMesh::Vertex& lhs, const SimplexMesh::Vertex& rhs)
{
    return lhs.pos == rhs.pos && lhs.self == rhs.self && lhs.edge == rhs.edge;
}

bool operator!=(const SimplexMesh::Vertex& lhs, const SimplexMesh::Vertex& rhs)
{
    return !(lhs == rhs);
}

bool operator==(const SimplexMesh::Edge& lhs, const SimplexMesh::Edge& rhs)
{
    return lhs.self == rhs.self && lhs.v0 == rhs.v0 && lhs.v1 == rhs.v1 &&
        lhs.triangle == rhs.triangle;
}

bool operator!=(const SimplexMesh::Edge& lhs, const SimplexMesh::Edge& rhs)
{
    return !(lhs == rhs);
}

bool operator==(const SimplexMesh::Triangle& lhs, const SimplexMesh::Triangle& rhs)
{
    return lhs.self == rhs.self && lhs.edges == rhs.edges &&
        lhs.t0 == rhs.t0 && lhs.t1 == rhs.t1;
}

bool operator!=(const SimplexMesh::Triangle& lhs, const SimplexMesh::Triangle& rhs)
{
    return !(lhs == rhs);
}

bool operator==(const SimplexMesh::Tetrahedron& lhs, const SimplexMesh::Tetrahedron& rhs)
{
    return lhs.self == rhs.self && lhs.tris == rhs.tris;
}

bool operator!=(const SimplexMesh::Tetrahedron& lhs, const SimplexMesh::Tetrahedron& rhs)
{
    return !(lhs == rhs);
}
