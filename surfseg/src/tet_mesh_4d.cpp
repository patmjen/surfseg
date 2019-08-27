#include "tet_mesh_4d.h"
#include <algorithm>
#include <fstream>
#include <queue>

#include "util.h"

Vec4f cross4(Vec4f t, Vec4f u, Vec4f v)
{
	// See: https://math.stackexchange.com/a/2371039
	float a1 = t[3] * u[2] * v[1] - t[2] * u[3] * v[1] - t[3] * u[1] * v[2] +
		t[1] * u[3] * v[2] + t[2] * u[1] * v[3] - t[1] * u[2] * v[3];
	float a2 = -t[3] * u[2] * v[0] + t[2] * u[3] * v[0] + t[3] * u[0] * v[2] -
		t[0] * u[3] * v[2] - t[2] * u[0] * v[3] + t[0] * u[2] * v[3];
	float a3 = t[3] * u[1] * v[0] - t[1] * u[3] * v[0] - t[3] * u[0] * v[1] +
		t[0] * u[3] * v[1] + t[1] * u[0] * v[3] - t[0] * u[1] * v[3];
	float a4 = -t[2] * u[1] * v[0] + t[1] * u[2] * v[0] + t[2] * u[0] * v[1] -
		t[0] * u[2] * v[1] - t[1] * u[0] * v[2] + t[0] * u[1] * v[2];
	return Vec4f(a1, a2, a3, a4);
}

TetMesh4d::TetMesh4d(const std::vector<Vec4f>& vertexPositions,
	const std::vector<std::array<int, 4>>& tetIdxs)
{
	build(vertexPositions, tetIdxs);
}

TetMesh4d::TetMesh4d(const TetMesh4d& other) :
	vertices(other.vertices),
	edges(other.edges),
	triangles(other.triangles),
	tets(other.tets)
{}

TetMesh4d::TetMesh4d(TetMesh4d&& other) :
	vertices(std::move(other.vertices)),
	edges(std::move(other.edges)),
	triangles(std::move(other.triangles)),
	tets(std::move(other.tets))
{}

TetMesh4d& TetMesh4d::operator=(const TetMesh4d& other)
{
	if (this != &other) {
		vertices = other.vertices;
		edges = other.edges;
		triangles = other.triangles;
		tets = other.tets;
	}
	return *this;
}

TetMesh4d& TetMesh4d::operator=(TetMesh4d&& other)
{
	if (this != &other) {
		vertices = std::move(other.vertices);
		edges = std::move(other.edges);
		triangles = std::move(other.triangles);
		tets = std::move(other.tets);
	}
	return *this;
}

void TetMesh4d::clear()
{
	vertices.clear();
	edges.clear();
	triangles.clear();
	tets.clear();
}

void TetMesh4d::build(const std::vector<Vec4f>& vertexPositions,
	const std::vector<std::array<int, 4>>& tetIdxs)
{
	clear();
	std::unordered_map<std::pair<int, int>, EdgeKey> edgeMap;
	std::unordered_map<std::array<int, 3>, TriKey> triMap;
	std::unordered_map<std::array<int, 4>, TetKey> tetMap;

	// Add vertices
	vertices.reserve(vertexPositions.size());
	for (const Vec4f& p : vertexPositions) {
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
		++i;
	}
}

TetMesh4d::EdgeKey TetMesh4d::addEdge_(int i1, int i2,
	std::unordered_map<std::pair<int, int>, EdgeKey>& edgeMap)
{
	auto it = edgeMap.find(std::make_pair(i1, i2));
	EdgeKey ek;
	if (it == edgeMap.end()) {
		// Edge has not been added previosly so add it now
		ek = edges.size();
		Edge e;
		e.self = ek;
		e.v0 = i1;
		e.v1 = i2;

		// Register new edge with its vertices - it doesn't matter that we overwrite every time
		vertices[i1].edge = ek;
		vertices[i2].edge = ek;

		edgeMap[std::make_pair(i1, i2)] = ek;
		edgeMap[std::make_pair(i2, i1)] = ek;
		edges.push_back(e);
	} else {
		ek = it->second;
	}
	return ek;
}

TetMesh4d::TriKey TetMesh4d::addTriangle_(int i1, int i2, int i3,
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
			edgeMap.at(std::make_pair(i2, i3)),
			edgeMap.at(std::make_pair(i1, i3)),
			edgeMap.at(std::make_pair(i1, i2))
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

TetMesh4d::TetKey TetMesh4d::addTetrahedron_(int i1, int i2, int i3, int i4,
	std::unordered_map<std::array<int, 4>, TetKey>& tetMap,
	const std::unordered_map<std::array<int, 3>, TriKey>& triMap)
{
	static auto makeTriMapKey = [](int i1, int i2, int i3) {
		std::array<int, 3> key = { i1, i2, i3 };
		std::sort(key.begin(), key.end());
		return key;
	};

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

		// We know all edges for this tet. have been added
		t.tris = {
			triMap.at(makeTriMapKey(i2, i3, i4)),
			triMap.at(makeTriMapKey(i1, i3, i4)),
			triMap.at(makeTriMapKey(i1, i2, i4)),
			triMap.at(makeTriMapKey(i1, i2, i3))
		};

		// We use the original ordering here so the orientation is consistent
		t.verts = { i1, i2, i3, i4 };

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

void TetMesh4d::laplaceSmooth(float lambda, int niter)
{
	if (lambda == 0.0f) {
		// No position changes will occur so return early
		return;
	}
	// Set each vertex position to weighted mean of its neighbors
	std::vector<Vec4f> newPositions;
	newPositions.reserve(vertices.size());
	for (int i = 0; i < niter; ++i) {
		newPositions.clear(); // Won't deallocate memory: https://en.cppreference.com/w/cpp/container/vector/clear
		// Compute all new positions
		for (const auto& v : vertices) {
			std::vector<VertKey> nborVerts;
			std::tie(nborVerts, std::ignore, std::ignore) = getVertLink(v.self);
			Vec4f newPos(0.0f);
			for (auto nk : nborVerts) {
				newPos += vertices[nk].pos;
			}
			newPos /= nborVerts.size();
			newPositions.push_back(newPos);
		}
		// Update positions
		int vi = 0;
		for (auto& v : vertices) {
			v.pos = (1.0f - lambda) * v.pos + lambda * newPositions[vi];
			++vi;
		}
	}
}

Vec4f TetMesh4d::computeCenter() const noexcept
{
	Vec4f center = Vec4f(0);
	for (const Vertex& v : vertices) {
		center += v.pos;
	}
	return center / vertices.size();
}

float TetMesh4d::computeBoundingSphere() const noexcept
{
	Vec4f tmp;
	return computeBoundingSphere(tmp);
}

float TetMesh4d::computeBoundingSphere(Vec4f& center) const noexcept
{
	// Compute centroid
	center = computeCenter();

	// Compute max. distance to centroid
	float r = 0;
	for (const Vertex& v : vertices) {
		r = std::max(r, length(v.pos - center));
	}
	return r;
}

void TetMesh4d::clearFlags()
{
	std::for_each(vertices.begin(), vertices.end(), [](auto& x) { x.flag = false; });
	std::for_each(edges.begin(), edges.end(), [](auto& x) { x.flag = false; });
	std::for_each(triangles.begin(), triangles.end(), [](auto& x) { x.flag = false; });
	std::for_each(tets.begin(), tets.end(), [](auto& x) { x.flag = false; });
}

void TetMesh4d::saveToDsc(const std::string& fname) const
{
	std::ofstream file(fname);
	assert(file); // Ensure file was opened

	// Write vertices
	for (const Vertex& v : vertices) {
		file << v.pos[0] << ' ' << v.pos[1] << ' ' << v.pos[2] << ' ' << v.pos[3] << '\n';
	}

	// Write tetrahedrons
	for (const Tetrahedron& t : tets) {
		file << 't';
		for (VertKey vk : t.verts) {
			file << ' ' << vk;
		}
		file << '\n';
	}
}

void TetMesh4d::loadFromDsc(const std::string& fname)
{
	std::ifstream file(fname);
	assert(file); // Ensure file was opened
	std::vector<Vec4f> vertexPositions;
	std::vector<std::array<int, 4>> tetIdxs;

	while (!file.eof()) {
		char c = '\0';
		file >> c;
		if (c == 'v') {
			float x, y, z, t;
			file >> x;
			file >> y;
			file >> z;
			file >> t;
			vertexPositions.push_back(Vec4f(x, y, z, t));
		} else if (c == 't') {
			int i1, i2, i3, i4, value;
			file >> i1;
			file >> i2;
			file >> i3;
			file >> i4;
			tetIdxs.push_back({ i1, i2, i3, i4 });
		}
	}
	build(vertexPositions, tetIdxs);
}

Vec4f TetMesh4d::tetNormal(TetKey tk) const
{
	return tetNormal(tets[tk]);
}

Vec4f TetMesh4d::tetNormal(const Tetrahedron& tet) const
{
	Vec4f v1 = vertices[tet.verts[0]].pos;
	Vec4f v2 = vertices[tet.verts[1]].pos;
	Vec4f v3 = vertices[tet.verts[2]].pos;
	Vec4f v4 = vertices[tet.verts[3]].pos;

	Vec4f t = v2 - v1;
	Vec4f u = v3 - v1;
	Vec4f v = v4 - v1;

	return cross4(t, u, v);
}

void TetMesh4d::computeTetNormals(bool normalize)
{
	for (auto& tet : tets) {
		Vec4f n = tetNormal(tet);
		if (normalize) {
			n = CGLA::cond_normalize(n);
		}
		tet.normal = n;
	}
}

Vec4f TetMesh4d::vertexNormal(VertKey vk, bool tetNormalsComputed)
{
	std::vector<TetKey> nborTets;
	std::tie(std::ignore, std::ignore, nborTets) = getVertCoboundary(vk);
	Vec4f n(0.0f);
	for (auto tk : nborTets) {
		Vec4f ni = tetNormalsComputed ? tets[tk].normal : tetNormal(tk);
		n += ni;
	}
	return n;
}

Vec4f TetMesh4d::vertexNormal(const Vertex& v, bool tetNormalsComputed)
{
	return vertexNormal(v.self, tetNormalsComputed);
}

void TetMesh4d::computeVertexNormals(bool normalize, bool tetNormalsComputed)
{
	if (!tetNormalsComputed) {
		computeTetNormals();
	}
	for (auto& v : vertices) {
		Vec4f n = vertexNormal(v, true);
		if (normalize) {
			n = CGLA::normalize(n);
		}
		v.normal = n;
	}
}

std::array<TetMesh4d::VertKey, 3> TetMesh4d::getTriVerts(TriKey tk) const
{
	const Triangle& t = triangles[tk];
	VertKey v0 = edges[t.edges[0]].v0;
	VertKey v1 = edges[t.edges[0]].v1;
	VertKey v2 = edges[t.edges[1]].v0;
	v2 = (v2 == v1 || v2 == v0) ? edges[t.edges[1]].v1 : v2;
	return { v0, v1, v2 };
}

std::array<TetMesh4d::EdgeKey, 6> TetMesh4d::getTetEdges(TetKey tk) const
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

std::tuple<std::vector<TetMesh4d::EdgeKey>,
	std::vector<TetMesh4d::TriKey>,
	std::vector<TetMesh4d::TetKey>> TetMesh4d::getVertCoboundary(VertKey vk)
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
				const auto& vs = tets[tk].verts;
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

std::tuple<std::vector<TetMesh4d::VertKey>,
	std::vector<TetMesh4d::EdgeKey>,
	std::vector<TetMesh4d::TriKey>> TetMesh4d::getVertLink(VertKey vk)
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
				const auto& vs = tets[tk].verts;
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

std::array<TetMesh4d::TetKey, 4> TetMesh4d::getAdjTets(TetKey tk) const
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

bool operator==(const TetMesh4d& lhs, const TetMesh4d& rhs)
{
	return lhs.vertices == rhs.vertices && lhs.edges == rhs.edges &&
		lhs.triangles == rhs.triangles && lhs.tets == rhs.tets;
}

bool operator!=(const TetMesh4d& lhs, const TetMesh4d& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const TetMesh4d::Vertex& lhs, const TetMesh4d::Vertex& rhs)
{
	return lhs.pos == rhs.pos && lhs.self == rhs.self && lhs.edge == rhs.edge;
}

bool operator!=(const TetMesh4d::Vertex& lhs, const TetMesh4d::Vertex& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const TetMesh4d::Edge& lhs, const TetMesh4d::Edge& rhs)
{
	return lhs.self == rhs.self && lhs.v0 == rhs.v0 && lhs.v1 == rhs.v1 &&
		lhs.triangle == rhs.triangle;
}

bool operator!=(const TetMesh4d::Edge& lhs, const TetMesh4d::Edge& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const TetMesh4d::Triangle& lhs, const TetMesh4d::Triangle& rhs)
{
	return lhs.self == rhs.self && lhs.edges == rhs.edges &&
		lhs.t0 == rhs.t0 && lhs.t1 == rhs.t1;
}

bool operator!=(const TetMesh4d::Triangle& lhs, const TetMesh4d::Triangle& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const TetMesh4d::Tetrahedron& lhs, const TetMesh4d::Tetrahedron& rhs)
{
	return lhs.self == rhs.self && lhs.tris == rhs.tris;
}

bool operator!=(const TetMesh4d::Tetrahedron& lhs, const TetMesh4d::Tetrahedron& rhs)
{
	return !(lhs == rhs);
}
