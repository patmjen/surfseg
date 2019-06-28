#include <unordered_map>
#include <utility>
#include <fstream>
#include <cassert>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <cmath>

#include <GL/glew.h>

#include "manifold_mesh.h"
#include "util.h"

const ManifoldMesh::VertKey ManifoldMesh::INVALID_VERT;
const ManifoldMesh::FaceKey ManifoldMesh::INVALID_FACE;
const ManifoldMesh::EdgeKey ManifoldMesh::INVALID_EDGE;

ManifoldMesh::ManifoldMesh(const ManifoldMesh& other) :
	vertices(other.vertices),
	faces(other.faces),
	edges(other.edges),
	vertPositions(other.vertPositions),
	vertNormals(other.vertNormals),
	faceNormals(other.faceNormals)
{}

ManifoldMesh::ManifoldMesh(ManifoldMesh&& other) :
	vertices(std::move(other.vertices)),
	faces(std::move(other.faces)),
	edges(std::move(other.edges)),
	vertPositions(std::move(other.vertPositions)),
	vertNormals(std::move(other.vertNormals)),
	faceNormals(std::move(other.faceNormals))
{}

ManifoldMesh& ManifoldMesh::operator=(const ManifoldMesh& other)
{
	if (this != &other) {
		vertices = other.vertices;
		faces = other.faces;
		edges = other.edges;
		vertPositions = other.vertPositions;
		vertNormals = other.vertNormals;
		faceNormals = other.faceNormals;
	}
	return *this;
}

ManifoldMesh& ManifoldMesh::operator=(ManifoldMesh&& other)
{
	if (this != &other) {
		vertices = std::move(other.vertices);
		faces = std::move(other.faces);
		edges = std::move(other.edges);
		vertPositions = std::move(other.vertPositions);
		vertNormals = std::move(other.vertNormals);
		faceNormals = std::move(other.faceNormals);
	}
	return *this;
}

void ManifoldMesh::clear()
{
	vertices.clear();
	faces.clear();
	edges.clear();
	vertPositions.clear();
	vertNormals.clear();
	faceNormals.clear();
}

void ManifoldMesh::build(std::vector<Vec3f> vertexPositions,
	const std::vector<std::vector<int>>& faceIdxLists, size_t expectedEdgeCount)
{
	typedef std::pair<int, int> HEKey; // Helper typedef

	// Declare these now so we can capture them in the lambda
	std::unordered_map<HEKey, EdgeKey> edgeMap;
	if (expectedEdgeCount > 0) {
		edgeMap.reserve(expectedEdgeCount);
	} else {
		// Since the bottleneck is adding stuff to edgeMap, we first compute the edge count
		size_t numEdges = std::accumulate(faceIdxLists.begin(), faceIdxLists.end(), 0,
			[](auto count, auto faceIdxs) { return count + faceIdxs.size(); });
		edgeMap.reserve(2 * numEdges);
	}
	EdgeKey prevEdge;
	FaceKey crntFace;

	// Helper function
	auto addEdge = [&](int i1, int i2, const std::vector<int>& faceIdxs) {
		int f1 = faceIdxs[i1], f2 = faceIdxs[i2];
		HEKey ekCrnt = std::make_pair(f1, f2);
		HEKey ekTwin = std::make_pair(f2, f1);

		// If the current edge is already registered the mesh is malformed
		assert(edgeMap.find(ekCrnt) == edgeMap.end());

		// Add new half-edge
		EdgeKey crntEdge = edges.size();
		edges.push_back(HalfEdge());
		edges[crntEdge].self = crntEdge;
		edges[crntEdge].prev = prevEdge;
		edges[crntEdge].vert = f1; // We know vertex indices corresponds to their keys
		edges[crntEdge].face = crntFace;
		if (prevEdge != INVALID_EDGE) {
			edges[prevEdge].next = crntEdge;
		}

		if (vertices[f1].edge == INVALID_EDGE) {
			vertices[f1].edge = crntEdge;
		}

		// Register edge
		edgeMap[ekCrnt] = crntEdge;

		auto twinEdgeRes = edgeMap.find(ekTwin);
		if (twinEdgeRes != edgeMap.end()) {
			// Twin is already registered so update fields
			EdgeKey twinEdge = twinEdgeRes->second;
			edges[twinEdge].twin = crntEdge;
			edges[crntEdge].twin = twinEdge;
		}

		return crntEdge;
	};

	// Add vertices
	vertPositions = vertexPositions;
	for (const Vec3f& pos : vertexPositions) {
		Vertex vert;
		vert.self = vertices.size();
		vertices.push_back(vert);
	}

	faces.reserve(faceIdxLists.size());
	for (const auto& faceIdxs : faceIdxLists) {
		assert(faceIdxs.size() >= 3);  // A face must have at least 3 vertices

		crntFace = faces.size();
		faces.push_back(Face());

		prevEdge = INVALID_EDGE;
		EdgeKey firstEdge = INVALID_EDGE;
		int i1, i2;
		for (i1 = 0, i2 = 1; i2 < faceIdxs.size(); i1++, i2++) {
			EdgeKey crntEdge = addEdge(i1, i2, faceIdxs);
			prevEdge = crntEdge;
			firstEdge = (firstEdge != INVALID_EDGE) ? firstEdge : crntEdge;
		}
		i2 = 0;
		EdgeKey crntEdge = addEdge(i1, i2, faceIdxs);
		edges[crntEdge].next = firstEdge;
		edges[firstEdge].prev = crntEdge;

		faces.back().edge = firstEdge;
		faces.back().self = crntFace;
	}
}

void ManifoldMesh::laplaceSmooth(float lambda, int niter)
{
	if (lambda == 0.0f) {
		// Nothing will happen for this lambda value, so just return now
		return;
	}
	// Set each vertex position to weighted mean of its neighbors
	for (int iter = 0; iter < niter; ++iter) {
		std::vector<Vec3f> newPositions;
		newPositions.reserve(vertices.size());
		// Compute all new positions
		for (const auto& v : vertices) {
			EdgeKey ek0 = twin(v.edge);
			EdgeKey ek = ek0;
			Vec3f newPos(0);
			float nn = 0.0f;
			do {
				const Vec3f p = vertPositions[edges[ek].vert];
				newPos += p;
				nn += 1.0f;
				ek = twin(next(ek)); // Next neighbor
			} while (ek != ek0);
			newPositions.push_back((1.0f - lambda) * vpos(v) + lambda * (newPos / nn));
		}
		// Update vertex positions
		vertPositions = std::move(newPositions);
	}
}

void ManifoldMesh::taubinSmooth(float lambda, float mu, int niter)
{
	for (int i = 0; i < niter; ++i) {
		laplaceSmooth(lambda);
		laplaceSmooth(-mu);
	}
}

void ManifoldMesh::splitEdge(EdgeKey ek, float a)
{
	// Get relevant parts of mesh
	EdgeKey ekt = twin(ek);
	HalfEdge& e = edges[ek];
	HalfEdge& et = edges[ekt];

	FaceKey fk1 = e.face;
	FaceKey fk2 = et.face;

	VertKey vk0 = e.vert;
	VertKey vk1 = et.vert;
	VertKey vk2 = edges[e.prev].vert;
	VertKey vk3 = edges[et.prev].vert;

	// Make new elements
	EdgeKey ek41 = edges.size() + 0;
	EdgeKey ek40 = edges.size() + 1;
	EdgeKey ek42 = edges.size() + 2;
	EdgeKey ek43 = edges.size() + 3;
	EdgeKey ek24 = edges.size() + 4;
	EdgeKey ek34 = edges.size() + 5;

	VertKey vkn = vertices.size();

	FaceKey fk3 = faces.size() + 0;
	FaceKey fk4 = faces.size() + 1;

	// Allocate extra memory
	// *** NOTE: From this point we can no longer use e and et! ***
	edges.resize(edges.size() + 6);
	vertices.resize(vertices.size() + 1);
	faces.resize(faces.size() + 2);

	// Reconnect
	std::array<EdgeKey, 4> inOut = { ek40, ek43, ek41, ek42 };
	std::array<EdgeKey, 4> outIn = { ek, ek34, ekt, ek24 };
	std::array<EdgeKey, 4> outer = { prev(ek), next(ekt), prev(ekt), next(ek) };
	std::array<FaceKey, 4> newFaces = { fk1, fk4, fk2, fk3 };
	std::array<VertKey, 4> outerVerts = { vk0, vk3, vk1, vk2 };
	for (int i = 0; i < 4; i++) {
		HalfEdge& o = edges[outer[i]];
		o.prev = inOut[(i + 3) % 4];
		o.next = outIn[i];
		o.face = newFaces[i];

		HalfEdge& oi = edges[outIn[i]];
		oi.self = outIn[i];
		oi.prev = outer[i];
		oi.next = inOut[(i + 3) % 4];
		oi.twin = inOut[i];
		oi.vert = outerVerts[i];
		oi.face = newFaces[i];

		HalfEdge& io = edges[inOut[i]];
		io.self = inOut[i];
		io.prev = outIn[(i + 1) % 4];
		io.next = outer[(i + 1) % 4];
		io.twin = outIn[i];
		io.vert = vkn;
		io.face = newFaces[(i + 1) % 4];

		Face& f = faces[newFaces[i]];
		f.self = newFaces[i];
		f.edge = outer[i];
	}
	Vertex& vn = vertices[vkn];
	vn.self = vkn;
	vn.edge = ek40;
	vpos(vn) = 0.25*(vpos(vk0) + vpos(vk1) + vpos(vk2) + vpos(vk3));
}

void ManifoldMesh::flipEdge(EdgeKey ek)
{
	// Get relevant parts of the mesh
	HalfEdge& e = edges[ek];
	HalfEdge& et = edges[twin(ek)];
	Vertex& v0 = vertices[e.vert];
	Vertex& v1 = vertices[et.vert];
	Vertex& vn0 = vertices[edges[prev(ek)].vert];
	Vertex& vn1 = vertices[edges[prev(twin(ek))].vert];

	Face& f1 = faces[e.face];
	Face& f2 = faces[et.face];

	// Disconnect vertices and edges
	v0.edge = et.next;
	v1.edge = e.next;

	// Reconnect
	edges[e.prev].next = et.next;
	edges[e.prev].prev = et.self;

	edges[e.next].next = e.self;
	edges[e.next].prev = et.prev;

	edges[et.prev].next = e.next;
	edges[et.prev].prev = e.self;

	edges[et.next].next = et.self;
	edges[et.next].prev = e.prev;

	e.vert = vn0.self;
	et.vert = vn1.self;
	std::swap(e.prev, e.next); // Now e.prev is correct
	std::swap(et.prev, et.next); // Now, et.prev is correct
	std::swap(e.next, et.next);

	// Update face data
	f1.edge = e.self;
	edges[e.next].face = f1.self;
	edges[e.prev].face = f1.self;

	f2.edge = et.self;
	edges[et.next].face = f2.self;
	edges[et.prev].face = f2.self;
}

void ManifoldMesh::removeEdge(EdgeKey ek)
{
	HalfEdge& e = edges[ek];
	HalfEdge& et = edges[e.twin];

	// Reconnect
	vertices[e.vert].edge = twin(prev(ek));
	vertices[edges[e.prev].vert].edge = twin(e.next);
	vertices[edges[et.prev].vert].edge = twin(et.next);

	edges[twin(next(ek))].twin = twin(prev(ek));
	edges[twin(prev(ek))].twin = twin(next(ek));
	edges[twin(next(twin(ek)))].twin = twin(prev(twin(ek)));
	edges[twin(prev(twin(ek)))].twin = twin(next(twin(ek)));
	EdgeKey ek0 = twin(prev(twin(ek)));
	EdgeKey eki = ek0;
	do {
		edges[eki].vert = e.vert;
		eki = next(twin(eki));
	} while (eki != ek0); // Next neighbor

	// Mark elements for deletion
	edges[prev(e.self)].self = INVALID_EDGE;
	edges[next(e.self)].self = INVALID_EDGE;
	e.self = INVALID_EDGE;

	edges[prev(et.self)].self = INVALID_EDGE;
	edges[next(et.self)].self = INVALID_EDGE;
	et.self = INVALID_EDGE;

	vertices[et.vert].self = INVALID_VERT;
	faces[e.face].self = INVALID_FACE;
	faces[et.face].self = INVALID_FACE;
}

void ManifoldMesh::removeInvalid()
{
	// Remove all elements with an invalid self key
	int i = 0, kNew = 0;
	std::vector<VertKey> vertMap(vertices.size(), INVALID_VERT);
	std::vector<FaceKey> faceMap(faces.size(), INVALID_FACE);
	std::vector<EdgeKey> edgeMap(edges.size(), INVALID_EDGE);
	vertices.erase(std::remove_if(vertices.begin(), vertices.end(), [&](const auto& v) {
		bool del = v.self == INVALID_VERT;
		if (!del) {
			vertMap[i] = kNew;
			kNew++;
		}
		i++;
		return del;
	}), vertices.end());
	i = 0;
	kNew = 0;
	faces.erase(std::remove_if(faces.begin(), faces.end(), [&](const auto& f) {
		bool del = f.self == INVALID_FACE;
		if (!del) {
			faceMap[i] = kNew;
			kNew++;
		}
		i++;
		return del;
	}), faces.end());
	i = 0;
	kNew = 0;
	edges.erase(std::remove_if(edges.begin(), edges.end(), [&](const auto& e) {
		bool del = e.self == INVALID_EDGE;
		if (!del) {
			edgeMap[i] = kNew;
			kNew++;
		}
		i++;
		return del;
	}), edges.end());
	// Update all keys
	for (Vertex& v : vertices) {
		v.self = vertMap[v.self];
		v.edge = edgeMap[v.edge];
	}
	for (Face& f : faces) {
		f.self = faceMap[f.self];
		f.edge = edgeMap[f.edge];
	}
	for (HalfEdge& e : edges) {
		e.self = edgeMap[e.self];
		e.prev = edgeMap[e.prev];
		e.next = edgeMap[e.next];
		e.twin = edgeMap[e.twin];
		e.vert = vertMap[e.vert];
		e.face = faceMap[e.face];
	}
}

float ManifoldMesh::computeBoundingSphere() const noexcept
{
	Vec3f tmp;
	return computeBoundingSphere(tmp);
}

float ManifoldMesh::computeBoundingSphere(Vec3f& center) const noexcept
{
	// Compute centroid
	center = Vec3f(0);
	for (const auto& p : vertPositions) {
		center += p;
	}
	center /= vertices.size();

	// Compute max. distance to centroid
	float r = 0;
	for (const auto& p : vertPositions) {
		r = std::max(r, length(p - center));
	}
	return r;
}

std::vector<float> ManifoldMesh::gaussCurvatures() const
{
	std::vector<float> gaussCurvs;
	gaussCurvs.reserve(vertices.size());

	for (const auto& v : vertices) {
		float Ag = 2.0f * M_PI;
		float As = 0.0f;
		EdgeKey ek0 = twin(v.edge);
		EdgeKey ek = ek0;
		do {
			EdgeKey ekn = twin(next(ek));  // Next neighbor
			const Vec3f& e1 = vpos(edges[ek].vert) - vpos(v);
			const Vec3f& e2 = vpos(edges[ekn].vert) - vpos(v);
			Ag -= acosf(dot(normalize(e1), normalize(e2)));
			As += length(cross(e1, e2));
			ek = ekn;
		} while (ek != ek0);
		As = As / 6.0f; // Divide by 2 and 3
		gaussCurvs.push_back(Ag / As);
	}
	return gaussCurvs;
}

Vec3f ManifoldMesh::faceNormal(FaceKey fk) const
{
	return faceNormal(faces[fk]);
}

Vec3f ManifoldMesh::faceNormal(const Face& f) const
{
	EdgeKey ek = f.edge;
	Vec3f v0 = vpos(edges[ek].vert);
	Vec3f v1 = vpos(edges[next(ek)].vert);
	Vec3f v2 = vpos(edges[next(next(ek))].vert);
	return cross(v1 - v0, v2 - v0);
}

void ManifoldMesh::computeFaceNormals(bool normalize)
{
	faceNormals.clear();
	faceNormals.reserve(faces.size());
	for (Face& f : faces) {
		Vec3f n = faceNormal(f);
		if (normalize) {
			n = CGLA::normalize(n);
		}
		faceNormals.push_back(n);
	}
}

Vec3f ManifoldMesh::vertexNormal(VertKey vk, bool faceNormalsComputed) const
{
	return vertexNormal(vertices[vk], faceNormalsComputed);
}

Vec3f ManifoldMesh::vertexNormal(const Vertex& v, bool faceNormalsComputed) const
{
	Vec3f n(0.0f);
	EdgeKey ek = v.edge;
	do {
		EdgeKey ekn = next(twin(ek));
		FaceKey fk = edges[ek].face;
		Vec3f e1 = vpos(edges[twin(ek)].vert) - vpos(v);
		Vec3f e2 = vpos(edges[twin(ekn)].vert) - vpos(v);
		float a = acosf(dot(normalize(e1), normalize(e2)));
		n += a * (faceNormalsComputed ? faceNormals[fk] : faceNormal(fk));
		ek = ekn;
	} while (ek != v.edge);
	return n;
}

void ManifoldMesh::computeVertexNormals(bool faceNormalsComputed, bool normalize)
{
	if (!faceNormalsComputed) {
		computeFaceNormals(normalize);
	}
	vertNormals.clear();
	vertNormals.reserve(vertices.size());
	for (Vertex& v : vertices) {
		Vec3f n = vertexNormal(v, true);
		if (normalize) {
			n = CGLA::normalize(n);
		}
		vertNormals.push_back(n);
	}
}

void ManifoldMesh::saveToObj(const std::string& fname) const
{
	std::ofstream file(fname);
	assert(file); // Ensure file was opened

	// Write vertices
	for (const auto& p : vertPositions) {
		file << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
	}

	// Write faces
	for (const auto& f : faces) {
		EdgeKey firstEdge = f.edge;
		file << "f " << edges[firstEdge].vert + 1;
		for (EdgeKey edge = next(firstEdge); edge != firstEdge; edge = next(edge)) {
			file << " " << edges[edge].vert + 1;
		}
		file << "\n";
	}
}

void ManifoldMesh::loadFromObj(const std::string& fname)
{
	std::ifstream file(fname);
	assert(file); // Ensure file was opened
	std::vector<Vec3f> vertexPositions;
	std::vector<std::vector<int>> faceIdxLists;

	// TODO: This is very slow, replace with faster method
	size_t numEdges = 0;
	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::istream_iterator<std::string> iter(iss);
		const std::istream_iterator<std::string> iterend;
		if (*iter == "#") {
			// Comment line so skip
			continue;
		}
		else if (*iter == "v") {
			// Vertex
			float x = std::stof(*(++iter));
			float y = std::stof(*(++iter));
			float z = std::stof(*(++iter));
			vertexPositions.push_back(Vec3f(x, y, z));
		}
		else if (*iter == "f") {
			// Face
			std::vector<int> faceIdxs;
			for (++iter; iter != iterend; ++iter) {
				faceIdxs.push_back(stoi(*iter) - 1);
				numEdges++;
			}
			faceIdxLists.push_back(faceIdxs);
		}
	}

	clear();
	build(vertexPositions, faceIdxLists, numEdges);
}

#ifndef NO_OPENGL
void ManifoldMesh::render(bool vertexNormals) const
{
	glBegin(GL_TRIANGLES);
	for (const Face& f : faces) {
		EdgeKey e0 = f.edge;
		EdgeKey e1 = next(e0);
		EdgeKey e2 = next(e1);

		Vertex v0 = vertices[edges[e0].vert];
		while (e2 != e0) {
			Vertex v1 = vertices[edges[e1].vert];
			Vertex v2 = vertices[edges[e2].vert];

			if (vertexNormals) {
				glNormal3fv(vnormal(v0).get());
				glVertex3fv(vpos(v0).get());
				glNormal3fv(vnormal(v1).get());
				glVertex3fv(vpos(v1).get());
				glNormal3fv(vnormal(v2).get());
				glVertex3fv(vpos(v2).get());
			}
			else {
				glNormal3fv(fnormal(f).get());
				glVertex3fv(vpos(v0).get());
				glVertex3fv(vpos(v1).get());
				glVertex3fv(vpos(v2).get());
			}

			e1 = e2;
			e2 = next(e2);
		}
	}
	glEnd();
}
#endif

ManifoldMesh::EdgeKey ManifoldMesh::next(EdgeKey ek) const
{
	return edges[ek].next;
}

ManifoldMesh::EdgeKey ManifoldMesh::prev(EdgeKey ek) const
{
	return edges[ek].prev;
}

ManifoldMesh::EdgeKey ManifoldMesh::twin(EdgeKey ek) const
{
	return edges[ek].twin;
}

Vec3f& ManifoldMesh::vpos(const Vertex& v)
{
	return vpos(v.self);
}

const Vec3f& ManifoldMesh::vpos(const Vertex& v) const
{
	return vpos(v.self);
}

Vec3f& ManifoldMesh::vpos(VertKey vk)
{
	return vertPositions[vk];
}

const Vec3f& ManifoldMesh::vpos(VertKey vk) const
{
	return vertPositions[vk];
}

Vec3f& ManifoldMesh::vnormal(const Vertex& v)
{
	return vnormal(v.self);
}

const Vec3f& ManifoldMesh::vnormal(const Vertex& v) const
{
	return vnormal(v.self);
}

Vec3f& ManifoldMesh::vnormal(VertKey vk)
{
	return vertNormals[vk];
}

const Vec3f& ManifoldMesh::vnormal(VertKey vk) const
{
	return vertNormals[vk];
}

Vec3f& ManifoldMesh::fnormal(const Face& f)
{
	return fnormal(f.self);
}

const Vec3f& ManifoldMesh::fnormal(const Face& f) const
{
	return fnormal(f.self);
}

Vec3f& ManifoldMesh::fnormal(FaceKey fk)
{
	return faceNormals[fk];
}

const Vec3f& ManifoldMesh::fnormal(FaceKey fk) const
{
	return faceNormals[fk];
}

void ManifoldMesh::assertConsistent() const
{
#ifndef NDEBUG
	for (const auto& e : edges) {
		if (e.self == INVALID_EDGE) {
			continue;
		}
		assert(e.prev != INVALID_EDGE && edges[e.prev].self != INVALID_EDGE);
		assert(e.next != INVALID_EDGE && edges[e.next].self != INVALID_EDGE);
		assert(e.vert != INVALID_VERT && vertices[e.vert].self != INVALID_VERT);
		assert(e.face != INVALID_FACE && faces[e.face].self != INVALID_FACE);
		assert(next(next(next(e.self))) == e.self);
		assert(prev(prev(prev(e.self))) == e.self);
		if (e.twin != INVALID_EDGE) {
			assert(edges[e.twin].twin == e.self);
		}
	}
	for (const auto& v : vertices) {
		if (v.self == INVALID_VERT) {
			continue;
		}
		assert(v.edge != INVALID_EDGE && edges[v.edge].self != INVALID_EDGE);
		assert(edges[v.edge].vert == v.self);
	}
	for (const auto& f : faces) {
		if (f.self == INVALID_FACE) {
			continue;
		}
		assert(f.edge != INVALID_EDGE && edges[f.edge].self != INVALID_EDGE);
		assert(edges[f.edge].face == f.self);
	}
#endif
}

bool operator==(const ManifoldMesh& lhs, const ManifoldMesh& rhs)
{
	return lhs.vertices == rhs.vertices && lhs.faces == rhs.faces && lhs.edges == rhs.edges;
}

bool operator!=(const ManifoldMesh& lhs, const ManifoldMesh& rhs)
{
	return !(lhs == rhs);
}

bool operator==(const ManifoldMesh::Vertex& lhs, const ManifoldMesh::Vertex& rhs) {
	// TODO: Use more robust float comparison
	return lhs.self == rhs.self && lhs.edge == rhs.edge;
}

bool operator!=(const ManifoldMesh::Vertex& lhs, const ManifoldMesh::Vertex& rhs) {
	return !(lhs == rhs);
}

bool operator==(const ManifoldMesh::Face& lhs, const ManifoldMesh::Face& rhs) {
	// TODO: Use more robust float comparison
	return lhs.self == rhs.self && lhs.edge == rhs.edge;
}

bool operator!=(const ManifoldMesh::Face& lhs, const ManifoldMesh::Face& rhs) {
	return !(lhs == rhs);
}

bool operator==(const ManifoldMesh::HalfEdge& lhs, const ManifoldMesh::HalfEdge& rhs) {
	return lhs.self == rhs.self && lhs.prev == rhs.prev && lhs.next == rhs.next &&
		lhs.twin == rhs.twin && lhs.vert == rhs.vert && lhs.face == rhs.face;
}

bool operator!=(const ManifoldMesh::HalfEdge& lhs, const ManifoldMesh::HalfEdge& rhs) {
	return !(lhs == rhs);
}
