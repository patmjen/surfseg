#ifndef SIMPLEX_MESH_H__
#define SIMPLEX_MESH_H__

#include <vector>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <string>

#include <GEL/CGLA/Vec3f.h>

using namespace CGLA;

/** Incidence simplicial data structure for storing 3D simplicial complexes */
class SimplexMesh {
public:
	typedef int VertKey;
	typedef int EdgeKey;
	typedef int TriKey;
	typedef int TetKey;

	static const VertKey INVALID_VERTEX = -1;
	static const EdgeKey INVALID_EDGE = -1;
	static const TriKey INVALID_TRIANGLE = -1;
	static const TetKey INVALID_TETRAHEDRON = -1;

	struct Vertex {
		bool interf;
		bool flag;
		Vec3f pos;
		VertKey self;
		EdgeKey edge;

		Vertex() :
			interf(false),
			flag(false),
			pos(0.0f),
			self(INVALID_VERTEX),
			edge(INVALID_EDGE) {}
	};

	struct Edge {
		bool interf;
		bool flag;
		EdgeKey self;
		VertKey v0;
		VertKey v1;
		TriKey triangle;

		Edge() :
			interf(false),
			flag(false),
			self(INVALID_EDGE),
			v0(INVALID_VERTEX),
			v1(INVALID_VERTEX),
			triangle(INVALID_TRIANGLE) {}
	};

	struct Triangle {
		bool interf;
		bool flag;
		TriKey self;
		std::array<EdgeKey, 3> edges;
		TetKey t0;
		TetKey t1;

		Triangle() :
			interf(false),
			flag(false),
			self(INVALID_TRIANGLE),
			edges({ INVALID_EDGE, INVALID_EDGE, INVALID_EDGE }),
			t0(INVALID_TETRAHEDRON),
			t1(INVALID_TETRAHEDRON) {}
	};

	struct Tetrahedron {
		int value;
		bool flag;
		TetKey self;
		std::array<TriKey, 4> tris;

		Tetrahedron() :
			value(0),
			flag(false),
			self(INVALID_TETRAHEDRON),
			tris({ INVALID_TRIANGLE, INVALID_TRIANGLE, INVALID_TRIANGLE, INVALID_TRIANGLE }) {}
	};

	std::vector<Vertex> vertices;
	std::vector<Edge> edges;
	std::vector<Triangle> triangles;
	std::vector<Tetrahedron> tets;

	explicit SimplexMesh() = default;
	explicit SimplexMesh(const SimplexMesh& other);
	explicit SimplexMesh(SimplexMesh&& other);

	SimplexMesh& operator=(const SimplexMesh& other);
	SimplexMesh& operator=(SimplexMesh&& other);

	void clear();
	void build(const std::vector<Vec3f>& vertexPositions,
		const std::vector<std::array<int, 4>>& tetIdxs,
		const std::vector<int> *tetValues = nullptr);

	float computeBoundingSphere() const noexcept;
	float computeBoundingSphere(Vec3f& center) const noexcept;

	bool isTriInterface(TriKey tk, int insideVal = 1) const;
	bool isTriInterface(const Triangle& t, int insideVal = 1) const;
	void setInterfaceFlags();

	void clearFlags();

	void saveToDsc(const std::string& fname) const;
	void loadFromDsc(const std::string& fname);

	void renderInterface(int insideVal = 1) const;
	void renderWireframe() const;

	std::array<VertKey, 3> getTriVerts(TriKey tk) const;
	std::array<VertKey, 4> getTetVerts(TetKey tk) const;

	std::array<EdgeKey, 6> getTetEdges(TetKey tk) const;

	std::tuple<std::vector<EdgeKey>, std::vector<TriKey>, std::vector<TetKey>>
		getVertCoboundary(VertKey vk);
	std::tuple<std::vector<VertKey>, std::vector<EdgeKey>, std::vector<TriKey>>
		getVertLink(VertKey vk);

	std::array<TetKey, 4> getAdjTets(TetKey tk) const;

private:
	EdgeKey addEdge_(int i1, int i2,
		std::unordered_map<std::pair<int, int>, EdgeKey>& edgeMap);

	TriKey addTriangle_(int i1, int i2, int i3,
		std::unordered_map<std::array<int, 3>, TriKey>& triMap,
		const std::unordered_map<std::pair<int, int>, EdgeKey>& edgeMap);

	TetKey addTetrahedron_(int i1, int i2, int i3, int i4,
		std::unordered_map<std::array<int, 4>, TetKey>& tetMap,
		const std::unordered_map<std::array<int, 3>, TriKey>& triMap);
};

bool operator==(const SimplexMesh& lhs, const SimplexMesh& rhs);
bool operator!=(const SimplexMesh& lhs, const SimplexMesh& rhs);

bool operator==(const SimplexMesh::Vertex& lhs, const SimplexMesh::Vertex& rhs);
bool operator!=(const SimplexMesh::Vertex& lhs, const SimplexMesh::Vertex& rhs);

bool operator==(const SimplexMesh::Edge& lhs, const SimplexMesh::Edge& rhs);
bool operator!=(const SimplexMesh::Edge& lhs, const SimplexMesh::Edge& rhs);

bool operator==(const SimplexMesh::Triangle& lhs, const SimplexMesh::Triangle& rhs);
bool operator!=(const SimplexMesh::Triangle& lhs, const SimplexMesh::Triangle& rhs);

bool operator==(const SimplexMesh::Tetrahedron& lhs, const SimplexMesh::Tetrahedron& rhs);
bool operator!=(const SimplexMesh::Tetrahedron& lhs, const SimplexMesh::Tetrahedron& rhs);

#endif // SIMPLEX_MESH_H__
