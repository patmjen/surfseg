#ifndef TET_MESH_4D_H__
#define TET_MESH_4D_H__

#include <vector>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <string>

#include <GEL/CGLA/Vec4f.h>

using namespace CGLA;

/** Incidence simplicial data structure for storing 3D simplicial complexes in 4D space.
  * NOTE: This code was adapted from the SimplexMesh class.
  */
class TetMesh4d {
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
		bool flag;
		Vec4f pos;
		Vec4f normal;
		VertKey self;
		EdgeKey edge;

		Vertex() :
			flag(false),
			pos(0.0f),
			normal(0.0f),
			self(INVALID_VERTEX),
			edge(INVALID_EDGE) {}
	};

	struct Edge {
		bool flag;
		EdgeKey self;
		VertKey v0;
		VertKey v1;
		TriKey triangle;

		Edge() :
			flag(false),
			self(INVALID_EDGE),
			v0(INVALID_VERTEX),
			v1(INVALID_VERTEX),
			triangle(INVALID_TRIANGLE) {}
	};

	struct Triangle {
		bool flag;
		TriKey self;
		std::array<EdgeKey, 3> edges;
		TetKey t0;
		TetKey t1;

		Triangle() :
			flag(false),
			self(INVALID_TRIANGLE),
			edges({ INVALID_EDGE, INVALID_EDGE, INVALID_EDGE }),
			t0(INVALID_TETRAHEDRON),
			t1(INVALID_TETRAHEDRON) {}
	};

	struct Tetrahedron {
		bool flag;
		TetKey self;
		std::array<TriKey, 4> tris;
		std::array<VertKey, 4> verts;
		Vec4f normal;

		Tetrahedron() :
			flag(false),
			self(INVALID_TETRAHEDRON),
			tris({ INVALID_TRIANGLE, INVALID_TRIANGLE, INVALID_TRIANGLE, INVALID_TRIANGLE }),
			verts({ INVALID_VERTEX, INVALID_VERTEX, INVALID_VERTEX, INVALID_VERTEX }),
			normal(0.0f) {}
	};

	// TODO: Store flags in seperate arrays so we can have proper const qualifiers
	// TODO: Store positions and normals in seperate vectors like ManifoldMesh
	// TODO: Profile speed-up of above suggestion
	std::vector<Vertex> vertices;
	std::vector<Edge> edges;
	std::vector<Triangle> triangles;
	std::vector<Tetrahedron> tets;

	explicit TetMesh4d() = default;
	explicit TetMesh4d(const std::vector<Vec4f>& vertexPositions,
		const std::vector<std::array<int, 4>>& tetIdxs);
	TetMesh4d(const TetMesh4d& other);
	TetMesh4d(TetMesh4d&& other);

	TetMesh4d& operator=(const TetMesh4d& other);
	TetMesh4d& operator=(TetMesh4d&& other);

	void clear();
	void build(const std::vector<Vec4f>& vertexPositions,
		const std::vector<std::array<int, 4>>& tetIdxs);

	void laplaceSmooth(float lambda, int niter = 1);

	Vec4f computeCenter() const noexcept;
	float computeBoundingSphere() const noexcept;
	float computeBoundingSphere(Vec4f& center) const noexcept;

	void clearFlags();

	void saveToDsc(const std::string& fname) const;
	void loadFromDsc(const std::string& fname);

	Vec4f tetNormal(TetKey tk) const;
	Vec4f tetNormal(const Tetrahedron& t) const;
	void computeTetNormals(bool normalize = true);

	Vec4f vertexNormal(VertKey vk, bool tetNormalsComputed = false);
	Vec4f vertexNormal(const Vertex& v, bool tetNormalsComputed = false);
	void computeVertexNormals(bool normalize = true, bool tetNormalsComputed = false);

	std::array<VertKey, 3> getTriVerts(TriKey tk) const;

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

bool operator==(const TetMesh4d& lhs, const TetMesh4d& rhs);
bool operator!=(const TetMesh4d& lhs, const TetMesh4d& rhs);

bool operator==(const TetMesh4d::Vertex& lhs, const TetMesh4d::Vertex& rhs);
bool operator!=(const TetMesh4d::Vertex& lhs, const TetMesh4d::Vertex& rhs);

bool operator==(const TetMesh4d::Edge& lhs, const TetMesh4d::Edge& rhs);
bool operator!=(const TetMesh4d::Edge& lhs, const TetMesh4d::Edge& rhs);

bool operator==(const TetMesh4d::Triangle& lhs, const TetMesh4d::Triangle& rhs);
bool operator!=(const TetMesh4d::Triangle& lhs, const TetMesh4d::Triangle& rhs);

bool operator==(const TetMesh4d::Tetrahedron& lhs, const TetMesh4d::Tetrahedron& rhs);
bool operator!=(const TetMesh4d::Tetrahedron& lhs, const TetMesh4d::Tetrahedron& rhs);

#endif // TET_MESH_4D_H__
