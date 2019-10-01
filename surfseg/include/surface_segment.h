#ifndef SURFACE_SEGMENT_H__
#define SURFACE_SEGMENT_H__

#include <vector>
#include <unordered_set>

#include <GEL/CGLA/Vec3f.h>

#include "manifold_mesh.h"
#include "subdivided_icosahedron.h"
#include "tet_mesh_4d.h"
#include "volume.h"
#include "volume4d.h"
#include "graph.h"

using namespace CGLA;

typedef Graph<float, float, float> FloatGraph;

enum CostType : int {
	ON_SURFACE = 0,
	IN_SURFACE = 1
};

ManifoldMesh surfaceCut(const Volume<float>& cost, const ManifoldMesh& init,
    int numSamples, float sampleStep, int maxDiff, CostType costType);
ManifoldMesh surfaceCut(const Volume<float>& cost, const ManifoldMesh&& init,
    int numSamples, float sampleStep, int maxDiff, CostType costType);

TetMesh4d surfaceCut4d(const Volume4d<float>& cost, TetMesh4d mesh,
	int numSamples, float sampleStep, int maxDiff, CostType costType,
	const std::unordered_set<int>& frozenVerts = std::unordered_set<int>());
template <class Func>
TetMesh4d surfaceCut4d(const Volume4d<float>& vol, TetMesh4d mesh,
	int numSamples, float sampleStep, int maxDiff, CostType costType, Func costFunc,
	const std::unordered_set<int>& frozenVerts = std::unordered_set<int>());

ManifoldMesh surfaceCutPlaneSep(const Volume<float>& cost, const ManifoldMesh& init,
	int numSamples, float sampleStep, int maxDiff, CostType costType, std::vector<Vec3f> planeNormals,
	std::vector<float> planeDists);
ManifoldMesh surfaceCutPlaneSep(const Volume<float>& cost, const ManifoldMesh&& init,
	int numSamples, float sampleStep, int maxDiff, CostType costType, std::vector<Vec3f> planeNormals,
	std::vector<float> planeDists);

std::vector<ManifoldMesh> surfaceCutPlaneSepQPBO(const Volume<float>& cost, std::vector<ManifoldMesh> meshes,
	int numSamples, float sampleStep, int maxDiff, CostType costType, const std::vector<Vec3f>& centers,
	const std::vector<std::vector<size_t>>& connections);

std::vector<ManifoldMesh> surfaceCutPlaneSepDual(const Volume<float>& cost, std::vector<ManifoldMesh> meshes,
	int numSamples, float sampleStep, int maxDiff, CostType costType, const std::vector<Vec3f>& centers,
	const std::vector<std::vector<size_t>>& connections);
std::vector<ManifoldMesh> surfaceCutPlaneSepDualParallel(const Volume<float>& cost,
	std::vector<ManifoldMesh> meshes, int numSamples, float sampleStep, int maxDiff, CostType costType,
	const std::vector<Vec3f>& centers, const std::vector<std::vector<size_t>>& connections, int numThreads);

std::vector<ManifoldMesh> kSurfaceCutNoOverlap(const std::vector<Volume<float>>& costs,
	const ManifoldMesh& init, int numSamples, float sampleStep, int maxDiff, int minSep, int maxSep,
	CostType costType);

#endif // SURFACE_SEGMENT_H__