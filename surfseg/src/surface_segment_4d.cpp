#include "surface_segment.h"
#include <tuple>
#include "util.h"

static float float_identity(const Volume4d<float>& vol, Vec4f pos)
{
	return vol.interp(pos);
}

FloatGraph& buildSurfaceGraph4d(FloatGraph& graph, const Volume<float>& costSamples,
	TetMesh4d& mesh, int maxDiff, CostType costType, const std::vector<int>& frozenVerts, size_t k = 0,
	size_t offset = 0);

TetMesh4d& updateVertices4d(const FloatGraph& graph, const Volume<float>& costSamples,
	TetMesh4d& mesh, float sampleStep, size_t k = 0, size_t offset = 0);

template <class Func>
Volume<float>& extractCostSamples4d(const Volume4d<float>& cost, const TetMesh4d& mesh,
	Volume<float>& samples, int numSamples, float sampleStep, Func costFunc, size_t k = 0);

TetMesh4d surfaceCut4d(const Volume4d<float>& cost, TetMesh4d mesh,
	int numSamples, float sampleStep, int maxDiff, CostType costType, const std::vector<int>& frozenVerts)
{
	return surfaceCut4d(cost, mesh, numSamples, sampleStep, maxDiff, costType, float_identity, frozenVerts);
}

template <class Func>
TetMesh4d surfaceCut4d(const Volume4d<float>& vol, TetMesh4d mesh,
	int numSamples, float sampleStep, int maxDiff, CostType costType, Func costFunc,
	const std::vector<int>& frozenVerts)
{
	mesh.computeVertexNormals(); // Ensure these are correct
	size_t numVerts = mesh.vertices.size();

	// Make cost sample volume
	Volume<float> costSamples(numVerts, numSamples, 1);
	costSamples.alloc();
	extractCostSamples4d(vol, mesh, costSamples, numSamples, sampleStep, costFunc);

	costSamples.saveToBin("costSamples.bin"); // Debug: save cost samples so we can look at them afterwards

	// Build min-cut graph and find optimal cut
	// TODO: Fix calculation of needed edges
	size_t totalSamples = costSamples.numElem();
	FloatGraph graph(totalSamples,
		mesh.edges.size() * (numSamples - maxDiff - 1) + totalSamples - costSamples.nx,
		graphErrFunc
	);
	graph.add_node(totalSamples);
	buildSurfaceGraph4d(graph, costSamples, mesh, maxDiff, costType, frozenVerts);

	graph.maxflow();

	// Update mesh vertex positions
	updateVertices4d(graph, costSamples, mesh, sampleStep);

	return mesh;
}

FloatGraph& buildSurfaceGraph4d(FloatGraph& graph, const Volume<float>& costSamples,
	TetMesh4d& mesh, int maxDiff, CostType costType, const std::vector<int>& frozenVerts, size_t k,
	size_t offset)
{
	using VertKey = TetMesh4d::VertKey;

	size_t numSamples = costSamples.ny;

	// Add edges
	for (const auto& v : mesh.vertices) {
		// Add intracolumn (downward) edges
		for (int i = numSamples - 1; i > 0; --i) {
			size_t ni = costSamples.idx(v.self, i, k) + offset;
			size_t nj = costSamples.idx(v.self, i - 1, k) + offset;
			graph.add_edge(ni, nj, infOrMax<float>, 0);
		}

		// Add intercolumn (neighbor) edges
		std::vector<VertKey> nborVerts;
		std::tie(nborVerts, std::ignore, std::ignore) = mesh.getVertLink(v.self);
		for (auto nk : nborVerts) {
			for (int i = numSamples - 1; i > maxDiff; --i) {
				size_t ni = costSamples.idx(v.self, i, k) + offset;
				size_t nj = costSamples.idx(nk, i - maxDiff, k) + offset;
				graph.add_edge(ni, nj, infOrMax<float>, 0);
			}
		}

		// Add source and sink edges
		float prevC = 0.0f;
		for (int i = 0; i < numSamples; ++i) {
			size_t ni = costSamples.idx(v.self, i, k) + offset;
			float c = costSamples.at(v.self, i, k);
			float w = costType == ON_SURFACE ? c - prevC : c;
			if (i == 0) {
				w = -1.0f;
			}
			if (w < 0) {
				graph.add_tweights(ni, -w, 0);
			}
			else {
				graph.add_tweights(ni, 0, w);
			}
			prevC = c;
		}
	}

	// Ensure all frozen vertices cannot move from their initial positions
	for (int vk : frozenVerts) {
		for (int i = 1; i < numSamples; ++i) {
			size_t ni = costSamples.idx(vk, i, k) + offset;
			graph.add_tweights(ni, 0, infOrMax<float>);
		}
	}

	return graph;
}

TetMesh4d& updateVertices4d(const FloatGraph& graph, const Volume<float>& costSamples,
	TetMesh4d& mesh, float sampleStep, size_t k, size_t offset)
{
	for (auto& v : mesh.vertices) {
		// Find upper position for this vertex
		for (int i = costSamples.ny - 1; i >= 0; --i) {
			size_t ni = costSamples.idx(v.self, i, k) + offset;
			if (graph.what_segment(ni) == SOURCE) {
				v.pos += i * sampleStep * v.normal;
				break;
			}
		}
	}

	return mesh;
}

template <class Func>
Volume<float>& extractCostSamples4d(const Volume4d<float>& vol, const TetMesh4d& mesh,
	Volume<float>& samples, int numSamples, float sampleStep, Func costFunc, size_t k)
{
	// Sample volume along vertex normals
	for (const auto& v : mesh.vertices) {
		Vec4f p = v.pos;
		Vec4f nrm = v.normal;
		for (int i = 0; i < numSamples; ++i) {
			//p[3] = fmaxf(0.0f, fminf(cost.nt - 1, p[3])); // Pad volume to repeat first and last time slice
			//p[3] = fmaxf(0.0f, p[3]); // Pad volume to repeat first slice
			//p[2] = fmaxf(0.0f, fminf(cost.nz - 1, p[2])); // Pad volume to repeat top and bottom slice
			samples.at(v.self, i, k) = vol.contains(p) ? costFunc(vol, p) : infOrMax<float>;
			p += sampleStep * nrm;
		}
	}

	return samples;
}
