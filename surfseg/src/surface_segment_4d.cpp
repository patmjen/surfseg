 #include "surface_segment.h"
#include <tuple>
#include <algorithm>
#include <thread>
#include <limits>
#include <iterator>
#include <GEL/CGLA/Mat3x3f.h>
#include <GEL/CGLA/Vec4d.h>
#include "util.h"
#include <GEL/Geometry/KDTree.h>

using Geometry::KDTree;

static float float_identity(const Volume4d<float>& vol, Vec4f pos)
{
	return vol.interp(pos);
}

static Vec4d f2dVec(Vec4f v)
{
	return Vec4d(
		static_cast<double>(v[0]),
		static_cast<double>(v[1]),
		static_cast<double>(v[2]),
		static_cast<double>(v[3])
	);
}

static Vec4f d2fVec(Vec4d v)
{
	return Vec4f(
		static_cast<float>(v[0]),
		static_cast<float>(v[1]),
		static_cast<float>(v[2]),
		static_cast<float>(v[3])
	);
}

FloatGraph& buildSurfaceGraph4d(FloatGraph& graph, const Volume<float>& costSamples,
	TetMesh4d& mesh, int maxDiff, CostType costType, const std::unordered_set<int>& frozenVerts, size_t k = 0,
	size_t offset = 0);

TetMesh4d& updateVertices4d(const FloatGraph& graph, const Volume<float>& costSamples,
	TetMesh4d& mesh, float sampleStep, size_t k = 0, size_t offset = 0);

TetMesh4d& updateVertices4d(const FloatGraph& graph, const Volume<float>& costSamples,
	TetMesh4d& mesh, const std::vector<Vec4f>& samplePositions, size_t k = 0, size_t offset = 0);

template <class Func>
void extractCostSamples4d(const Volume4d<float>& cost, const TetMesh4d& mesh,
	Volume<float>& samples, std::vector<Vec4f>& samplePositions, int numSamples, float sampleStep,
    Func costFunc, size_t k = 0);

template <class Func>
void extractCostSamples4dElf(const Volume4d<float>& vol, TetMesh4d& mesh, Volume<float>& samples,
	std::vector<Vec4f>& samplePositions, int numSamples, int numReverseSamples, double sampleStep,
	Func costFunc, const std::unordered_set<int>& frozenVerts, size_t k = 0, double tolP = 1e-4);

TetMesh4d surfaceCut4d(const Volume4d<float>& cost, TetMesh4d mesh,
	int numSamples, float sampleStep, int maxDiff, CostType costType,
	const std::unordered_set<int>& frozenVerts)
{
	return surfaceCut4d(cost, mesh, numSamples, sampleStep, maxDiff, costType, float_identity, frozenVerts);
}

template <class Func>
TetMesh4d surfaceCut4d(const Volume4d<float>& vol, TetMesh4d mesh,
	int numSamples, float sampleStep, int maxDiff, CostType costType, Func costFunc,
	const std::unordered_set<int>& frozenVerts)
{
	mesh.computeVertexNormals(); // Ensure these are correct
	size_t numVerts = mesh.vertices.size();

	const int numReverseSamples = 0; // TODO: make this an input

	// Make cost sample volume
	Volume<float> costSamples(numVerts, numSamples + numReverseSamples, 1);
	std::vector<Vec4f> samplePositions(numVerts * (numSamples + numReverseSamples));
	costSamples.alloc();
    //extractCostSamples4d(vol, mesh, costSamples, samplePositions, numSamples, sampleStep, costFunc);
	extractCostSamples4dElf(vol, mesh, costSamples, samplePositions, numSamples, numReverseSamples,
		sampleStep, costFunc, frozenVerts);

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
	updateVertices4d(graph, costSamples, mesh, samplePositions);

	return mesh;
}

FloatGraph& buildSurfaceGraph4d(FloatGraph& graph, const Volume<float>& costSamples,
	TetMesh4d& mesh, int maxDiff, CostType costType, const std::unordered_set<int>& frozenVerts, size_t k,
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

TetMesh4d& updateVertices4d(const FloatGraph& graph, const Volume<float>& costSamples,
	TetMesh4d& mesh, const std::vector<Vec4f>& samplePositions, size_t k, size_t offset)
{
	for (auto& v : mesh.vertices) {
		// Find upper position for this vertex
		for (int i = costSamples.ny - 1; i >= 0; --i) {
			size_t ni = costSamples.idx(v.self, i, k) + offset;
			if (graph.what_segment(ni) == SOURCE) {
				v.pos = samplePositions[ni];
				break;
			}
		}
	}

	return mesh;
}

template <class Func>
void extractCostSamples4d(const Volume4d<float>& vol, const TetMesh4d& mesh,
	Volume<float>& samples, std::vector<Vec4f>& samplePositions, int numSamples, float sampleStep,
    Func costFunc, size_t k)
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
            samplePositions[samples.idx(v.self, i, k)] = p;
			p += sampleStep * nrm;
		}
	}
}

template <class Func>
void extractCostSamples4dElf(const Volume4d<float>& vol, TetMesh4d& mesh, Volume<float>& samples,
	std::vector<Vec4f>& samplePositions, int numSamples, int numReverseSamples, double sampleStep,
	Func costFunc, const std::unordered_set<int>& frozenVerts, size_t k, double tolP)
{
	// Use sampling method from:
	// Yin et al., LOGISMOS-layered optimal graph image segmentation of multiple objects and surfaces...
	// 2010, IEEE Transactions on Medical Imaging
	using TetKey = TetMesh4d::TetKey;
	using VertKey = TetMesh4d::VertKey;
	using Vertex = TetMesh4d::Vertex;
	constexpr float div6 = 1.0f / 6.0f;

	if (numSamples == 0) {
		// Nothing to do, so just return now
		return;
	}

    // Precompute non-frozen vertices
    std::vector<int> nonFrozenVerts;
    nonFrozenVerts.reserve(mesh.vertices.size() - frozenVerts.size());
    for (const auto& v : mesh.vertices) {
        if (frozenVerts.find(v.self) == frozenVerts.end()) {
            nonFrozenVerts.push_back(v.self);
        }
    }

	// Precompute tet. volumes
	std::vector<float> tetVolumes;
	tetVolumes.reserve(mesh.tets.size());
	for (const auto& tet : mesh.tets) {
		tetVolumes.push_back(mesh.tetVolume(tet));
	}

    // Make kD-tree with vertex positions and their charges
	std::vector<float> charges;
	charges.reserve(mesh.vertices.size());
    KDTree<Vec4f, float> vertexKdTree;
    for (const auto& v : mesh.vertices) {
        std::vector<TetKey> nborTets;
        std::tie(std::ignore, std::ignore, nborTets) = mesh.getVertCoboundary(v.self);
        float q = 0;
        for (auto tk : nborTets) {
            q += tetVolumes[tk];
        }
        vertexKdTree.insert(v.pos, q);
        charges.push_back(q);
    }
    vertexKdTree.build();

	// Lambda function for evaluating field
	auto evalE = [&](Vec4d p) -> Vec4d
	{
        // Just using floats instead of GEL vectors is a lot faster, since std algorithms
        // are (apparently) not optimized out properly.
		Vec4f pf = d2fVec(p);
        float e0 = 0, e1 = 0, e2 = 0, e3 = 0;
		for (const auto& v : mesh.vertices) {
            const float q = charges[v.self];
            const float rv0 = pf[0] - v.pos[0];
            const float rv1 = pf[1] - v.pos[1];
            const float rv2 = pf[2] - v.pos[2];
            const float rv3 = pf[3] - v.pos[3];
			
            // NOTE: Using fast inverse square root seems to be ~10% faster, but is less accurate
            const float r2 = rv0 * rv0 + rv1 * rv1 + rv2 * rv2 + rv3 * rv3;
			const float r = sqrtf(r2);
			const float r4 = r2 * r2;
            const float r5 = r4 * r;
            const float f = q / r5;
            e0 += f * rv0;
            e1 += f * rv1;
            e2 += f * rv2;
            e3 += f * rv3;
		}
        return Vec4d(e0, e1, e2, e3);
	};

    // Lambda function for sampling single EFL
	auto sampleEfl = [&](const Vertex& v, int startSample, int numSamples, bool reverse)
	{
		Vec4d nrm = f2dVec(v.normal);
		Vec4d p = f2dVec(v.pos);
		Vec4f pf = v.pos;

		// Perturp position along the normal direction so the EFL points in the right direction
		const double perturb = sampleStep * 0.002;
		p += nrm * (reverse ? -1.0 : 1.0) * perturb;

		// Initialize h so the first step will have length close to 0.001 * sampleStep
		// Empirically, this allows us to accept a step in 1-2 attempts
		double perturb4 = perturb * perturb;
		perturb4 *= perturb4;
		double h = 2 * 0.001 * sampleStep * perturb4 / charges[v.self];
        int maxIter = 1000; // Max integration steps (usually needs less than 50)
        int iter; // Current integration iteration
		const double minStep = h * 6e-8; // We allow approx. 24 halvings from the initial step
		for (int i = 0; i < numSamples; ++i) {
			double crntLen = 0;
			// Integrate field equations using RKF45 with adaptive step (step doubling)
            for (iter = 0; iter < maxIter; ++iter) {
				// Define all RKF45 coefficients
				constexpr double a21 = 0.25;
				constexpr double a31 = 3.0 / 32.0;
				constexpr double a32 = 9.0 / 32.0;
				constexpr double a41 = 1932.0 / 2197.0;
				constexpr double a42 = -7200.0 / 2197.0;
				constexpr double a43 = 7296.0 / 2197.0;
				constexpr double a51 = 439.0 / 216.0;
				constexpr double a52 = -8.0;
				constexpr double a53 = 3680.0 / 513.0;
				constexpr double a54 = -845.0 / 4104.0;
				constexpr double a61 = -8.0 / 27.0;
				constexpr double a62 = 2.0;
				constexpr double a63 = -3544.0 / 2565.0;
				constexpr double a64 = 1859.0 / 4104.0;
				constexpr double a65 = -11.0 / 40.0;

				constexpr double s41 = 25.0 / 216.0;
				constexpr double s43 = 1408.0 / 2565.0;
				constexpr double s44 = 2197.0 / 4104.0;
				constexpr double s45 = -1.0 / 5.0;

				constexpr double s51 = 16.0 / 135.0;
				constexpr double s53 = 6656.0 / 12825.0;
				constexpr double s54 = 28561.0 / 56430.0;
				constexpr double s55 = -9.0 / 50.0;
				constexpr double s56 = 2.0 / 55.0;

				// Compute 4th and 5th order steps (dp4 and dp5, respectively)
				const Vec4d k1 = h * evalE(p);
				const Vec4d k2 = h * evalE(p + a21 * k1);
				const Vec4d k3 = h * evalE(p + a31 * k1 + a32 * k2);
				const Vec4d k4 = h * evalE(p + a41 * k1 + a42 * k2 + a43 * k3);
				const Vec4d k5 = h * evalE(p + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4);
				const Vec4d k6 = h * evalE(p + a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5);

				const Vec4d dp4 = s41 * k1 + s43 * k3 + s44 * k4 + s45 * k5;
				const Vec4d dp5 = s51 * k1 + s53 * k3 + s54 * k4 + s55 * k5 + s56 * k6;

				// If the estimated error is too large we halve the step and try again
				if (length(dp4 - dp5) > tolP && h > minStep) {
					h *= 0.5;
					continue;
				}

				const double dpLen = length(dp4);

				/*if (!vol.contains(d2fVec(p + dp4))) {
					// Truncate step if we would exit the volume
					Vec4d a0 = -p / dp4;
					Vec4d a1 = (Vec4d(vol.nx - 1, vol.ny - 1, vol.nz - 1, vol.nt - 1) - p) / dp4;
					double a = std::min(1.0, v_max(a0, a1).min_coord());

					// We can't take more steps, so stop iterations now
					p += a * dp4;
					break;
				}*/

				if (crntLen + dpLen >= sampleStep) {
					// Truncate step so we don't go too far
					// TODO: Maybe find better way to acheive this
					double a = (sampleStep - crntLen)/dpLen;
					p += a * dp4;

					// We are now at the next sample position, so stop iterations and reset
					break;
				}
				p += dp4;
				crntLen += dpLen;
				// Step worked so double for next iteration.
				h *= 2;
			}
            if (iter == maxIter) {
                // If we maxed out the iterations, then we should just stay here for the remaining steps
                maxIter = 0;
            }

			pf = d2fVec(p);
			const int idx = reverse ? startSample - i : startSample + i;
			samplePositions[samples.idx(v.self, idx, k)] = pf;
            pf[3] = std::max(0.0, std::min(static_cast<double>(vol.nt) - 1.0, p[3])); // Pad time
            //pf[2] = std::max(0.0, std::min(static_cast<double>(vol.nz) - 1.0, p[2])); // Pad z
			samples.at(v.self, idx, k) = vol.contains(pf) ? costFunc(vol, pf) : infOrMax<float>;
		}
	};

	// Sample volume along electrical field lines (EFLs)
	auto sampleSubset = [&](int begin, int end)
	{
		for (int vi = begin; vi < end; ++vi) {
			// We use doubles in this loop for extra precision
			const Vertex& v = mesh.vertices[nonFrozenVerts[vi]];
			Vec4d nrm = f2dVec(v.normal);
			Vec4d p = f2dVec(v.pos);
			Vec4f pf = v.pos;

			// Take first sample now, before perturbing p
			samples.at(v.self, numReverseSamples, k) = vol.contains(pf) ? costFunc(vol, pf) : infOrMax<float>;
			samplePositions[samples.idx(v.self, numReverseSamples, k)] = pf;

			sampleEfl(v, numReverseSamples + 1, numSamples - 1, false);
			sampleEfl(v, numReverseSamples - 1, numReverseSamples, true);
		}
	};

    int numThreads = std::thread::hardware_concurrency();
	int numVert = nonFrozenVerts.size();
	if (numThreads == 1) {
		// Just run on this thread
		sampleSubset(0, numVert);
	} else {
		// Split sampling for non-frozen vertices among threads
		int verticesPerThread = (numVert + (numVert % numThreads == 0 ? 0 : numThreads)) / numThreads;
		std::vector<std::thread> threads;
		for (int i = 0; i < numThreads; ++i) {
			int begin = i * verticesPerThread;
			int end = std::min((i + 1) * verticesPerThread, numVert);
			threads.push_back(std::thread(sampleSubset,begin, end));
		}
		// Wait for all threads to finish
		for (auto& th : threads) {
			if (th.joinable()) {
				th.join();
			}
		}
	}
    // Handle all frozen vertices in this thread, since they require almost no work
    for (int vi : frozenVerts) {
        const Vertex& v = mesh.vertices[vi];
        Vec4d p = f2dVec(v.pos);
        Vec4f pf = v.pos;

        const float val = vol.contains(pf) ? costFunc(vol, pf) : infOrMax<float>;
        for (int i = 0; i < numSamples + numReverseSamples; ++i) {
            samples.at(v.self, i, k) = val;
            samplePositions[samples.idx(v.self, i, k)] = pf;
        }
    }

    // DEBUG BEGIN
    Volume<float> samplePositionsVol(samples.nx, samples.ny, 4);
    samplePositionsVol.alloc();
    for (const auto& v : mesh.vertices) {
        for (int i = 0; i < numSamples + numReverseSamples; ++i) {
            Vec4f p = samplePositions[samples.idx(v.self, i, k)];
            samplePositionsVol.at(v.self, i, 0) = p[0];
            samplePositionsVol.at(v.self, i, 1) = p[1];
            samplePositionsVol.at(v.self, i, 2) = p[2];
            samplePositionsVol.at(v.self, i, 3) = p[3];
        }
    }
	samplePositionsVol.saveToBin("samplePositionsVol.bin");
    // DEBUG END
}
  