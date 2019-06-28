#include "surface_segment.h"
#include "util.h"
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <tuple>

ManifoldMesh& updateVertices(const FloatGraph& graph, const Volume<float>& costSamples,
	ManifoldMesh& mesh, float sampleStep, size_t k = 0, size_t offset = 0);

ManifoldMesh& updateVerticesDual(const FloatGraph& graph, const Volume<float>& costSamples,
	ManifoldMesh& mesh, float sampleStep, size_t k = 0, size_t offset = 0);

ManifoldMesh& updateVerticesQPBO(const FloatGraph& graph, const Volume<float>& costSamples,
	ManifoldMesh& mesh, float sampleStep, size_t kp = 0, size_t kd = 1, size_t offset = 0);

FloatGraph& buildSurfaceGraph(FloatGraph& graph, const Volume<float>& costSamples,
	const ManifoldMesh& mesh, int maxDiff, CostType costType, size_t k = 0, size_t offset = 0);

FloatGraph& buildDualSurfaceGraph(FloatGraph& graph, const Volume<float>& costSamples,
	const ManifoldMesh& mesh, int maxDiff, CostType costType, size_t k = 0, size_t offset = 0);

FloatGraph& buildQPBOSurfaceGraph(FloatGraph& graph, const Volume<float>& costSamples,
	const ManifoldMesh& mesh, int maxDiff, CostType costType, size_t kp = 0, size_t kd = 1, size_t offset = 0);

Volume<float>& extractCostSamples(const Volume<float>& cost, const ManifoldMesh& mesh,
	Volume<float>& samples, int numSamples, float sampleStep, size_t k = 0);

ManifoldMesh surfaceCut(const Volume<float>& cost, const ManifoldMesh& init,
    int numSamples, float sampleStep, int maxDiff, CostType costType)
{
    ManifoldMesh mesh(init);
    return surfaceCut(cost, std::move(mesh), numSamples, sampleStep, maxDiff, costType);
}

ManifoldMesh surfaceCut(const Volume<float>& cost, const ManifoldMesh&& init,
    int numSamples, float sampleStep, int maxDiff, CostType costType)
{
    ManifoldMesh mesh(std::move(init));
    mesh.computeVertexNormals(false, true); // Ensure these are correct
    size_t numNodes = init.vertices.size();

    // Make cost sample volume
    Volume<float> costSamples(numNodes, numSamples, 1);
    costSamples.alloc();
    extractCostSamples(cost, mesh, costSamples, numSamples, sampleStep);

    // Build min-cut graph and find optimal cut
    size_t totalSamples = costSamples.numElem();
    FloatGraph graph(totalSamples,
        mesh.edges.size() * (numSamples - maxDiff - 1) + totalSamples - costSamples.nx, 
		graphErrFunc
	);
    graph.add_node(totalSamples);
    buildSurfaceGraph(graph, costSamples, mesh, maxDiff, costType);

    graph.maxflow();

    // Update mesh vertex positions
	updateVertices(graph, costSamples, mesh, sampleStep);

    return mesh;
}


ManifoldMesh surfaceCutPlaneSep(const Volume<float>& cost, const ManifoldMesh& init,
	int numSamples, float sampleStep, int maxDiff, CostType costType, std::vector<Vec3f> planeNormals,
	std::vector<float> planeDists)
{
	ManifoldMesh mesh(init);
	return surfaceCutPlaneSep(cost, std::move(mesh), numSamples, sampleStep, maxDiff, costType,
		planeNormals, planeDists);
}

ManifoldMesh surfaceCutPlaneSep(const Volume<float>& cost, const ManifoldMesh&& init,
	int numSamples, float sampleStep, int maxDiff, CostType costType, std::vector<Vec3f> planeNormals,
	std::vector<float> planeDists)
{
	assert(planeNormals.size() == planeDists.size());

	// TODO: This is very similar to surfaceCut - need to refactor
	ManifoldMesh mesh(std::move(init));
	mesh.computeVertexNormals(); // Ensure these are correct
	size_t numNodes = init.vertices.size();

	// Make cost sample volume
	Volume<float> costSamples(numNodes, numSamples, 1);
	costSamples.alloc();
	extractCostSamples(cost, mesh, costSamples, numSamples, sampleStep);

	// Build min-cut graph
	size_t totalSamples = costSamples.numElem();
	FloatGraph graph(totalSamples,
		mesh.edges.size() * (numSamples - maxDiff - 1) + totalSamples - costSamples.nx, 
		graphErrFunc
	);
	graph.add_node(totalSamples);
	buildSurfaceGraph(graph, costSamples, mesh, maxDiff, costType);

	// Add plane constraints
	for (const auto& v : mesh.vertices) {
		for (int pi = 0; pi < planeNormals.size(); ++pi) {
			Vec3f nrm = planeNormals[pi];
			float d = planeDists[pi];
			Vec3f p = mesh.vpos(v);
			for (int i = 0; i < numSamples; ++i) {
				if (dot(nrm, p) > d) {
					size_t ni = costSamples.idx(v.self, i, 0);
					graph.add_tweights(ni, 0, infOrMax<float>);
				}
				p += sampleStep * mesh.vnormal(v);
			}
		}
	}

	graph.maxflow();

	// Update mesh vertex positions
	updateVertices(graph, costSamples, mesh, sampleStep);

	return mesh;
}

std::vector<ManifoldMesh> surfaceCutPlaneSepQPBO(const Volume<float>& cost, std::vector<ManifoldMesh> meshes,
	int numSamples, float sampleStep, int maxDiff, CostType costType, const std::vector<Vec3f>& centers,
	const std::vector<std::vector<size_t>>& connections)
{
	assert(meshes.size() == centers.size());
	assert(meshes.size() == connections.size());
	for (auto& m : meshes) {
		m.computeVertexNormals(false, true); // Ensure these are correct
	}

	// Make cost sample volumes
	std::vector<Volume<float>> costSamples;
	costSamples.reserve(meshes.size());
	size_t totalSamples = 0;
	size_t totalEdges = 0;
	for (const auto& m : meshes) {
		size_t numNodes = m.vertices.size();
		Volume<float> cs(numNodes, numSamples, 2);
		cs.alloc();
		extractCostSamples(cost, m, cs, numSamples, sampleStep, 0);
		extractCostSamples(cost, m, cs, numSamples, sampleStep, 1);
		size_t layerLen = numNodes * numSamples;

		totalSamples += cs.numElem();
		totalEdges += m.edges.size() * (numSamples - maxDiff - 1) + cs.numElem() - cs.nx;
		costSamples.push_back(std::move(cs));
	}

	// Compute number of needed nodes for the surface to plane stuff
	size_t extraNodes = 0;
	for (int i = 0; i < meshes.size(); ++i) {
		Vec3f ceni = centers[i];
		for (int j : connections[i]) {
			Vec3f cenj = centers[j];
			extraNodes += 2 * roundf(length(ceni - cenj) / sampleStep);
		}
	}

	// Build surface graphs
	std::vector<size_t> offsets;
	offsets.reserve(meshes.size());
	// Adding two extra totalEdges is a heuristic for how many surface to plane position edges we are gonna need
	// It seems to be a good upper bound, but there are no guarantees
	const size_t heuristic = totalEdges;
	FloatGraph graph(totalSamples + extraNodes, heuristic + 2 * totalEdges + 2 * extraNodes, graphErrFunc);

	graph.add_node(totalSamples + extraNodes);
	for (size_t i = 0, offset = 0; i < meshes.size(); ++i) {
		auto& cs = costSamples[i];
		const auto& m = meshes[i];

		buildQPBOSurfaceGraph(graph, cs, m, maxDiff, costType, 0, 1, offset);

		offsets.push_back(offset);
		offset += cs.numElem();
		// From this point we don't need the actual samples anymore, but we want to keep the Volume
		// for computing indices
		cs.data = nullptr;
	}

	// Add plane edges
	for (int i = 0, planeNodeOffset = totalSamples; i < meshes.size(); ++i) {
		Vec3f ceni = centers[i];
		for (int j : connections[i]) {
			Vec3f cenj = centers[j];
			int numPlaneNodes = roundf(length(ceni - cenj) / sampleStep);
			Vec3f normal = normalize(cenj - ceni);
			// Add intracolumn plane edges
			for (int n = 0; n < numPlaneNodes; ++n) {
				if (n > 0) {
					size_t nip = planeNodeOffset + n;
					size_t nid = planeNodeOffset + numPlaneNodes + n;
					size_t njp = planeNodeOffset + n - 1;
					size_t njd = planeNodeOffset + numPlaneNodes + n - 1;
					graph.add_edge(nip, njp, infOrMax<float>, 0);
					graph.add_edge(njd, nid, infOrMax<float>, 0);
				}
			}
			// Add edges from surface to plane
			for (const auto& v : meshes[i].vertices) {
				Vec3f p = meshes[i].vpos(v);
				Vec3f nrm = meshes[i].vnormal(v);
				if (dot(normal, p) <= dot(normal, ceni) && dot(normal, nrm) <= 0) {
					// Vertex is behind the plane and normal points away, so we can skip it
					continue;
				}
				for (int n = 0, si = 0; n < numPlaneNodes; ++n) {
					// We don't reset si, since if a point was behind the plane for some n, then it will also
					// be behind it for n + 1
					float dist = dot(normal, ceni + n * sampleStep * normal);
					for (; si < numSamples; ++si) {
						if (dot(normal, p) >= dist) {
							size_t nip = costSamples[i].idx(v.self, si, 0) + offsets[i];
							size_t nid = costSamples[i].idx(v.self, si, 1) + offsets[i];
							size_t njp = planeNodeOffset + n;
							size_t njd = planeNodeOffset + numPlaneNodes + n;
							graph.add_edge(nip, njp, infOrMax<float>, 0);
							graph.add_edge(njd, nid, infOrMax<float>, 0);
							break;
						}
						p += sampleStep * nrm;
					}
				}
			}
			// Add edges from plane to surface
			for (const auto& v : meshes[j].vertices) {
				Vec3f p = meshes[j].vpos(v);
				Vec3f nrm = meshes[j].vnormal(v);
				if (dot(normal, p) >= dot(normal, ceni) && dot(normal, nrm) >= 0) {
					// Vertex is in front of the plane and normal points away, so we can skip it
					continue;
				}
				for (int n = numPlaneNodes - 1, si = 0; n >= 0; --n) {
					// We don't reset si, since if a point was in front of the plane for some n, then it will
					// also be in front of it for n + 1
					float dist = dot(normal, ceni + n * sampleStep * normal);
					for (; si < numSamples; ++si) {
						if (dot(normal, p) <= dist) {
							size_t nip = costSamples[j].idx(v.self, si, 0) + offsets[j];
							size_t nid = costSamples[j].idx(v.self, si, 1) + offsets[j];
							size_t njp = planeNodeOffset + n;
							size_t njd = planeNodeOffset + numPlaneNodes + n;
							graph.add_edge(nip, njd, infOrMax<float>, 0);
							graph.add_edge(njp, nid, infOrMax<float>, 0);
							break;
						}
						p += sampleStep * nrm;
					}
				}
			}
			planeNodeOffset += 2 * numPlaneNodes;
		}
	}
	
	graph.maxflow();

	// Update mesh vertex positions
	for (size_t i = 0; i < meshes.size(); ++i) {
		const auto& cs = costSamples[i];
		auto& m = meshes[i];
		updateVerticesQPBO(graph, cs, m, sampleStep, 0, 1, offsets[i]);
	}
	return meshes;
}

std::vector<ManifoldMesh> surfaceCutPlaneSepDual(const Volume<float>& cost, std::vector<ManifoldMesh> meshes,
	int numSamples, float sampleStep, int maxDiff, CostType costType, const std::vector<Vec3f>& centers,
	const std::vector<std::vector<size_t>>& connections)
{
	assert(meshes.size() == centers.size());
	assert(meshes.size() == connections.size());
	for (auto& m : meshes) {
		m.computeVertexNormals(false, true); // Ensure these are correct
	}

	// Make cost sample volumes
	std::vector<Volume<float>> costSamples;
	costSamples.reserve(meshes.size());
	size_t totalSamples = 0;
	size_t totalEdges = 0;
	for (const auto& m : meshes) {
		size_t numNodes = m.vertices.size();
		Volume<float> cs(numNodes, numSamples, 2);
		cs.alloc();
		extractCostSamples(cost, m, cs, numSamples, sampleStep, 0);
		extractCostSamples(cost, m, cs, numSamples, sampleStep, 1);
		size_t layerLen = numNodes * numSamples;

		totalSamples += cs.numElem();
		totalEdges += m.edges.size() * (numSamples - maxDiff - 1) + cs.numElem() - cs.nx;
		costSamples.push_back(std::move(cs));
	}

	// Compute number of needed nodes for the surface to plane stuff
	size_t extraNodes = 0;
	for (int i = 0; i < meshes.size(); ++i) {
		Vec3f ceni = centers[i];
		for (int j : connections[i]) {
			Vec3f cenj = centers[j];
			extraNodes += roundf(length(ceni - cenj) / sampleStep);
		}
	}

	// Build surface graphs
	std::vector<size_t> offsets;
	offsets.reserve(meshes.size());
	// Adding an extra totalEdges is a heuristic for how many surface to plane position edges we are gonna need
	// It seems to be a good upper bound, but there are no guarantees
	FloatGraph graph(totalSamples + extraNodes, 3 * totalEdges + 2 * extraNodes, graphErrFunc);

	graph.add_node(totalSamples + extraNodes);
	for (size_t i = 0, offset = 0; i < meshes.size(); ++i) {
		auto& cs = costSamples[i];
		const auto& m = meshes[i];
		buildSurfaceGraph(graph, cs, m, maxDiff, costType, 0, offset);
		buildDualSurfaceGraph(graph, cs, m, maxDiff, costType, 1, offset);
		offsets.push_back(offset);
		offset += cs.numElem();
		// From this point we don't need the actual samples anymore, but we want to keep the Volume
		// for computing indices
		cs.data = nullptr;
	}

	// Add plane edges
	for (int i = 0, planeNodeOffset = totalSamples; i < meshes.size(); ++i) {
		Vec3f ceni = centers[i];
		for (int j : connections[i]) {
			Vec3f cenj = centers[j];
			int numNodes = roundf(length(ceni - cenj) / sampleStep);
			Vec3f normal = normalize(cenj - ceni);
			// Add intracolumn plane edges
			for (int n = 0; n < numNodes; ++n) {
				if (n > 0) {
					size_t ni = planeNodeOffset + n;
					size_t nj = planeNodeOffset + n - 1;
					graph.add_edge(ni, nj, infOrMax<float>, 0);
				}
			}
			// Add edges from primary surface i to plane positions
			for (const auto& v : meshes[i].vertices) {
				Vec3f p = meshes[i].vpos(v);
				Vec3f nrm = meshes[i].vnormal(v);
				if (dot(normal, p) <= dot(normal, ceni) && dot(normal, nrm) <= 0) {
					// Vertex is behind the plane and normal points away, so we can skip it
					continue;
				}
				for (int n = 0, si = 0; n < numNodes; ++n) {
					// We don't reset si, since if a point was behind the plane for some n, then it will also
					// be behind it for n + 1
					float dist = dot(normal, ceni + n * sampleStep * normal);
					for (; si < numSamples; ++si) {
						if (dot(normal, p) > dist) {
							size_t ni = costSamples[i].idx(v.self, si, 0) + offsets[i];
							size_t nj = planeNodeOffset + n;
							graph.add_edge(ni, nj, infOrMax<float>, 0);
							break;
						}
						p += sampleStep * nrm;
					}
				}
			}
			// Add edges from dual surface j to plane positions
			for (const auto& v : meshes[j].vertices) {
				Vec3f p = meshes[j].vpos(v);
				Vec3f nrm = meshes[j].vnormal(v);
				if (dot(normal, p) >= dot(normal, ceni) && dot(normal, nrm) >= 0) {
					// Vertex is in front of the plane and normal points away, so we can skip it
					continue;
				}
				for (int n = numNodes - 1, si = 0; n >= 0; --n) {
					// We don't reset si, since if a point was in front of the plane for some n, then it will
					// also be in front of it for n + 1
					float dist = dot(normal, ceni + n * sampleStep * normal);
					for (; si < numSamples; ++si) {
						if (dot(normal, p) < dist) {
							size_t ni = planeNodeOffset + n;
							size_t nj = costSamples[j].idx(v.self, si, 1) + offsets[j];
							graph.add_edge(ni, nj, infOrMax<float>, 0);
							break;
						}
						p += sampleStep * nrm;
					}
				}
			}
			planeNodeOffset += numNodes;
		}
	}

	graph.maxflow();

	// Update mesh vertex positions
	for (size_t i = 0; i < meshes.size(); ++i) {
		const auto& cs = costSamples[i];
		auto& m = meshes[i];
		updateVertices(graph, cs, m, sampleStep, 0, offsets[i]);
		//updateVerticesDual(graph, cs, m, sampleStep, 1, offsets[i]);
	}
	return meshes;
}

std::vector<ManifoldMesh> surfaceCutPlaneSepDualParallel(const Volume<float>& cost,
	std::vector<ManifoldMesh> meshes, int numSamples, float sampleStep, int maxDiff, CostType costType,
	const std::vector<Vec3f>& centers, const std::vector<std::vector<size_t>>& connections, int numThreads)
{
	std::ofstream flog("prog-log.txt");
	assert(flog);

	if (numThreads == 1) {
		// For one thread just run the single threaded version
		return surfaceCutPlaneSepDual(cost, meshes, numSamples, sampleStep, maxDiff, costType, centers, connections);
	}

	flog << "Building graphs" << std::endl;
	assert(meshes.size() == centers.size());
	assert(meshes.size() == connections.size());
	int numMesh = meshes.size();
	for (auto& m : meshes) {
		m.computeVertexNormals(false, true); // Ensure these are correct
	}

	size_t totalExtraNodes = 0;

	std::unordered_map<std::pair<int, int>, std::pair<int, int>> sharedNodes;

	std::vector<Volume<float>> costSamples;
	std::vector<FloatGraph> primalGraphs;
	std::vector<FloatGraph> dualGraphs;
	costSamples.reserve(numMesh);
	primalGraphs.reserve(numMesh);
	dualGraphs.reserve(numMesh);
	// Make surface graphs
	for (int i = 0; i < numMesh; ++i) {
		const auto& mesh = meshes[i];
		// Make cost sample volume
		size_t numNodes = mesh.vertices.size();
		Volume<float> samples(numNodes, numSamples, 1);
		samples.alloc();
		extractCostSamples(cost, mesh, samples, numSamples, sampleStep);
		extractCostSamples(cost, mesh, samples, numSamples, sampleStep);

		// Compute number of nodes and edges
		Vec3f ceni = centers[i];
		size_t extraNodes = 0;
		for (int j : connections[i]) {
			Vec3f cenj = centers[j];
			extraNodes += roundf(length(ceni - cenj) / sampleStep);
		}
		totalExtraNodes += extraNodes;
		size_t surfaceNodes = samples.numElem();
		size_t surfaceEdges = mesh.edges.size() * (numSamples - maxDiff - 1) + surfaceNodes - numNodes;
		// We must use emplace_back so the Graph is constructed in-place in the vector, otherwise a temporary
		// object is made which will have its destructor called which frees the memory
		primalGraphs.emplace_back(surfaceNodes + extraNodes, surfaceEdges + surfaceEdges / 2, graphErrFunc);
		dualGraphs.emplace_back(surfaceNodes + extraNodes, surfaceEdges + surfaceEdges / 2, graphErrFunc);

		// Build graphs
		FloatGraph& primalGraph = primalGraphs.back();
		FloatGraph& dualGraph = dualGraphs.back();
		primalGraph.add_node(surfaceNodes + extraNodes);
		dualGraph.add_node(surfaceNodes + extraNodes);

		buildSurfaceGraph(primalGraph, samples, mesh, maxDiff, costType);
		buildDualSurfaceGraph(dualGraph, samples, mesh, maxDiff, costType);

		// Add plane position edges
		int planeNodeOffset = surfaceNodes;
		for (int j : connections[i]) {
			Vec3f cenj = centers[j];
			int numNodes = roundf(length(ceni - cenj) / sampleStep); // TODO: Rename this as it shadows a variable 
			Vec3f normal = normalize(cenj - ceni);
			// Add intracolumn plane edges and Lagrange edges
			for (int n = 0; n < numNodes; ++n) {
				// Initially, all Lagrange multipliers are 0
				primalGraph.add_tweights(planeNodeOffset + n, 0.0f, 0.0f);
				dualGraph.add_tweights(planeNodeOffset + n, 0.0f, 0.0f);
				if (n > 0) {
					size_t ni = planeNodeOffset + n;
					size_t nj = planeNodeOffset + n - 1;
					primalGraph.add_edge(ni, nj, infOrMax<float>, 0);

					ni = planeNodeOffset + numNodes - n - 2;
					nj = planeNodeOffset + numNodes - n - 1;
					dualGraph.add_edge(ni, nj, infOrMax<float>, 0);
				}
			}
			// Add edges from primary surface to plane positions
			for (const auto& v : mesh.vertices) {
				Vec3f p = mesh.vpos(v);
				Vec3f nrm = mesh.vnormal(v);
				if (dot(normal, p) <= dot(normal, ceni) && dot(normal, nrm) <= 0) {
					// Vertex is behind the plane and normal points away, so we can skip it
					continue;
				}
				for (int n = 0, si = 0; n < numNodes; ++n) {
					// We don't reset si, since if a point was behind the plane for some n, then it will also
					// be behind it for n + 1
					float dist = dot(normal, ceni + n * sampleStep * normal);
					for (; si < numSamples; ++si) {
						if (dot(normal, p) > dist) {
							size_t ni = samples.idx(v.self, si, 0);
							size_t nj = planeNodeOffset + n;
							primalGraph.add_edge(ni, nj, infOrMax<float>, 0);
							break;
						}
						p += sampleStep * nrm;
					}
				}
			}
			// Add edges from dual surface to plane positions
			for (const auto& v : mesh.vertices) {
				Vec3f p = mesh.vpos(v);
				Vec3f nrm = mesh.vnormal(v);
				if (dot(normal, p) >= dot(normal, cenj) && dot(normal, nrm) >= 0) {
					// Vertex is in front of the plane and normal points away, so we can skip it
					continue;
				}
				for (int n = 0, si = 0; n < numNodes; ++n) {
					// We don't reset si, since if a point was in front of the plane for some n, then it will
					// also be in front of it for n + 1
					float dist = dot(normal, ceni + n * sampleStep * normal);
					for (; si < numSamples; ++si) {
						if (dot(normal, p) < dist) {
							size_t ni = planeNodeOffset + n;
							size_t nj = samples.idx(v.self, si, 0);
							dualGraph.add_edge(ni, nj, infOrMax<float>, 0);
							break;
						}
						p += sampleStep * nrm;
					}
				}
			}
			sharedNodes[std::make_pair(i, j)] = std::make_pair(planeNodeOffset, numNodes);
			planeNodeOffset += numNodes;
		}
		costSamples.push_back(std::move(samples));
	}

	flog << "Compute pairings" << std::endl;

	// Precompute indices for all pairs
	// TODO: This seems quite wasteful -- try to remove
	std::vector<std::tuple<int, int, int, int>> pairings;
	pairings.reserve(totalExtraNodes);
	for (int i = 0; i < numMesh; ++i) {
		for (int j : connections[i]) {
			int starti, numi, startj, numj;
			std::tie(starti, numi) = sharedNodes[std::make_pair(i, j)];
			std::tie(startj, numj) = sharedNodes[std::make_pair(j, i)];
			assert(numi == numj);
			flog << "Pairing " << i << " <-> " << j << std::endl;
			for (int n = 0; n < numi; ++n) {
				flog << "i: " << starti + n << ", j: " << startj + n << std::endl;
				pairings.push_back(std::make_tuple(i, j, starti + n, startj + n));
			}
		}
	}

	flog << "Asserting sizes " << pairings.size() << ", " << totalExtraNodes << std::endl;
	assert(pairings.size() == totalExtraNodes);

	std::vector<float> lagrangeMults(totalExtraNodes, 0.0f);
	std::vector<float> steps(totalExtraNodes, 0.5f);
	std::vector<bool> signs(totalExtraNodes, false); // false = negative, true = positive
	
	flog << "Run opti loop" << std::endl;
	// Run optimization loop
	for (int iter = 0; iter < 1000; ++iter) {
		// Solve subproblems
		flog << "Solve subgraphs" << std::endl;
		for (auto& g : primalGraphs) {
			g.maxflow();
		}
		for (auto& g : dualGraphs) {
			g.maxflow();
		}
		flog << "Check convergence" << std::endl;
		// Check for convergence and do step
		bool converged = true;
		for (int i = 0; i < pairings.size(); ++i) {
			const auto& p = pairings[i];
			size_t si = std::get<0>(p);
			size_t sj = std::get<1>(p);
			size_t ni = std::get<2>(p);
			size_t nj = std::get<3>(p);
			int xi = primalGraphs[si].what_segment(ni);
			int yi = dualGraphs[sj].what_segment(nj);
			flog << "Check si: " << si << ", sj: " << sj << ", ni: " << ni << ", nj: " << nj << ", xi: " << xi << ", yi: " << yi << std::endl;
			if (xi != yi) {
				converged = false;
				float oldMult = lagrangeMults[i];
				lagrangeMults[i] += steps[i] * (xi - yi);
				bool positive = xi - yi > 0;
				if (positive != signs[i]) {
					steps[i] *= 0.5;
				}
				signs[i] = positive;
				if (lagrangeMults[i] < 0) {
					primalGraphs[si].add_tweights(ni, 0, oldMult - lagrangeMults[ni]);
					dualGraphs[sj].add_tweights(nj, -lagrangeMults[nj], oldMult);
				} else {
					primalGraphs[si].add_tweights(nj, lagrangeMults[ni], oldMult);
					dualGraphs[sj].add_tweights(ni, 0, lagrangeMults[nj] + oldMult);
				}
				flog << "Caps: " << primalGraphs[si].get_trcap(ni) << dualGraphs[sj].get_trcap(nj) << std::endl;
				flog << "Lagrange: " << lagrangeMults[i] << std::endl;
				primalGraphs[si].mark_node(ni);
				dualGraphs[sj].mark_node(nj);
			}
		}
		if (converged) {
			break;
		}
	}

	for (size_t i = 0; i < meshes.size(); ++i) {
		const auto& cs = costSamples[i];
		auto& m = meshes[i];
		//updateVertices(primalGraphs[i], cs, m, sampleStep, 0);
		updateVerticesDual(dualGraphs[i], cs, m, sampleStep, 0);
	}

	return meshes;
}

std::vector<ManifoldMesh> kSurfaceCutNoOverlap(const std::vector<Volume<float>>& costs,
    const ManifoldMesh& init, int numSamples, float sampleStep, int maxDiff, int minSep, int maxSep,
	CostType costType)
{
    // TODO: Maybe generalize this to different initial meshes
    size_t k = costs.size();
    size_t numNodes = init.vertices.size();
    size_t numEdges = init.edges.size();

    std::vector<ManifoldMesh> meshes(k, init);
    for (auto& m : meshes) {
        m.computeVertexNormals(false, true);
    }

    // Make cost sample volume
    Volume<float> costSamples(numNodes, numSamples, k);
    costSamples.alloc();
    for (int i = 0; i < k; ++i) {
        extractCostSamples(costs[i], meshes[i], costSamples, numSamples, sampleStep, i);
    }

    // Build min-cut graph and find optimal cut
    size_t totalSamples = costSamples.numElem();
    FloatGraph graph(totalSamples,
        k * numEdges * (numSamples - maxDiff - 1) + totalSamples - k * costSamples.nx +
		(k - 1) * (1 + numNodes * (2 * numSamples - maxSep - minSep)), 
		graphErrFunc
	);
    graph.add_node(totalSamples);
    for (int i = 0; i < k; ++i) {
        buildSurfaceGraph(graph, costSamples, meshes[i], maxDiff, costType, i);
    }
    // Add intersurface edges
    for (int m = 0; m < k - 1; ++m) {
        for (int vk = 0; vk < numNodes; ++vk) {
            for (int i = maxSep; i < numSamples; ++i) {
                size_t ni = costSamples.idx(vk, i, m);
                size_t nj = costSamples.idx(vk, i - maxSep, m + 1);
                graph.add_edge(ni, nj, infOrMax<float>, 0);
            }
            for (int i = 0; i < numSamples - minSep; ++i) {
                size_t ni = costSamples.idx(vk, i, m + 1);
                size_t nj = costSamples.idx(vk, i + minSep, m);
                graph.add_edge(ni, nj, infOrMax<float>, 0);
            }
        }        
        size_t ni = costSamples.idx(0, 0, m);
        size_t nj = costSamples.idx(0, 0, m + 1);
        graph.add_edge(ni, nj, infOrMax<float>, 0);
    }

    graph.maxflow();

    for (int m = 0; m < k; ++m) {
		updateVertices(graph, costSamples, meshes[m], sampleStep, m);
    }

	return meshes;
}

ManifoldMesh& updateVertices(const FloatGraph& graph, const Volume<float>& costSamples, 
	ManifoldMesh& mesh, float sampleStep, size_t k, size_t offset)
{
	for (auto& v : mesh.vertices) {
		// Find upper position for this vertex
		for (int i = costSamples.ny - 1; i >= 0; --i) {
			size_t ni = costSamples.idx(v.self, i, k) + offset;
			if (graph.what_segment(ni) == SOURCE) {
				mesh.vpos(v) += i * sampleStep * mesh.vnormal(v);
				break;
			}
		}
	}

	return mesh;
}

ManifoldMesh& updateVerticesDual(const FloatGraph& graph, const Volume<float>& costSamples,
	ManifoldMesh& mesh, float sampleStep, size_t k, size_t offset)
{
	// TODO: Merge with updateVertices
	for (auto& v : mesh.vertices) {
		// Find upper position for this vertex
		for (int i = costSamples.ny - 1; i >= 0; --i) {
			size_t ni = costSamples.idx(v.self, i, k) + offset;
			if (graph.what_segment(ni, SINK) == SINK) {
				mesh.vpos(v) += i * sampleStep * mesh.vnormal(v);
				break;
			}
		}
	}

	return mesh;
}

ManifoldMesh& updateVerticesQPBO(const FloatGraph& graph, const Volume<float>& costSamples,
	ManifoldMesh& mesh, float sampleStep, size_t kp, size_t kd, size_t offset)
{
	// TODO: Merge with updateVertices
	for (auto& v : mesh.vertices) {
		// Find upper position for this vertex
		int ip = 0;
		for (int i = costSamples.ny - 1; i >= 0; --i) {
			size_t ni = costSamples.idx(v.self, i, kp) + offset;
			if (graph.what_segment(ni) == SOURCE) {
				ip = i;
				break;
			}
		}
		int id = 0;
		for (int i = costSamples.ny - 1; i >= 0; --i) {
			size_t ni = costSamples.idx(v.self, i, kd) + offset;
			if (graph.what_segment(ni, SINK) == SINK) {
				id = i;
				break;
			}
		}
		if (ip == id) {
			// Primal and dual graph agree so just update position
			mesh.vpos(v) += ip * sampleStep * mesh.vnormal(v);
		} else {
			throw std::runtime_error("Primal and dual graph did not agree");
		}
	}

	return mesh;
}

FloatGraph& buildSurfaceGraph(FloatGraph& graph, const Volume<float>& costSamples,
    const ManifoldMesh& mesh, int maxDiff, CostType costType, size_t k, size_t offset)
{
    using VertKey = ManifoldMesh::VertKey;
    using EdgeKey = ManifoldMesh::EdgeKey;

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
        EdgeKey ek0 = mesh.twin(v.edge);
        EdgeKey ek = ek0;
        do {
            assert(mesh.edges[ek].vert != v.self && mesh.edges[mesh.twin(ek)].vert == v.self);
            VertKey nk = mesh.edges[ek].vert;
            for (int i = numSamples - 1; i > maxDiff; --i) {
                size_t ni = costSamples.idx(v.self, i, k) + offset;
                size_t nj = costSamples.idx(nk, i - maxDiff, k) + offset;
                graph.add_edge(ni, nj, infOrMax<float>, 0);
            }
            ek = mesh.twin(mesh.next(ek)); // Next neighbor
        } while (ek != ek0);

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
            } else {
                graph.add_tweights(ni, 0, w);
            }
            prevC = c;
        }
    }

    return graph;
}

FloatGraph& buildDualSurfaceGraph(FloatGraph& graph, const Volume<float>& costSamples,
	const ManifoldMesh& mesh, int maxDiff, CostType costType, size_t k, size_t offset)
{
	// TODO: Merge with buildSurfaceGraph
	using VertKey = ManifoldMesh::VertKey;
	using EdgeKey = ManifoldMesh::EdgeKey;

	size_t numSamples = costSamples.ny;
	size_t totalSamples = costSamples.numElem();

	// Add edges
	for (const auto& v : mesh.vertices) {
		// Add intracolumn (upward) edges
		for (int i = 0; i < numSamples - 1; ++i) {
			size_t ni = costSamples.idx(v.self, i, k) + offset;
			size_t nj = costSamples.idx(v.self, i + 1, k) + offset;
			graph.add_edge(ni, nj, infOrMax<float>, 0);
		}

		// Add intercolumn (neighbor) edges
		EdgeKey ek0 = mesh.twin(v.edge);
		EdgeKey ek = ek0;
		do {
			assert(mesh.edges[ek].vert != v.self && mesh.edges[mesh.twin(ek)].vert == v.self);
			VertKey nk = mesh.edges[ek].vert;
			for (int i = 0; i < numSamples - maxDiff - 1; ++i) {
				size_t ni = costSamples.idx(v.self, i, k) + offset;
				size_t nj = costSamples.idx(nk, i + maxDiff, k) + offset;
				graph.add_edge(ni, nj, infOrMax<float>, 0);
			}
			ek = mesh.twin(mesh.next(ek)); // Next neighbor
		} while (ek != ek0);

		// Add source and sink edges
		float prevC = 0.0f;
		for (int i = 0; i < numSamples; ++i) {
			size_t ni = costSamples.idx(v.self, i, k) + offset;
			float c = -costSamples.at(v.self, i, k);
			float w = costType == ON_SURFACE ? c - prevC : c;
			if (i == 0) {
				w = -1.0f;
			}
			if (w < 0) {
				graph.add_tweights(ni, -w, 0);
			} else {
				graph.add_tweights(ni, 0, w);
			}
			prevC = c;
		}
	}

	return graph;
}

FloatGraph& buildQPBOSurfaceGraph(FloatGraph& graph, const Volume<float>& costSamples,
	const ManifoldMesh& mesh, int maxDiff, CostType costType, size_t kp, size_t kd, size_t offset)
{
	// TODO: Merge with buildSurfaceGraph
	using VertKey = ManifoldMesh::VertKey;
	using EdgeKey = ManifoldMesh::EdgeKey;

	size_t numSamples = costSamples.ny;
	size_t totalSamples = costSamples.numElem();

	// Add edges
	for (const auto& v : mesh.vertices) {
		// Add intracolumn (downward/upward) edges
		for (int i = numSamples - 1; i > 0; --i) {
			size_t nip = costSamples.idx(v.self, i, kp) + offset;
			size_t nid = costSamples.idx(v.self, i, kd) + offset;
			size_t njp = costSamples.idx(v.self, i - 1, kp) + offset;
			size_t njd = costSamples.idx(v.self, i - 1, kd) + offset;
			graph.add_edge(nip, njp, infOrMax<float>, 0);
			graph.add_edge(njd, nid, infOrMax<float>, 0);
		}

		// Add intercolumn (neighbor) edges
		EdgeKey ek0 = mesh.twin(v.edge);
		EdgeKey ek = ek0;
		do {
			assert(mesh.edges[ek].vert != v.self && mesh.edges[mesh.twin(ek)].vert == v.self);
			VertKey nk = mesh.edges[ek].vert;
			for (int i = numSamples - 1; i > maxDiff; --i) {
				size_t nip = costSamples.idx(v.self, i, kp) + offset;
				size_t nid = costSamples.idx(v.self, i, kd) + offset;
				size_t njp = costSamples.idx(nk, i - maxDiff, kp) + offset;
				size_t njd = costSamples.idx(nk, i - maxDiff, kd) + offset;
				graph.add_edge(nip, njp, infOrMax<float>, 0);
				graph.add_edge(njd, nid, infOrMax<float>, 0);
			}
			ek = mesh.twin(mesh.next(ek)); // Next neighbor
		} while (ek != ek0);

		// Add source and sink edges
		float prevC = 0.0f;
		for (int i = 0; i < numSamples; ++i) {
			size_t nip = costSamples.idx(v.self, i, kp) + offset;
			size_t nid = costSamples.idx(v.self, i, kd) + offset;
			float c = costSamples.at(v.self, i, kp);
			float w = 0.5f * (costType == ON_SURFACE ? c - prevC : c);
			if (i == 0) {
				w = -0.5f;
			}
			if (w < 0) {
				graph.add_tweights(nip, -w, 0);
				graph.add_tweights(nid, 0, -w);
			} else {
				graph.add_tweights(nip, 0, w);
				graph.add_tweights(nid, w, 0);
			}
			prevC = c;
		}
	}

	return graph;
}

Volume<float>& extractCostSamples(const Volume<float>& cost, const ManifoldMesh& mesh,
    Volume<float>& samples, int numSamples, float sampleStep, size_t k)
{
    // Sample volume along vertex normals
    for (const auto& v : mesh.vertices) {
        Vec3f p = mesh.vpos(v);
		Vec3f nrm = mesh.vnormal(v);
        for (int i = 0; i < numSamples; ++i) {
			//p[2] = fmaxf(0.0f, fminf(cost.nz - 1, p[2])); // Pad volume to repeat top and bottom slice
			samples.at(v.self, i, k) = cost.contains(p) ? cost.interp(p) : infOrMax<float>;
            p += sampleStep * nrm;
        }
    }

    return samples;
}
