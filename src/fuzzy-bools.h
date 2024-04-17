#pragma once

#include "geometry.h"
#include "shared-position.h"
#include "clip-mesh.h"

namespace fuzzybools
{
	std::vector<Geometry> GetMeshesByNeighbourTriangles(const Geometry &mesh)
	{
		fuzzybools::SharedPosition sp;
		for (size_t i = 0; i < mesh.numFaces; i++)
		{
			Face tri = mesh.GetFace(i);
			auto a = mesh.GetPoint(tri.i0);
			auto b = mesh.GetPoint(tri.i1);
			auto c = mesh.GetPoint(tri.i2);

			glm::dvec3 norm;
			if (computeSafeNormal(a, b, c, norm, EPS_SMALL))
			{
				auto ia = sp.AddPoint(a);
				auto ib = sp.AddPoint(b);
				auto ic = sp.AddPoint(c);

				double da = glm::dot(norm, a);
				double db = glm::dot(norm, b);
				double dc = glm::dot(norm, c);

				size_t planeId = sp.AddPlane(norm, da);

				sp.planes[planeId].AddPoint(a);
				sp.planes[planeId].AddPoint(b);
				sp.planes[planeId].AddPoint(c);

				sp.A.AddFace(planeId, ia, ib, ic);
			}
		}

		std::vector<Geometry> newMeshes;

		std::vector<Triangle> meshTriangles = sp.A.triangles;
		std::vector<bool> visited(meshTriangles.size(), false);
		for (size_t triangleId = 0; triangleId < meshTriangles.size(); triangleId++)
		{
			if (visited[triangleId])
				continue;

			Geometry newMesh;
			std::queue<size_t> q;
			q.push(triangleId);

			while (!q.empty())
			{
				size_t currentId = q.front();
				q.pop();
				if (visited[currentId])
					continue;

				visited[currentId] = true;

				Triangle triangle = meshTriangles[currentId];
				glm::dvec3 a = sp.points[triangle.a].location3D;
				glm::dvec3 b = sp.points[triangle.b].location3D;
				glm::dvec3 c = sp.points[triangle.c].location3D;
				newMesh.AddFace(a, b, c);

				for (const auto &item : sp.A.GetNeighbourTriangles(triangle))
				{
					for (int i = 0; i < item.second.size(); i++)
					{
						size_t neighbourId = item.second[i];
						if (visited[neighbourId])
							continue;

						q.push(neighbourId);
					}
				}
			}

			newMeshes.push_back(newMesh);
		}

		return newMeshes;
	}

	inline Geometry Subtract(const Geometry &A, const Geometry &B)
	{
		fuzzybools::SharedPosition sp;
		sp.Construct(A, B);

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);

		auto geom = Normalize(sp);

		DumpGeometry(geom, L"Post-normalize.obj");

		return fuzzybools::clipSubtract(geom, bvh1, bvh2);
	}

	inline Geometry Union(const Geometry &A, const Geometry &B)
	{
		fuzzybools::SharedPosition sp;
		sp.Construct(A, B);

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);

		auto geom = Normalize(sp);

		return fuzzybools::clipJoin(geom, bvh1, bvh2);
	}
}