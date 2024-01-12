#pragma once

#include "geometry.h"
#include "shared-position.h"
#include "clip-mesh.h"

namespace fuzzybools
{

	inline void CopyPlanes(std::vector<Plane>& plane_list, Geometry geom)
	{
		for (auto& plane : geom.planesData)
		{
			for (auto& plane2 : plane_list)
			{
				if (plane2.IsEqualTo(plane.normal, plane.distance))
				{
					return;
				}
			}
			
			Plane p;
			p.id = plane_list.size();
			p.normal = plane.normal;
			p.distance = plane.distance;
			plane_list.push_back(p);
		}
	}

	inline Geometry Subtract(const Geometry& A, const Geometry& B, double scale = 1)
	{
		EPS_SCALE = scale;
		fuzzybools::SharedPosition sp;
		sp.Construct(A, B);

		std::vector<Plane> plane_list;
		CopyPlanes(plane_list, A);
		CopyPlanes(plane_list, B);

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);

		auto geom = Normalize(sp, plane_list);

		DumpGeometry(geom, L"Post-normalize.obj");

		Geometry result = fuzzybools::clipSubtract(geom, bvh1, bvh2);

		result.planesData = plane_list;

		return result;
	}

	inline Geometry Union(const Geometry& A, const Geometry& B, double scale = 1)
	{
		EPS_SCALE = scale;
		fuzzybools::SharedPosition sp;
		sp.Construct(A, B);

		std::vector<Plane> plane_list;
		CopyPlanes(plane_list, A);
		CopyPlanes(plane_list, B);

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);

		auto geom = Normalize(sp, plane_list);

		return fuzzybools::clipJoin(geom, bvh1, bvh2);
	}
}