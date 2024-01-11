#pragma once

#include "geometry.h"
#include "shared-position.h"
#include "clip-mesh.h"

namespace fuzzybools
{
	inline Geometry Subtract(const Geometry& A, const Geometry& B, double scale = 1)
	{
		EPS_SCALE = scale;
		fuzzybools::SharedPosition sp;
		sp.Construct(A, B);

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);

		auto geom = Normalize(sp);

		DumpGeometry(geom, L"Post-normalize.obj");

		Geometry result = fuzzybools::clipSubtract(geom, bvh1, bvh2);

		result.CopyPlanes(A);
		result.CopyPlanes(B);

		return result;
	}

	inline Geometry Union(const Geometry& A, const Geometry& B, double scale = 1)
	{
		EPS_SCALE = scale;
		fuzzybools::SharedPosition sp;
		sp.Construct(A, B);

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);

		auto geom = Normalize(sp);

		return fuzzybools::clipJoin(geom, bvh1, bvh2);
	}
}