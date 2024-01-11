#pragma once

#include "geometry.h"
#include "shared-position.h"
#include "clip-mesh.h"

namespace fuzzybools
{
	inline Geometry Subtract(Geometry& A, Geometry& B)
	{
		A.GeneratePlanes();
		B.GeneratePlanes();

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

	inline Geometry Union(Geometry& A, Geometry& B)
	{
		A.GeneratePlanes();
		B.GeneratePlanes();

		fuzzybools::SharedPosition sp;
		sp.Construct(A, B);

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);

		auto geom = Normalize(sp);

		return fuzzybools::clipJoin(geom, bvh1, bvh2);
	}
}