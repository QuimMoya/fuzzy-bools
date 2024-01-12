#pragma once

#include <vector>
#include <glm/glm.hpp>

#include "math.h"
#include "aabb.h"

namespace fuzzybools
{
	constexpr int VERTEX_FORMAT_SIZE_FLOATS = 6;

	struct PlaneBasis
	{
		glm::dvec3 origin;

		glm::dvec3 up;
		glm::dvec3 left;
		glm::dvec3 right;

		glm::dvec2 project(const glm::dvec3& pt)
		{
			auto relative = pt - origin;
			return glm::dvec2(glm::dot(relative, left), glm::dot(relative, right));
		}
	};

	struct ReferencePlane
	{
		size_t planeID;
		size_t pointID;
		size_t lineID;
		glm::dvec2 location;
	};

	struct ReferenceLine
	{
		size_t lineID;
		size_t pointID;
		double location;
	};

	struct Line
	{
		size_t id;
		size_t globalID;
		glm::dvec3 origin;
		glm::dvec3 direction;

		Line()
		{
			static size_t idcounter = 0;
			idcounter++;
			globalID = idcounter;
		}

		/*
		bool IsPointOnLine(const glm::dvec3& pos) const
		{
			double t = glm::dot(pos - origin, direction);
			auto closestPoint = origin + t * direction;

			double dist = glm::distance(closestPoint, pos);
			return dist < EPS_BIG;
		}
		*/
		bool IsPointOnLine(const glm::dvec3& pos) const
		{
			glm::dvec3 A = pos;
			glm::dvec3 B = origin;
			glm::dvec3 C = origin + direction;

			glm::dvec3 d = (C - B) / glm::distance(C, B);
			glm::dvec3 v = A - B;
			double t = glm::dot(v, d);
			glm::dvec3 P = B + t * d;
			return glm::distance(P, A) < EPS_BIG();
		}

		double GetPosOnLine(const glm::dvec3& pos) const
		{
			return glm::dot(pos - origin, direction);
		}

		glm::dvec3 GetPosOnLine(const double dist) const
		{
			return origin + direction * dist;
		}

		bool IsColinear(const Line& other) const
		{
			return (equals(other.direction, direction, EPS_SMALL()) || equals(other.direction, -direction, EPS_SMALL()));
		}

		bool IsEqualTo(const glm::dvec3& pos, const glm::dvec3& dir) const
		{
			// check dir
			if (!(equals(dir, direction, EPS_SMALL()) || equals(dir, -direction, EPS_SMALL())))
			{
				return false;
			}

			// check pos
			if (!IsPointOnLine(pos))
			{
				return false;
			}

			return true;
		}

		void AddPointToLine(double dist, size_t id)
		{
			// check existing
			for (auto& p : points)
			{
				if (p.second == id) return;
			}

			// add new point
			points.push_back(std::make_pair(dist, id));

			// re-sort all
			std::sort(points.begin(), points.end(),
				[&](const std::pair<double, size_t>& left, const std::pair<double, size_t>& right) {
					return left.first < right.first;
				});
		}

		std::vector<std::pair<size_t, size_t>> GetSegments() const
		{
			std::vector<std::pair<size_t, size_t>> retval;

			for (size_t i = 1; i < points.size(); i++)
			{
				retval.push_back(std::make_pair(points[i - 1].second, points[i].second));
			}

			return retval;
		}

		std::vector<std::pair<double, size_t>> points;

		std::vector<ReferencePlane> planes;
	};

	struct Point
	{
		size_t id;
		size_t globalID;
		glm::dvec3 location3D;

		Point()
		{
			static size_t idcounter = 0;
			idcounter++;
			globalID = idcounter;
		}

		bool operator==(const glm::dvec3& pt)
		{
			return equals(location3D, pt, EPS_BIG());
		}

		std::vector<ReferenceLine> lines;
		std::vector<ReferencePlane> planes;
	};

	struct Plane
	{
		size_t id;
		size_t globalID;
		double distance;
		glm::dvec3 normal;

		std::vector<Line> lines;
		AABB aabb;

		void AddPoint(const glm::dvec3& pt)
		{
			aabb.merge(pt);
		}

		Plane()
		{
			static size_t idcounter = 0;
			idcounter++;
			globalID = idcounter;

			if (globalID == 320 || globalID == 321)
			{
				printf("adsf");
			}
		}

		double round(double input)
		{
			input = std::fabs(input) < EPS_BIG() ? 0.0 : input;
			input = std::fabs(input) < (1 - EPS_BIG()) ? input :
				input > 0 ? 1.0 : -1.0;
			return input;
		}

		glm::dvec3 round(glm::dvec3 in)
		{
			in.x = round(in.x);
			in.y = round(in.y);
			in.z = round(in.z);

			return in;
		}

		glm::dvec3 GetDirection(glm::dvec3 a, glm::dvec3 b)
		{
			auto dir = b - a;
			return glm::normalize(dir);
		}

		std::pair<size_t, bool> AddLine(const Point& a, const Point& b)
		{
			glm::dvec3 pos = a.location3D;
			glm::dvec3 dir = GetDirection(pos, b.location3D);

			auto lineId = AddLine(pos, dir);

			if (!lines[lineId.first].IsPointOnLine(a.location3D))
			{
				printf("bad point");
			}
			if (!lines[lineId.first].IsPointOnLine(b.location3D))
			{
				printf("bad point");
			}

			if (!aabb.contains(a.location3D))
			{
				printf("bad points");
			}
			if (!aabb.contains(b.location3D))
			{
				printf("bad points");
			}

			lines[lineId.first].AddPointToLine(lines[lineId.first].GetPosOnLine(a.location3D), a.id);
			lines[lineId.first].AddPointToLine(lines[lineId.first].GetPosOnLine(b.location3D), b.id);

			return lineId;
		}

		std::pair<size_t, bool> AddLine(const glm::dvec3& pos, const glm::dvec3& dir)
		{
			for (auto& line : lines)
			{
				if (line.IsEqualTo(pos, dir))
				{
					return { line.id, false };
				}
			}

			Line l;
			l.id = lines.size();
			l.origin = pos;
			l.direction = dir;

			lines.push_back(l);

			return { l.id, true };
		}

		void RemoveLastLine()
		{
			lines.pop_back();
		}

		bool IsEqualTo(const glm::dvec3& n, double d)
		{
			// TODO: this EPS_BIG2 is too large, 1 mm
			return (equals(normal, n, EPS_BIG2()) && equals(distance, d, EPS_BIG2())) ||
				(equals(normal, -n, EPS_BIG2()) && equals(distance, -d, EPS_BIG2()));
		}

		glm::dvec2 GetPosOnPlane(const glm::dvec3& pos)
		{
			return {};
		}

		bool HasOverlap(const std::pair<size_t, size_t>& A, const std::pair<size_t, size_t>& B)
		{
			return (A.first == B.first || A.first == B.second || A.second == B.first || A.second == B.second);
		}

		void PutPointOnLines(Point& p)
		{
			for (auto& l : lines)
			{
				if (l.IsPointOnLine(p.location3D))
				{
					ReferenceLine ref;
					ref.pointID = p.id;
					ref.lineID = l.id;
					ref.location = l.GetPosOnLine(p.location3D);
					p.lines.push_back(ref);
				}
			}
		}

		bool IsPointOnPlane(const glm::dvec3& pos)
		{
			double d = glm::dot(normal, pos);
			return equals(distance, d, EPS_BIG());
		}

		PlaneBasis MakeBasis()
		{
			glm::dvec3 origin = normal * distance;
			glm::dvec3 up = normal;

			glm::dvec3 worldUp = glm::dvec3(0, 1, 0);
			glm::dvec3 worldRight = glm::dvec3(1, 0, 0);

			bool normalIsUp = equals(up, worldUp, EPS_SMALL()) || equals(-up, worldUp, EPS_SMALL());
			glm::dvec3 left = normalIsUp ? glm::cross(up, worldRight) : glm::cross(up, worldUp);
			glm::dvec3 right = glm::cross(left, up);

			PlaneBasis basis;

			basis.origin = origin;
			basis.up = glm::normalize(up);
			basis.left = glm::normalize(left);
			basis.right = glm::normalize(right);

			return basis;
		}
	};

	struct Face
	{
		int i0;
		int i1;
		int i2;
		int ip;
	};

	struct Geometry
	{
		std::vector<float> fvertexData;
		std::vector<double> vertexData;
		std::vector<uint32_t> indexData;
		std::vector<uint32_t> indexPlane;
		std::vector<Plane> planesData;


		void BuildFromVectors(std::vector<double>& d, std::vector<uint32_t>& i)
		{
			vertexData = d;
			indexData = i;

			numPoints = indexData.size();
			numFaces = indexData.size() / 3;
		}

		uint32_t numPoints = 0;
		uint32_t numFaces = 0;

		inline void AddPoint(glm::dvec4& pt, glm::dvec3& n)
		{
			glm::dvec3 p = pt;
			AddPoint(p, n);
		}

		AABB GetAABB() const
		{
			AABB aabb;

			for (uint32_t i = 0; i < numPoints; i++)
			{
				aabb.min = glm::min(aabb.min, GetPoint(i));
				aabb.max = glm::max(aabb.max, GetPoint(i));
			}

			return aabb;
		}

		inline void AddPoint(glm::dvec3& pt, glm::dvec3& n)
		{
			//vertexData.reserve((numPoints + 1) * VERTEX_FORMAT_SIZE_FLOATS);
			//vertexData[numPoints * VERTEX_FORMAT_SIZE_FLOATS + 0] = pt.x;
			//vertexData[numPoints * VERTEX_FORMAT_SIZE_FLOATS + 1] = pt.y;
			//vertexData[numPoints * VERTEX_FORMAT_SIZE_FLOATS + 2] = pt.z;
			vertexData.push_back(pt.x);
			vertexData.push_back(pt.y);
			vertexData.push_back(pt.z);

			vertexData.push_back(n.x);
			vertexData.push_back(n.y);
			vertexData.push_back(n.z);

			if (std::isnan(pt.x) || std::isnan(pt.y) || std::isnan(pt.z))
			{
				printf("NaN in geom!\n");
			}

			if (std::isnan(n.x) || std::isnan(n.y) || std::isnan(n.z))
			{
				printf("NaN in geom!\n");
			}

			//vertexData[numPoints * VERTEX_FORMAT_SIZE_FLOATS + 3] = n.x;
			//vertexData[numPoints * VERTEX_FORMAT_SIZE_FLOATS + 4] = n.y;
			//vertexData[numPoints * VERTEX_FORMAT_SIZE_FLOATS + 5] = n.z;

			numPoints += 1;
		}

		inline void AddFace(glm::dvec3 a, glm::dvec3 b, glm::dvec3 c, uint32_t p)
		{
			glm::dvec3 normal;

			double area = areaOfTriangle(a, b, c);

			if (!computeSafeNormal(a, b, c, normal, EPS_SMALL()))
			{
				// bail out, zero area triangle
				printf("zero tri");
				return;
			}

			AddPoint(a, normal);
			AddPoint(b, normal);
			AddPoint(c, normal);

			AddFace(numPoints - 3, numPoints - 2, numPoints - 1, p);
		}

		inline void AddFace(uint32_t a, uint32_t b, uint32_t c, uint_32 p)
		{
			//indexData.reserve((numFaces + 1) * 3);
			//indexData[numFaces * 3 + 0] = a;
			//indexData[numFaces * 3 + 1] = b;
			//indexData[numFaces * 3 + 2] = c;
			indexData.push_back(a);
			indexData.push_back(b);
			indexData.push_back(c);
			facePlane.push_back(p);

			double area = areaOfTriangle(GetPoint(a), GetPoint(b), GetPoint(c));

			glm::dvec3 normal;
			if (!computeSafeNormal(GetPoint(a), GetPoint(b), GetPoint(c), normal))
			{
				// bail out, zero area triangle
				printf("zero tri");
			}

			numFaces++;
		}

		inline Face GetFace(size_t index) const
		{
			Face f;
			f.i0 = indexData[index * 3 + 0];
			f.i1 = indexData[index * 3 + 1];
			f.i2 = indexData[index * 3 + 2];
			f.ip = indexPlane[index];
			return f;
		}

		inline AABB GetFaceBox(size_t index) const
		{
			AABB aabb;
			aabb.index = static_cast<uint32_t>(index);

			glm::dvec3 a = GetPoint(indexData[index * 3 + 0]);
			glm::dvec3 b = GetPoint(indexData[index * 3 + 1]);
			glm::dvec3 c = GetPoint(indexData[index * 3 + 2]);

			aabb.min = glm::min(a, aabb.min);
			aabb.min = glm::min(b, aabb.min);
			aabb.min = glm::min(c, aabb.min);

			aabb.max = glm::max(a, aabb.max);
			aabb.max = glm::max(b, aabb.max);
			aabb.max = glm::max(c, aabb.max);

			aabb.center = (aabb.max + aabb.min) / 2.0;

			return aabb;
		}

		inline glm::dvec3 GetPoint(size_t index) const
		{
			return glm::dvec3(
				vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 0],
				vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 1],
				vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 2]
			);
		}

		void GetCenterExtents(glm::dvec3& center, glm::dvec3& extents) const
		{
			glm::dvec3 min(DBL_MAX, DBL_MAX, DBL_MAX);
			glm::dvec3 max(DBL_MIN, DBL_MIN, DBL_MIN);

			for (size_t i = 0; i < numPoints; i++)
			{
				auto pt = GetPoint(i);
				min = glm::min(min, pt);
				max = glm::max(max, pt);
			}

			extents = (max - min);
			center = min + extents / 2.0;
		}

		Geometry Normalize(glm::dvec3 center, glm::dvec3 extents) const
		{
			Geometry newGeom;

			double scale = std::max(extents.x, std::max(extents.y, extents.z)) / 10.0;

			for (size_t i = 0; i < numFaces; i++)
			{
				auto face = GetFace(i);
				auto pa = GetPoint(face.i0);
				auto pb = GetPoint(face.i1);
				auto pc = GetPoint(face.i2);

				auto a = (pa - center) / scale;
				auto b = (pb - center) / scale;
				auto c = (pc - center) / scale;

				// std::cout << areaOfTriangle(pa, pb, pc) << std::endl;

				newGeom.AddFace(pa, pb, pc, face.ip);
			}

			return newGeom;
		}

		size_t AddPlane(const glm::dvec3& normal, double d)
		{
			for (auto& plane : planesData)
			{
				if (plane.IsEqualTo(normal, d))
				{
					return plane.id;
				}
			}

			Plane p;
			p.id = planesData.size();
			p.normal = normal;
			p.distance = d;
			planesData.push_back(p);

			return p.id;
		}

		void GeneratePlanes()
		{
			if(planesData.size() == 0)
			{
				for (size_t i = 0; i < numFaces; i++)
				{
					Face f = GetFace(i);

					auto faceBox = GetFaceBox(i);

					auto a = GetPoint(f.i0);
					auto b = GetPoint(f.i1);
					auto c = GetPoint(f.i2);

					glm::dvec3 norm;
					if (computeSafeNormal(a, b, c, norm, EPS_SMALL()))
					{
						double da = glm::dot(norm, a);
						double db = glm::dot(norm, b);
						double dc = glm::dot(norm, c);

						size_t planeId = AddPlane(norm, da);
						f.ip = planeId;
					}
				}
			}
		}

		Geometry DeNormalize(glm::dvec3 center, glm::dvec3 extents) const
		{
			Geometry newGeom;

			double scale = std::max(extents.x, std::max(extents.y, extents.z)) / 10.0;

			for (size_t i = 0; i < numFaces; i++)
			{
				auto face = GetFace(i);
				auto pa = GetPoint(face.i0);
				auto pb = GetPoint(face.i1);
				auto pc = GetPoint(face.i2);

				// std::cout << areaOfTriangle(pa, pb, pc) << std::endl;

				auto a = pa * scale + center;
				auto b = pb * scale + center;
				auto c = pc * scale + center;

				newGeom.AddFace(a, b, c, face.ip);
			}


			return newGeom;
		}
		
		bool IsEmpty()
		{
			return vertexData.empty();
		}

		double Volume(const glm::dmat4& trans = glm::dmat4(1))
		{
			double totalVolume = 0;

			for (uint32_t i = 0; i < numFaces; i++)
			{
				Face f = GetFace(i);

				glm::dvec3 a = trans * glm::dvec4(GetPoint(f.i0), 1);
				glm::dvec3 b = trans * glm::dvec4(GetPoint(f.i1), 1);
				glm::dvec3 c = trans * glm::dvec4(GetPoint(f.i2), 1);

				glm::dvec3 norm;

				if (computeSafeNormal(a, b, c, norm))
				{
					double area = areaOfTriangle(a, b, c);
					double height = glm::dot(norm, a);

					double tetraVolume = area * height / 3;

					totalVolume += tetraVolume;
				}
			}

			return totalVolume;
		}
	};
}