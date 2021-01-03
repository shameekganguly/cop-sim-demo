// ContactPatch.h

#ifndef CONTACT_PATCH_H
#define CONTACT_PATCH_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>

namespace Sai2COPSim {

namespace ContactPatchAlgorithmicConstants {
	const double POINT_ON_BOUNDARY_DISTANCE_THRESHOLD = 1e-4;
};

// NOTE: line points should be in patch coordinates. patch interior point is at (0,0)
// in patch coordinates
struct LineSegment {
	Eigen::Vector2d point1;
	Eigen::Vector2d point2;
	double length;
	Eigen::Vector2d dir; // from point 1 to point 2
	Eigen::Vector2d interior_pt_proj;

	// ctor
	LineSegment(const Eigen::Vector2d& point1, const Eigen::Vector2d& point2)
	: point1(point1), point2(point2)
	{
		Eigen::Vector2d temp = point2 - point1;
		length = temp.norm();
		assert(length > 1e-1); // TODO: this should be relative to some global object feature
		dir = temp/length;
		interior_pt_proj = point1 + (-point1.dot(dir))*dir;
	}

	// note that the line segment is treated as a line for this.
	// i.e. min distance from the line itself is considered
	double distanceToPoint(const Eigen::Vector2d& point) const {
		Eigen::Vector2d temp = point - point1;
		Eigen::Vector2d temp2 = temp - (temp.dot(dir))*dir;
		return (temp.dot(interior_pt_proj) < 0)? temp2.norm() : -temp2.norm();
	}

	// assumes that the test point in inside the contact patch
	// note that the line segment is treated as a line for this.
	// NOTE: returned distance is positive if curve intersects in the direction of the ray.
	// if it intersects in the opposite direction, returned distance is negative.
	double distanceFromPointAlongRay(const Eigen::Vector2d& point, const Eigen::Vector2d& direction) const {
		Eigen::Vector2d temp = point - point1;
		Eigen::Vector2d temp2 = temp - (temp.dot(dir))*dir;
		double temp3 = direction.dot(dir);
		double dist_along_ray = temp2.norm()/sqrt(1 - temp3*temp3); // will be inf when the two directions are parallel
		return (direction.dot(temp2) < 0)? dist_along_ray: -dist_along_ray;
	}
};

//TODO: these should be only convex curves in the sense that the interior point
// should be on the convex side of the curve
struct Curve {
	// distance to closest point on the curve from the test point
	// Note: distance is positive if point is on the same side of the curve as the
	// interior point, else it is negative
	virtual double distanceToPoint(const Eigen::Vector2d& point) const {
		return 0;
	}

	// this function assumes that the point is INSIDE the curve, that is, it is on the
	// same side of the curve as the interior point
	// NOTE: returned distance is positive if curve intersects in the direction of the ray.
	// if it intersects in the opposite direction, returned distance is negative.
	virtual double distanceFromPointAlongRay(const Eigen::Vector2d& point, const Eigen::Vector2d& direction) const {
		return 0;
	}
};

struct Circle: public Curve {
	Eigen::Vector2d center;
	double distanceToPoint(const Eigen::Vector2d& point) const {
		return radius - (point - center).norm();
	}
	double distanceFromPointAlongRay(const Eigen::Vector2d& point, const Eigen::Vector2d& direction) const {
		double x_delta_times_cos_alp = (point - center).dot(direction);
		double x_delta = (point - center).norm();
		return sqrt(radius*radius
						+ x_delta_times_cos_alp*x_delta_times_cos_alp
						- x_delta*x_delta
		) - x_delta_times_cos_alp;
	}
	double radius;
};

struct PointTestResult {
	bool f_is_in_patch;
	bool f_is_on_line_seg;
	uint line_seg_id;
	bool f_is_on_vertex;
	double min_dist_to_boundary;
	double max_cross_extent_at_point;
	PointTestResult()
	: f_is_in_patch(false), f_is_on_line_seg(false), line_seg_id(0), f_is_on_vertex(false), min_dist_to_boundary(0.0) {
		// nothing to do
	}
};

class ContactPatch {
public:
	// ctor

	// dtor
	~ContactPatch() {
		clear();
	}

	Eigen::Vector3d patchToLineSegmentCoordinates(uint line_segment_id);
	//TODO: we need the rotation matrix as well unless we use arbitrary line segment directions
	// in the line COP solver

	// test where given point is with respect to the boundary
	PointTestResult testPoint(const Eigen::Vector2d& point) const;

	// get distance to boundary in a given ray direction
	// this function assumes that the point is INSIDE the patch, that is, it is on the
	// same side of the patch boundary as the interior point
	// returned distance is always positive for a well defined contact patch
	double distanceFromBoundaryAlongRay(const Eigen::Vector2d& point, const Eigen::Vector2d& direction) const;

	// clear the data
	void clear() {
		// clear line segments
		_line_segments.clear();

		// deallocate curves memory
		for(auto c: _intersection_curves) {
			if(c != NULL) delete c;
			c = NULL;
		}

		// clear curves
		_intersection_curves.clear();

		// clear max extent
		max_extent = -1;
	}

public:
	double max_extent; // used to compute mu_rotation
	std::vector<Curve*> _intersection_curves; //TODO: what about unions? think about a box with rounded corners
	Eigen::Vector3d _interior_point; // origin of patch coordinates in global frame
	std::vector<LineSegment> _line_segments;
};

}

#endif //CONTACT_PATCH_H