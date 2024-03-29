// PrimPrimContactInfo.h

#ifndef PRIM_PRIM_CONTACT_INFO_H
#define PRIM_PRIM_CONTACT_INFO_H

#include <Eigen/Dense>
#include <vector>
#include "ContactPatch.h"

namespace Sai2COPSim {

namespace PrimitiveAlgorithmicConstants {
static constexpr double DISTANCE_FULL_PENETRATION = -1.0;
}

enum ContactType {
	UNDEFINED,
	POINT,
	LINE,
	SURFACE,
	// Currently, CONCAVE implies that the contact consists of a set of points, each with
	// a different normal, curvatures and tangent directions
	CONCAVE
};

// TODO: reorganize this to reflect the core assumption that two primitives can be in
// contact across a set of contact patches, where each contact patch is *planar*.
class PrimPrimContactInfo {
public:
	double min_distance;
	Eigen::Vector3d normal_dir; // in world frame.
	// ^ Note that this is directed from prim A to prim B, in the order that the
	// prim-prim distance computation function was called
	Eigen::Vector3d constraint_dir1; // in world frame. for line contacts and surface contacts
	Eigen::Vector3d constraint_dir2; // in world frame. for surface contacts only
	// TODO: for geometries with different min and max curvatures at the contact point, we need
	// to ensure that constraint_dir1 and constraint_dir2 are aligned with the max and min
	// curvature planes respectively of primA.
	// TODO: for concave objects, the each contact point can have its own normal and
	// constraint dirs. Consider if we want to support that, or assume that each primitive
	// is strictly convex
	std::vector<Eigen::Vector3d> contact_points; // closest points on either prim A or prim B in world frame
	//TODO: do we need to be consistent about which primitive the points lie on?
	// or do we need to return points on both bodies?

	// members used for CONCAVE contact type
	std::vector<Eigen::Vector3d> normal_dirs; // in world frame, corresponding to each contact point
	std::vector<Eigen::Vector3d> constraint_dir1s; // in world frame, corresponding to each contact point
	std::vector<Eigen::Vector3d> constraint_dir2s; // in world frame, corresponding to each contact point
	std::vector<double> distances; // distance at each contact point

	// signed curvature for each body
	// sign is positive is center of curvature lies in the positive normal direction
	// and negative otherwise
	// TODO: this works only for point contacts. need to extend to line contacts and surface contacts
	double primA_max_radius;
	double primA_min_radius;
	double primB_max_radius;
	double primB_min_radius;

	// angle from max curvature plane of primA to max curvature plane of primB
	double inter_prim_max_curvature_plane_angle;

	// Note: contact_points might not be set for a surface-surface contact. e.g. cylinder
	// on plane
	// contact patch info
	ContactPatch contact_patch;

	ContactType type;

public:
	PrimPrimContactInfo():
		min_distance(0),
		primA_max_radius(0),
		primA_min_radius(0),
		primB_max_radius(0),
		primB_min_radius(0),
		inter_prim_max_curvature_plane_angle(0),
		type(ContactType::UNDEFINED) { }
	PrimPrimContactInfo(double adist, ContactType atype):
		min_distance(adist),
		primA_max_radius(0),
		primA_min_radius(0),
		primB_max_radius(0),
		primB_min_radius(0),
		inter_prim_max_curvature_plane_angle(0),
		type(atype) { }

	// clear
	void clear() {
		contact_points.clear();
		type = ContactType::UNDEFINED;
		contact_patch.clear();
		primA_max_radius = 0;
		primA_min_radius = 0;
		primB_max_radius = 0;
		primB_min_radius = 0;
		inter_prim_max_curvature_plane_angle = 0;
		normal_dirs.clear();
		constraint_dir1s.clear();
		constraint_dir2s.clear();
		distances.clear();
	}

	// flip normal
	void flipNormal();

	// filter contact points, normal_dirs, constraint_dir1s, constraint_dir2s and
	// distances by a max distance threshold.
	// Only for CONCAVE contact type
	void filterContactPoints(double max_distance);

	// filter contact points, normal_dirs, constraint_dir1s, constraint_dir2s and
	// distances by a boolean condition on the contact point
	// Only for CONCAVE contact type
	// If qualifier function returns true for the test point, then that point is retained
	// otherwise the point is dropped
	void filterContactPoints(std::function<bool(const Eigen::Vector3d& pt_in_world)>& qualifier);

	void filterContactPoints(const std::vector<uint>& keep_indices);

	// Only for CONCAVE contact type. distances should be set already
	void setMinDistanceFromDistances();

	// Add other PrimPrimContactInfo and set concave
	void addOtherContactInfo(const PrimPrimContactInfo& other_contact_info);
};

}

#endif // PRIM_PRIM_CONTACT_INFO_H