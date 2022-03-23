// PrimPrimContactInfo.cpp

#include "PrimPrimContactInfo.h"

namespace Sai2COPSim {

void PrimPrimContactInfo::flipNormal() {
	normal_dir *= -1.0;
	if(abs(primB_max_radius - primB_min_radius) > 1e-5) {
		// align constraint dir1 and constraint dir2 with primB max and min curvature
		// planes
		double ca = cos(inter_prim_max_curvature_plane_angle);
		double sa = sin(inter_prim_max_curvature_plane_angle);
		Eigen::Matrix3d rotz;
		rotz << ca, -sa,  0,
				sa,  ca,  0,
				0,	  0,  1;
		constraint_dir1 = rotz * constraint_dir1;
		constraint_dir2 = rotz * constraint_dir2;
	}
	constraint_dir2 *= -1.0;
	std::swap(primA_max_radius, primB_max_radius);
	std::swap(primA_min_radius, primB_min_radius);
}

}