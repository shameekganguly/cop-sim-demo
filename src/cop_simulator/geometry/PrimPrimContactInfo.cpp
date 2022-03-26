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

	if (type == ContactType::CONCAVE) {
		for(auto& tnormal: normal_dirs) {
			tnormal *= -1.0;
		}
		for(auto& tconstraint_dir2: constraint_dir2s) {
			tconstraint_dir2 *= -1.0;
		}
	}
}

void PrimPrimContactInfo::filterContactPoints(double max_distance) {
	if(type != ContactType::CONCAVE) {
		return;
	}
	// TODO: this is inefficient. change the data structures (e.g. to a priority list)
	// to avoid all this copying.
	std::vector<uint> filtered_inds;
	for(uint i = 0; i < distances.size(); i++) {
		if(distances[i] <= max_distance) {
			filtered_inds.push_back(i);
		}
	}
	auto resize_lambda = [&filtered_inds](std::vector<Eigen::Vector3d>& input) {
		for(uint j = 0; j < filtered_inds.size(); j++) {
			input[j] = input[filtered_inds[j]];
		}
		input.resize(filtered_inds.size());
	};
	resize_lambda(contact_points);
	resize_lambda(normal_dirs);
	resize_lambda(constraint_dir1s);
	resize_lambda(constraint_dir2s);
	for(uint j = 0; j < filtered_inds.size(); j++) {
		distances[j] = distances[filtered_inds[j]];
	}
	distances.resize(filtered_inds.size());
}

void PrimPrimContactInfo::setMinDistanceFromDistances() {
	min_distance = *(std::min_element(distances.begin(), distances.end()));
}

}