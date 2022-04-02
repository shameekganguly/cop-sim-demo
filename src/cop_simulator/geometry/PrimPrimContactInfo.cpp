// PrimPrimContactInfo.cpp

#include "PrimPrimContactInfo.h"

namespace Sai2COPSim {

using namespace Eigen;

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
	std::vector<uint> filtered_inds;
	for(uint i = 0; i < distances.size(); i++) {
		if(distances[i] <= max_distance) {
			filtered_inds.push_back(i);
		}
	}
	filterContactPoints(filtered_inds);
}

void PrimPrimContactInfo::filterContactPoints(
	std::function<bool(const Eigen::Vector3d& pt_in_world)> &qualifier) {
	if(type != ContactType::CONCAVE) {
		return;
	}
	std::vector<uint> filtered_inds;
	for(uint i = 0; i < contact_points.size(); i++) {
		if(qualifier(contact_points[i])) {
			filtered_inds.push_back(i);
		}
	}
	filterContactPoints(filtered_inds);
}

void PrimPrimContactInfo::filterContactPoints(const std::vector<uint>& keep_indices) {
	// TODO: this is inefficient. change the data structures (e.g. to a priority list)
	// to avoid all this copying.
	auto resize_lambda = [&keep_indices](std::vector<Vector3d>& input) {
		for(uint j = 0; j < keep_indices.size(); j++) {
			input[j] = input[keep_indices[j]];
		}
		input.resize(keep_indices.size());
	};
	resize_lambda(contact_points);
	resize_lambda(normal_dirs);
	resize_lambda(constraint_dir1s);
	resize_lambda(constraint_dir2s);
	for(uint j = 0; j < keep_indices.size(); j++) {
		distances[j] = distances[keep_indices[j]];
	}
	distances.resize(keep_indices.size());
	if(!distances.empty()) {
		setMinDistanceFromDistances();
	}
}

void PrimPrimContactInfo::setMinDistanceFromDistances() {
	min_distance = *(std::min_element(distances.begin(), distances.end()));
}

void PrimPrimContactInfo::addOtherContactInfo(const PrimPrimContactInfo& other_contact_info) {
	if(other_contact_info.type != ContactType::CONCAVE) {
		throw(std::runtime_error("Unimplemented addOtherContactInfo case: " + std::to_string(other_contact_info.type)));
	}
	if(type != ContactType::CONCAVE) {
		uint pts_size = contact_points.size();
		normal_dirs.resize(pts_size, normal_dir);
		constraint_dir1s.resize(pts_size, constraint_dir1);
		constraint_dir2s.resize(pts_size, constraint_dir2);
		distances.resize(pts_size, min_distance);
	}

	type = ContactType::CONCAVE;

	auto append_vectorB_to_vectorA_double = [](std::vector<double>& vectorA, const std::vector<double>& vectorB) {
		vectorA.insert(vectorA.end(), vectorB.begin(), vectorB.end());
	};

	auto append_vectorB_to_vectorA_Vector3d = [](std::vector<Vector3d>& vectorA, const std::vector<Vector3d>& vectorB) {
		vectorA.insert(vectorA.end(), vectorB.begin(), vectorB.end());
	};

	append_vectorB_to_vectorA_Vector3d(contact_points, other_contact_info.contact_points);
	append_vectorB_to_vectorA_Vector3d(normal_dirs, other_contact_info.normal_dirs);
	append_vectorB_to_vectorA_Vector3d(constraint_dir1s, other_contact_info.constraint_dir1s);
	append_vectorB_to_vectorA_Vector3d(constraint_dir2s, other_contact_info.constraint_dir2s);
	append_vectorB_to_vectorA_double(distances, other_contact_info.distances);

	setMinDistanceFromDistances();
}

}