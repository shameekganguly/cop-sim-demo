/* ContactCurvature.cpp */

#include "ContactSpaceModel.h"

namespace Sai2COPSim {

double ContactPairState::getAccelerationDueToCurvedContact(
	const Eigen::Vector3d& omegaA, // in contact patch frame
	const Eigen::Vector3d& omegaB, // in contact patch frame
	const Eigen::Vector3d& linear_contact_velocity // in contact patch frame
) const {
	double ret_curve_acc = 0.0;
	if(_active_points.size() == 1) {
		const auto& prim_info = *_geom_prim_pair->info;
		if(prim_info.type == ContactType::POINT) {
			if(abs(prim_info.primA_max_radius) > 1e-5 && abs(prim_info.primA_max_radius - prim_info.primA_min_radius) < 1e-5) {
				return (omegaA.segment<2>(0) - omegaB.segment<2>(0)).squaredNorm()*prim_info.primA_max_radius -
						2.0*(omegaB.cross(linear_contact_velocity))[2];
			} //else if (abs(prim_info.primB_max_radius) > 1e-5 && abs(prim_info.primB_max_radius - prim_info.primB_min_radius) < 1e-5)
		}
	}
	return ret_curve_acc;
}

}