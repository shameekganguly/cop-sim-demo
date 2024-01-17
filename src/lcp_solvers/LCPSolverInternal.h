#ifndef LCPSOLVER_INTERNAL_H
#define LCPSOLVER_INTERNAL_H

#include "LCPSolver.h"

#define LCP_LOG_DEBUG false

namespace Sai2LCPSolver {

CollLCPPointSolution solveCollLCPPoint(
	uint Npoints,
	const Eigen::MatrixXd& A,
	const Eigen::VectorXd& b,
	const Eigen::VectorXd& pre_v,
	double epsilon,
	double mu,
	bool is_redundant_x = true
);

CollLCPPointSolution solveCollLCPOnePoint(
	const Eigen::Matrix3d& A,
	const Eigen::Vector3d& b,
	const Eigen::Vector3d& pre_v,
	double epsilon,
	double mu
);

}  // namespace Sai2LCPSolver

#endif  // LCPSOLVER_INTERNAL_H
