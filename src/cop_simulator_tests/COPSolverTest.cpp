// COPSolverTests.cpp

#include <iostream>
#include <Eigen/Dense>
#include "cop_simulator/COPSolverExtended.h"

using namespace std;
using namespace Eigen;
using namespace Sai2COPSim;

void testLinePtCOP();

int main (int argc, char** argv) {
	testLinePtCOP();
	return 0;
}

void printCOPSolution(ContactCOPSolution& sol) {
	cout << "Result: " << uint(sol.result) << "\n";
	cout << "Forces: " << sol.force_sol.transpose() << "\n";
	cout << "Local COP pos: " << sol.local_cop_pos.transpose() << "\n";
	cout << "COP type: " << uint(sol.cop_type) << "\n";
	cout << endl;
}

void testLinePtCOP() {
	{
		cout << "--- TEST 1 ---" << endl;
		cout << "No omegas" << endl;
		cout << "No lin vels, rolling/impending" << endl;
		MatrixXd A_constraint = MatrixXd::Identity(8, 8);
		MatrixXd rhs_constraint = VectorXd::Zero(8);
		rhs_constraint(2) = -0.1;
		rhs_constraint(7) = -0.1;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d::Zero());
		omegaBs.push_back(Vector3d::Zero());
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d::Zero());
		linear_contact_vels.push_back(Vector3d::Zero());
		COPSolver solver;
		vector<vector<Vector3d>> boundary_points;
		vector<Vector3d> line_bdry_pts;
		line_bdry_pts.push_back(Vector3d::Zero());
		line_bdry_pts.push_back(Vector3d(0.2, 0, 0));
		boundary_points.push_back(line_bdry_pts);
		vector<Vector3d> pt_bdry_pts;
		boundary_points.push_back(pt_bdry_pts);
		vector<uint> patch_indices;
		patch_indices.push_back(0); // first line contact
		patch_indices.push_back(5); // then point contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::LINE);
		contact_types.push_back(ContactType::POINT);
		ContactCOPSolution sol = solver.solveStartWithPatchCentroidOnePointOneLine(
			0.1,
			A_constraint,
			rhs_constraint,
			boundary_points,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 2 ---" << endl;
		cout << "No omegas" << endl;
		cout << "No lin vels, rolling/impending" << endl;
		MatrixXd A_constraint = MatrixXd::Identity(8, 8);
		MatrixXd rhs_constraint = VectorXd::Zero(8);
		rhs_constraint(2) = -0.1;
		rhs_constraint(5) = -0.1;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d::Zero());
		omegaBs.push_back(Vector3d::Zero());
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d::Zero());
		linear_contact_vels.push_back(Vector3d::Zero());
		COPSolver solver;
		vector<vector<Vector3d>> boundary_points;
		vector<Vector3d> pt_bdry_pts;
		boundary_points.push_back(pt_bdry_pts);
		vector<Vector3d> line_bdry_pts;
		line_bdry_pts.push_back(Vector3d::Zero());
		line_bdry_pts.push_back(Vector3d(0.2, 0, 0));
		boundary_points.push_back(line_bdry_pts);
		vector<uint> patch_indices;
		patch_indices.push_back(0); // first point contact
		patch_indices.push_back(3); // then line contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::POINT);
		contact_types.push_back(ContactType::LINE);
		ContactCOPSolution sol = solver.solveStartWithPatchCentroidOnePointOneLine(
			0.1,
			A_constraint,
			rhs_constraint,
			boundary_points,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 3 ---" << endl;
		cout << "No omegas" << endl;
		cout << "Line contact sliding: lin vel, No point contact lin vel: rolling, " << endl;
		MatrixXd A_constraint = MatrixXd::Identity(8, 8);
		MatrixXd rhs_constraint = VectorXd::Zero(8);
		rhs_constraint(2) = -0.1;
		rhs_constraint(7) = -0.1;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d::Zero());
		omegaBs.push_back(Vector3d::Zero());
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d(0.1, 0.0, 0.0));
		linear_contact_vels.push_back(Vector3d::Zero());
		COPSolver solver;
		vector<vector<Vector3d>> boundary_points;
		vector<Vector3d> line_bdry_pts;
		line_bdry_pts.push_back(Vector3d::Zero());
		line_bdry_pts.push_back(Vector3d(0.2, 0, 0));
		boundary_points.push_back(line_bdry_pts);
		vector<Vector3d> pt_bdry_pts;
		boundary_points.push_back(pt_bdry_pts);
		vector<uint> patch_indices;
		patch_indices.push_back(0); // first line contact
		patch_indices.push_back(5); // then point contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::LINE);
		contact_types.push_back(ContactType::POINT);
		ContactCOPSolution sol = solver.solveStartWithPatchCentroidOnePointOneLine(
			0.1,
			A_constraint,
			rhs_constraint,
			boundary_points,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 4 ---" << endl;
		cout << "No omegas" << endl;
		cout << "Point contact sliding: lin vel, No line contact lin vel: rolling, " << endl;
		MatrixXd A_constraint = MatrixXd::Identity(8, 8);
		MatrixXd rhs_constraint = VectorXd::Zero(8);
		rhs_constraint(2) = -0.1;
		rhs_constraint(7) = -0.1;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d::Zero());
		omegaBs.push_back(Vector3d::Zero());
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d::Zero());
		linear_contact_vels.push_back(Vector3d(0.1, 0.1, 0.0));
		COPSolver solver;
		vector<vector<Vector3d>> boundary_points;
		vector<Vector3d> line_bdry_pts;
		line_bdry_pts.push_back(Vector3d::Zero());
		line_bdry_pts.push_back(Vector3d(0.2, 0, 0));
		boundary_points.push_back(line_bdry_pts);
		vector<Vector3d> pt_bdry_pts;
		boundary_points.push_back(pt_bdry_pts);
		vector<uint> patch_indices;
		patch_indices.push_back(0); // first line contact
		patch_indices.push_back(5); // then point contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::LINE);
		contact_types.push_back(ContactType::POINT);
		ContactCOPSolution sol = solver.solveStartWithPatchCentroidOnePointOneLine(
			0.1,
			A_constraint,
			rhs_constraint,
			boundary_points,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 5 ---" << endl;
		cout << "No omegas" << endl;
		cout << "No lin vels, rolling/impending. Point contact separating" << endl;
		MatrixXd A_constraint = MatrixXd::Identity(8, 8);
		MatrixXd rhs_constraint = VectorXd::Zero(8);
		rhs_constraint(2) = 0.1;
		rhs_constraint(5) = -0.1;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d::Zero());
		omegaBs.push_back(Vector3d::Zero());
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d::Zero());
		linear_contact_vels.push_back(Vector3d::Zero());
		COPSolver solver;
		vector<vector<Vector3d>> boundary_points;
		vector<Vector3d> pt_bdry_pts;
		boundary_points.push_back(pt_bdry_pts);
		vector<Vector3d> line_bdry_pts;
		line_bdry_pts.push_back(Vector3d::Zero());
		line_bdry_pts.push_back(Vector3d(0.2, 0, 0));
		boundary_points.push_back(line_bdry_pts);
		vector<uint> patch_indices;
		patch_indices.push_back(0); // first point contact
		patch_indices.push_back(3); // then line contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::POINT);
		contact_types.push_back(ContactType::LINE);
		ContactCOPSolution sol = solver.solveStartWithPatchCentroidOnePointOneLine(
			0.1,
			A_constraint,
			rhs_constraint,
			boundary_points,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 6 ---" << endl;
		cout << "No omegas" << endl;
		cout << "No lin vels, rolling/impending. Line contact separating" << endl;
		MatrixXd A_constraint = MatrixXd::Identity(8, 8);
		MatrixXd rhs_constraint = VectorXd::Zero(8);
		rhs_constraint(2) = -0.1;
		rhs_constraint(5) = 0.1;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d::Zero());
		omegaBs.push_back(Vector3d::Zero());
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d::Zero());
		linear_contact_vels.push_back(Vector3d::Zero());
		COPSolver solver;
		vector<vector<Vector3d>> boundary_points;
		vector<Vector3d> pt_bdry_pts;
		boundary_points.push_back(pt_bdry_pts);
		vector<Vector3d> line_bdry_pts;
		line_bdry_pts.push_back(Vector3d::Zero());
		line_bdry_pts.push_back(Vector3d(0.2, 0, 0));
		boundary_points.push_back(line_bdry_pts);
		vector<uint> patch_indices;
		patch_indices.push_back(0); // first point contact
		patch_indices.push_back(3); // then line contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::POINT);
		contact_types.push_back(ContactType::LINE);
		ContactCOPSolution sol = solver.solveStartWithPatchCentroidOnePointOneLine(
			0.1,
			A_constraint,
			rhs_constraint,
			boundary_points,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 7 ---" << endl;
		cout << "No omegas" << endl;
		cout << "Line contact impending: lin vel, No point contact lin vel: rolling, " << endl;
		MatrixXd A_constraint = MatrixXd::Identity(8, 8);
		MatrixXd rhs_constraint = VectorXd::Zero(8);
		rhs_constraint(0) = 0.3;
		rhs_constraint(2) = -0.1;
		rhs_constraint(7) = -0.1;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d::Zero());
		omegaBs.push_back(Vector3d::Zero());
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d::Zero());
		linear_contact_vels.push_back(Vector3d::Zero());
		COPSolver solver;
		vector<vector<Vector3d>> boundary_points;
		vector<Vector3d> line_bdry_pts;
		line_bdry_pts.push_back(Vector3d::Zero());
		line_bdry_pts.push_back(Vector3d(0.2, 0, 0));
		boundary_points.push_back(line_bdry_pts);
		vector<Vector3d> pt_bdry_pts;
		boundary_points.push_back(pt_bdry_pts);
		vector<uint> patch_indices;
		patch_indices.push_back(0); // first line contact
		patch_indices.push_back(5); // then point contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::LINE);
		contact_types.push_back(ContactType::POINT);
		ContactCOPSolution sol = solver.solveStartWithPatchCentroidOnePointOneLine(
			0.1,
			A_constraint,
			rhs_constraint,
			boundary_points,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 8 - Patch ---" << endl;
		cout << "No omegas" << endl;
		cout << "No lin vels, rolling/impending" << endl;
		MatrixXd A_constraint = MatrixXd::Identity(6, 6);
		MatrixXd rhs_constraint = VectorXd::Zero(6);
		rhs_constraint(2) = -0.1;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d::Zero());
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d::Zero());
		COPSolver solver;
		ContactPatch patch;
		patch.max_extent = 0.2;
		Circle* c = new Circle();
		c->center.setZero();
		c->radius = 0.1;
		patch._intersection_curves.push_back(c);
		patch._interior_point.setZero();

		vector<uint> patch_indices;
		patch_indices.push_back(0); // first line contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::SURFACE);
		ContactCOPSolution sol = solver.solveSurfaceContact(
			0.1,
			A_constraint,
			rhs_constraint,
			patch,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 9 - Patch ---" << endl;
		cout << "Separating" << endl;
		MatrixXd A_constraint = MatrixXd::Identity(6, 6);
		MatrixXd rhs_constraint = VectorXd::Zero(6);
		rhs_constraint(2) = 0.5;
		rhs_constraint(3) = 1.0;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d::Zero());
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d::Zero());
		COPSolver solver;
		ContactPatch patch;
		patch.max_extent = 0.2;
		Circle* c = new Circle();
		c->center.setZero();
		c->radius = 0.1;
		patch._intersection_curves.push_back(c);
		patch._interior_point.setZero();

		vector<uint> patch_indices;
		patch_indices.push_back(0); // first line contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::SURFACE);
		ContactCOPSolution sol = solver.solveSurfaceContact(
			0.1,
			A_constraint,
			rhs_constraint,
			patch,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 10 - Patch ---" << endl;
		cout << "Tilting X moment" << endl;
		MatrixXd A_constraint = MatrixXd::Identity(6, 6);
		MatrixXd rhs_constraint = VectorXd::Zero(6);
		rhs_constraint(2) = -0.1;
		rhs_constraint(3) = 1.0;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d::Zero());
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d::Zero());
		COPSolver solver;
		ContactPatch patch;
		patch.max_extent = 0.2;
		Circle* c = new Circle();
		c->center.setZero();
		c->radius = 0.1;
		patch._intersection_curves.push_back(c);
		patch._interior_point.setZero();

		vector<uint> patch_indices;
		patch_indices.push_back(0); // first line contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::SURFACE);
		ContactCOPSolution sol = solver.solveSurfaceContact(
			0.1,
			A_constraint,
			rhs_constraint,
			patch,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 11 - Patch ---" << endl;
		cout << "Tilting Y moment" << endl;
		MatrixXd A_constraint = MatrixXd::Identity(6, 6);
		MatrixXd rhs_constraint = VectorXd::Zero(6);
		rhs_constraint(2) = -0.1;
		rhs_constraint(4) = 1.0;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d::Zero());
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d::Zero());
		COPSolver solver;
		ContactPatch patch;
		patch.max_extent = 0.2;
		Circle* c = new Circle();
		c->center.setZero();
		c->radius = 0.1;
		patch._intersection_curves.push_back(c);
		patch._interior_point.setZero();

		vector<uint> patch_indices;
		patch_indices.push_back(0); // first line contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::SURFACE);
		ContactCOPSolution sol = solver.solveSurfaceContact(
			0.1,
			A_constraint,
			rhs_constraint,
			patch,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 12 - Patch ---" << endl;
		cout << "Linear slip" << endl;
		MatrixXd A_constraint = MatrixXd::Identity(6, 6);
		MatrixXd rhs_constraint = VectorXd::Zero(6);
		rhs_constraint(2) = -0.1;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d::Zero());
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d(0.1, 0, 0));
		COPSolver solver;
		ContactPatch patch;
		patch.max_extent = 0.2;
		Circle* c = new Circle();
		c->center.setZero();
		c->radius = 0.1;
		patch._intersection_curves.push_back(c);
		patch._interior_point.setZero();

		vector<uint> patch_indices;
		patch_indices.push_back(0); // first line contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::SURFACE);
		ContactCOPSolution sol = solver.solveSurfaceContact(
			0.1,
			A_constraint,
			rhs_constraint,
			patch,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 13 - Patch ---" << endl;
		cout << "Rotational slip" << endl;
		MatrixXd A_constraint = MatrixXd::Identity(6, 6);
		MatrixXd rhs_constraint = VectorXd::Zero(6);
		rhs_constraint(2) = -0.1;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d(0, 0, 0.1));
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d(0.0, 0, 0));
		COPSolver solver;
		ContactPatch patch;
		patch.max_extent = 0.2;
		Circle* c = new Circle();
		c->center.setZero();
		c->radius = 0.1;
		patch._intersection_curves.push_back(c);
		patch._interior_point.setZero();

		vector<uint> patch_indices;
		patch_indices.push_back(0); // first line contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::SURFACE);
		ContactCOPSolution sol = solver.solveSurfaceContact(
			0.1,
			A_constraint,
			rhs_constraint,
			patch,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 14 - Patch ---" << endl;
		cout << "Impending translational slip" << endl;
		MatrixXd A_constraint = MatrixXd::Identity(6, 6);
		MatrixXd rhs_constraint = VectorXd::Zero(6);
		rhs_constraint(2) = -0.1;
		rhs_constraint(1) = 0.1;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d(0, 0, 0.0));
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d(0.0, 0, 0));
		COPSolver solver;
		ContactPatch patch;
		patch.max_extent = 0.2;
		Circle* c = new Circle();
		c->center.setZero();
		c->radius = 0.1;
		patch._intersection_curves.push_back(c);
		patch._interior_point.setZero();

		vector<uint> patch_indices;
		patch_indices.push_back(0); // first line contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::SURFACE);
		ContactCOPSolution sol = solver.solveSurfaceContact(
			0.1,
			A_constraint,
			rhs_constraint,
			patch,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
	{
		cout << "--- TEST 15 - Patch ---" << endl;
		cout << "Impending rotational slip" << endl;
		MatrixXd A_constraint = MatrixXd::Identity(6, 6);
		MatrixXd rhs_constraint = VectorXd::Zero(6);
		rhs_constraint(2) = -0.1;
		rhs_constraint(5) = 0.1;
		vector<Vector3d> omegaAs;
		omegaAs.push_back(Vector3d::Zero());
		vector<Vector3d> omegaBs;
		omegaBs.push_back(Vector3d(0, 0, 0.0));
		vector<Vector3d> linear_contact_vels;
		linear_contact_vels.push_back(Vector3d(0.0, 0, 0));
		COPSolver solver;
		ContactPatch patch;
		patch.max_extent = 0.2;
		Circle* c = new Circle();
		c->center.setZero();
		c->radius = 0.1;
		patch._intersection_curves.push_back(c);
		patch._interior_point.setZero();

		vector<uint> patch_indices;
		patch_indices.push_back(0); // first line contact
		vector<ContactType> contact_types;
		contact_types.push_back(ContactType::SURFACE);
		ContactCOPSolution sol = solver.solveSurfaceContact(
			0.1,
			A_constraint,
			rhs_constraint,
			patch,
			patch_indices,
			contact_types,
			omegaAs,
			omegaBs,
			linear_contact_vels
		);
		printCOPSolution(sol);
	}
}
