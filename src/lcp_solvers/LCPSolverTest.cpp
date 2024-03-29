// LCP solver test file

#include <iostream>
#include "LCPSolver.h"
#include "LCPSolverInternal.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace Sai2LCPSolver;

int main (int argc, char** argv) {
	// tests are based on a loaded stick hitting the ground at two extremes.
	const double load = 5; // kg
	const double side = 1; // m, distance from load point to end 1 and end 2, which are also contact 1 and contact 2
	const double stick_end_mass = 0.2; // kg, on each end
	// In tests 1 - 8, stick is assumed to be aligned with X axis.
	MatrixXd M_inv(6, 6), Lambda_inv(6, 6), J1(3, 6), J2(3, 6);
	M_inv.setZero();
	M_inv.block(0, 0, 3, 3) = Matrix3d::Identity()* 1.0/(load + 2*stick_end_mass);
	M_inv(3, 3) = 1.0/0.00001; // rotation about x axis which is along length of rod
	M_inv(4, 4) = 1.0/(2*stick_end_mass*side*side);
	M_inv(5, 5) = 1.0/(2*stick_end_mass*side*side);
	J1 << 1, 0, 0, 0, 0, 0,
		0, 1, 0, 0, 0, side,
		0, 0, 1, 0, -side, 0;
	J2 << 1, 0, 0, 0, 0, 0,
		0, 1, 0, 0, 0, -side,
		0, 0, 1, 0, side, 0;
	Lambda_inv.block(0,0,3,3) = J1 * M_inv * J1.transpose();
	Lambda_inv.block(0,3,3,3) = J1 * M_inv * J2.transpose();
	Lambda_inv.block(3,3,3,3) = J2 * M_inv * J2.transpose();
	Lambda_inv.block(3,0,3,3) = J2 * M_inv * J1.transpose();

	// cout << "Lambda_inv" << endl;

	// cout << Lambda_inv << endl;

	// test 1 - two contacts, rolling - rolling, no tangential pre-velocity
	{
		cout << " --- Test 1 ---" << endl;
		VectorXd pre_v(6);
		const double epsilon = 1.0; // elastic
		const double mu = 0.1;
		pre_v << 0, 0, -0.1, 0, 0, -0.1;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		cout << "Vsol: " << (Lambda_inv*test1_sol.p_sol + pre_v).transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol new: " << (Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 2 - two contacts, rolling - rolling, no tangential pre-velocity
	{
		cout << endl;
		cout << " --- Test 2 ---" << endl;
		VectorXd pre_v(6);
		const double epsilon = 1.0; // elastic
		const double mu = 0.1;
		pre_v << 0, 0, -0.1, 0, 0, 0;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		cout << "Vsol: " << (Lambda_inv*test1_sol.p_sol + pre_v).transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol new: " << (Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 3 - two contacts, rolling - rolling, x-tangential pre-velocity
	{
		cout << endl;
		cout << " --- Test 3 ---" << endl;
		VectorXd pre_v(6);
		const double epsilon = 1.0; // elastic
		const double mu = 0.1;
		pre_v << 0.02, 0, -0.1, 0.02, 0, -0.1;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		cout << "Vsol: " << (Lambda_inv*test1_sol.p_sol + pre_v).transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol new: " << (Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 3.2 - two contacts, rolling - rolling, x-tangential pre-velocity
	{
		cout << endl;
		cout << " --- Test 3.2 ---" << endl;
		VectorXd pre_v(6);
		const double epsilon = 1.0; // elastic
		const double mu = 0.1;
		pre_v << -0.02, 0, -0.1, -0.02, 0, -0.1;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		cout << "Vsol: " << (Lambda_inv*test1_sol.p_sol + pre_v).transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol new: " << (Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 3.3 - two contacts, sliding - sliding, x-tangential pre-velocity
	{
		cout << endl;
		cout << " --- Test 3.3 ---" << endl;
		VectorXd pre_v(6);
		const double epsilon = 1.0; // elastic
		const double mu = 0.0;
		pre_v << 0.02, 0, -0.1, 0.02, 0, -0.1;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		cout << "Vsol: " << (Lambda_inv*test1_sol.p_sol + pre_v).transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol new: " << (Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 4 - two contacts, sliding - sliding, x-tangential pre-velocity
	{
		cout << endl;
		cout << " --- Test 4 ---" << endl;
		VectorXd pre_v(6);
		const double epsilon = 1.0; // elastic
		const double mu = 0.1;
		pre_v << 0.1, 0, -0.1, 0.1, 0, -0.1;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		cout << "Vsol: " << (Lambda_inv*test1_sol.p_sol + pre_v).transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol new: " << (Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 5 - two contacts, contact 1 sliding - contact 2 rolling, y-tangential pre-velocity
	{
		cout << endl;
		cout << " --- Test 5 ---" << endl;
		VectorXd pre_v(6);
		const double epsilon = 1.0; // elastic
		const double mu = 0.1;
		pre_v << 0.0, 0.1, -0.1, 0.0, 0.01, -0.1;
		// cout << "pre energy: " << stick_end_mass*(pre_v.head(3).norm() + pre_v.tail(3).norm()) + load*((pre_v.head(3)+pre_v.tail(3)).norm()/2) << endl;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		cout << "Vsol: " << (Lambda_inv*test1_sol.p_sol + pre_v).transpose() << endl;
		// VectorXd post_v = (Lambda_inv*test1_sol.p_sol + pre_v).transpose();
		// cout << "post energy: " << stick_end_mass*(post_v.head(3).norm() + post_v.tail(3).norm()) + load*((post_v.head(3)+post_v.tail(3)).norm()/2) << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol new: " << (Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 5.2 - two contacts, contact 1 rolling - contact 2 sliding, y-tangential pre-velocity
	{
		cout << endl;
		cout << " --- Test 5.2 ---" << endl;
		VectorXd pre_v(6);
		const double epsilon = 1.0; // elastic
		const double mu = 0.1;
		pre_v << 0.0, 0.01, -0.1, 0.0, 0.1, -0.1;
		// cout << "pre energy: " << stick_end_mass*(pre_v.head(3).norm() + pre_v.tail(3).norm()) + load*((pre_v.head(3)+pre_v.tail(3)).norm()/2) << endl;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		cout << "Vsol: " << (Lambda_inv*test1_sol.p_sol + pre_v).transpose() << endl;
		// VectorXd post_v = (Lambda_inv*test1_sol.p_sol + pre_v).transpose();
		// cout << "post energy: " << stick_end_mass*(post_v.head(3).norm() + post_v.tail(3).norm()) + load*((post_v.head(3)+post_v.tail(3)).norm()/2) << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol new: " << (Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}

	// test 6 - contact 1 only, rolling, no tangential pre-velocity
	{
		cout << endl;
		cout << " --- Test 6 ---" << endl;
		VectorXd pre_v(6);
		const double epsilon = 0.0; // elastic
		const double mu = 0.4;
		pre_v << 0.001, 0, -0.1, 0.001, 0, 0.2;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		cout << "Vsol: " << (Lambda_inv*test1_sol.p_sol + pre_v).transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol new: " << (Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 6.2 - contact 1 only, sliding, no tangential pre-velocity
	{
		cout << endl;
		cout << " --- Test 6.2 ---" << endl;
		VectorXd pre_v(6);
		const double epsilon = 0.0; // elastic
		const double mu = 0.1;
		pre_v << 0.001, 0, -0.1, 0.001, 0, 0.2;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		cout << "Vsol: " << (Lambda_inv*test1_sol.p_sol + pre_v).transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol new: " << (Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 7 - contact 2 only, rolling, no tangential pre-velocity
	{
		cout << endl;
		cout << " --- Test 7 ---" << endl;
		VectorXd pre_v(6);
		const double epsilon = 0.0; // elastic
		const double mu = 0.1;
		pre_v << 0, 0, 0.2, 0, 0, -0.1;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		cout << "Vsol: " << (Lambda_inv*test1_sol.p_sol + pre_v).transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol new: " << (Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 7.2 - contact 2 only, sliding, no tangential pre-velocity
	{
		cout << endl;
		cout << " --- Test 7.2 ---" << endl;
		VectorXd pre_v(6);
		const double epsilon = 0.0; // elastic
		const double mu = 0.1;
		pre_v << 0.001, 0, 0.2, 0.001, 0, -0.1;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		cout << "Vsol: " << (Lambda_inv*test1_sol.p_sol + pre_v).transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Lambda_inv, pre_v, pre_v, epsilon, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol new: " << (Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 8 - one contact only
	{
		cout << endl;
		cout << " --- Test 8 ---" << endl;
		VectorXd pre_v(6);
		const double mu = 0.01;
		// pre_v << -1.31951e-08, 0, -0.000262909, -2.44319e-11, 0, -0.000131114;
		pre_v << -6.6098e-09, 0, -0.000262909, -6.6098e-09, 0, -0.000131114;
		Eigen::MatrixXd Test8Lambda_inv(6,6);
		Test8Lambda_inv<< 1.33099, 0, 0.332461, 1.33092, 0, -0.334272,
        				0, 2.44252, 0, 0, 1.77578, 0,
 					0.332461, 0, 1.33175, 0.332394, 0, 0.665086,
  					1.33092, 0, 0.332394, 1.33086, 0, -0.334206,
        			0, 1.77578, 0, 0, 2.44238, 0,
					-0.334272, 0, 0.665086, -0.334206, 0, 1.33175;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Test8Lambda_inv, pre_v, pre_v, 0.0, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Test8Lambda_inv, pre_v, pre_v, 0.0, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 9
	{
		cout << endl;
		cout << " --- Test 9 ---" << endl;
		VectorXd pre_v(6);
		const double mu = 0.1;
		pre_v << 0.00118824, 0.00191278, -0.000855545, 0.00118824, -0.00117223, -0.000874077;
		Eigen::MatrixXd Test9Lambda_inv(6,6);
		Test9Lambda_inv<< 1.12296, -5.48628e-06, -0.123986, 1.12303, -5.48628e-06, 0.125942,
					-5.48628e-06, 2.23404, -0.000184815, -5.48628e-06, 1.98412, -0.000184815,
   					-0.123986, -0.000184815, 1.12298, -0.124058, -0.000184815, 0.872976,
     				1.12303, -5.48628e-06, -0.124058, 1.12311, -5.48628e-06, 0.126014,
				-5.48628e-06, 1.98412, -0.000184815, -5.48628e-06, 2.23419, -0.000184815,
    				0.125942, -0.000184815, 0.872976, 0.126014, -0.000184815, 1.12298;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Test9Lambda_inv, pre_v, pre_v, 0.0, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Test9Lambda_inv, pre_v, pre_v, 0.0, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 10
	{
		cout << endl;
		cout << " --- Test 10 ---" << endl;
		VectorXd pre_v(6);
		const double mu = 0.1;
		pre_v << -7.01934e-08, -0.00129741, -4.55199e-08, -7.02126e-08, -0.00386053, -1.12357e-07;
		Eigen::MatrixXd Test9Lambda_inv(6,6);
		Test9Lambda_inv<< 1.12296,  -5.4884e-06,    -0.123986,      1.12303,  -5.4884e-06,     0.125942,
 -5.4884e-06,      2.23404, -0.000184889,  -5.4884e-06,      1.98412, -0.000184889,
   -0.123986, -0.000184889,      1.12298,    -0.124058, -0.000184889,     0.872976,
     1.12303,  -5.4884e-06,    -0.124058,      1.12311,  -5.4884e-06,     0.126014,
 -5.4884e-06,      1.98412, -0.000184889,  -5.4884e-06,      2.23419, -0.000184889,
    0.125942, -0.000184889,     0.872976,     0.126014, -0.000184889,      1.12298;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Test9Lambda_inv, pre_v, pre_v, 0.0, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Test9Lambda_inv, pre_v, pre_v, 0.0, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 11
	{
		cout << endl;
		cout << " --- Test 11 ---" << endl;
		VectorXd pre_v(6);
		const double mu = 0.1;
		pre_v << 0,        0, -1.33997,        0,        0, -1.33997;
		Eigen::MatrixXd Test9Lambda_inv(6,6);
		Test9Lambda_inv<< 1.12201,        0,   -0.125,  1.12201,        0,    0.125,
       0,  2.23412,        0,        0,  1.98412,        0,
  -0.125,        0,    1.124,   -0.125,        0, 0.874001,
 1.12201,        0,   -0.125,  1.12201,        0,    0.125,
       0,  1.98412,        0,        0,  2.23412,        0,
   0.125,        0, 0.874001,    0.125,        0,    1.124;
		CollLCPPointSolution test1_sol = solveCollLCPPoint (2, Test9Lambda_inv, pre_v, pre_v, 0.3, mu);
		cout << static_cast<int>(test1_sol.result) << endl;
		cout << "Psol: " << test1_sol.p_sol.transpose() << endl;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Test9Lambda_inv, pre_v, pre_v, 0.3, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;

		cout << "Psol diff: " << (test1_sol.p_sol - test1_sol_new.p_sol).transpose() << endl;
	}
	// test 12
	{
		cout << endl;
		cout << " --- Test 12 ---" << endl;
		VectorXd pre_v(9);
		const double mu = 0.1;
		pre_v << -1.26939e-05,    -0.055485, -3.05108e-09, -1.26939e-05 ,  -0.0554849 ,-3.05108e-09 ,  5.6805e-05 ,     0.15985,  8.67281e-05;
		Eigen::MatrixXd Test9Lambda_inv(9,9);
		Test9Lambda_inv<< 1.12201,   1.0721e-09,       -0.125  ,    1.12201 ,  1.0721e-09 ,       0.125,    -0.888859, -3.82536e-05, -6.97898e-05,
  1.0721e-09,      2.23412,  1.38778e-17,   1.0721e-09,      1.98412,  1.99493e-17,   -0.0626027,     0.247576,     -0.50048,
      -0.125,  1.38778e-17,        1.124,       -0.125,  1.56125e-17,     0.874001,    -0.108018,      0.50097,    -0.864415,
     1.12201,   1.0721e-09,       -0.125,      1.12201,   1.0721e-09,        0.125,    -0.888859, -3.82536e-05, -6.97898e-05,
  1.0721e-09,      1.98412,  1.47451e-17,   1.0721e-09,      2.23412,  1.73472e-17,    0.0627556,     0.247708,    -0.500384,
       0.125,  2.08167e-17,     0.874001,        0.125,  1.73472e-17,        1.124,     0.108282,     0.500894,     -0.86425,
   -0.888859,   -0.0626027,    -0.108018,   -0.888859 ,   0.0627556 ,    0.108282 ,     3.43311, -0.000834629, -1.73007e-05,
-3.82536e-05,     0.247576,      0.50097, -3.82536e-05,     0.247708,     0.500894, -0.000834629,      3.29645,    -0.152076,
-6.97898e-05,     -0.50048,    -0.864415, -6.97898e-05,    -0.500384,     -0.86425, -1.73007e-05,    -0.152076 ,     2.11717;
		auto solver = Sai2LCPSolver::LCPSolver();
		CollLCPPointSolution test1_sol_new = solver.solve(Test9Lambda_inv, pre_v, pre_v, 0.0, mu);
		cout << static_cast<int>(test1_sol_new.result) << endl;
		cout << "Psol new: " << test1_sol_new.p_sol.transpose() << endl;
		cout << "Vsol: " << (Test9Lambda_inv*test1_sol_new.p_sol + pre_v).transpose() << endl;
	}
	return 0;
}
