// LCP solver test file

#include <iostream>
#include "LCPSolver.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

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
	}
	// test 8 - one contact only
	{
		cout << endl;
		cout << " --- Test 8 ---" << endl;
		VectorXd pre_v(6);
		const double mu = 0.01;
		pre_v << -1.31951e-08, 0, -0.000262909, -2.44319e-11, 0, -0.000131114;
		//-6.6098e-09, 0, -0.000262909, -6.6098e-09, 0, -0.000131114;
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
	}
	return 0;
}
