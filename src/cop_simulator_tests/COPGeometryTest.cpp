// COPGeometryTest.cpp

#include <iostream>
#include <Eigen/Dense>
#include "cop_simulator/geometry/Primitive.h"

using namespace std;
using namespace Eigen;
using namespace Sai2COPSim;

void printPrimPrimInfo(const PrimPrimContactInfo& info) {
	cout << "Min dist: " << info.min_distance << "\n";
	if(info.type == ContactType::POINT) cout << "Type: Point" << "\n";
	if(info.type == ContactType::LINE) cout << "Type: Line" << "\n";
	if(info.type == ContactType::SURFACE) cout << "Type: Surface" << "\n";
	if(info.type == ContactType::UNDEFINED) cout << "Type: Undefined" << "\n";
	cout << "Normal: " << info.normal_dir.transpose() << "\n";
	cout << "Constraint dir 1: " << info.constraint_dir1.transpose() << "\n";
	cout << "Constraint dir 2: " << info.constraint_dir2.transpose() << "\n";
	for(const auto& pt: info.contact_points) {
		cout << "Point: " << pt.transpose() << "\n";
	}
	cout << endl;
}

void testPlaneCapsuleDistance();
void testCapsuleCapsuleDistance();
void testPlaneCylinderDistance();

int main (int argc, char** argv) {
	// testPlaneCapsuleDistance();
	// testCapsuleCapsuleDistance();
	testPlaneCylinderDistance();
	return 0;
}

void testPlaneCapsuleDistance() {
	cout << "000000000000 ----- PLANE - CAPSULE ----- 000000000000" << endl;
	{
		// test distance between z-plane and capsule
		cout << "--- TEST 1 ---" << endl; 
		PlanePrimitive* plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		CapsulePrimitive* capsule = new CapsulePrimitive("tc", 0.1, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d capTF = Affine3d::Identity();
		capTF.translation() << 0.0, 0.0, 0.5;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(plane, planeTF, capsule, capTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between z-plane and capsule, with z-plane displacement
		cout << "--- TEST 2 ---" << endl; 
		PlanePrimitive* plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		CapsulePrimitive* capsule = new CapsulePrimitive("tc", 0.1, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		planeTF.translation() << 0.1, 0.4, -0.3;
		Affine3d capTF = Affine3d::Identity();
		capTF.translation() << 0.0, 0.0, 0.5;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(plane, planeTF, capsule, capTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between z-plane and capsule, with plane rotation about X
		cout << "--- TEST 3 ---" << endl; 
		PlanePrimitive* plane = new PlanePrimitive("tp", Vector3d(0.0, cos(M_PI/4), sin(M_PI/4)), Vector3d(0.0, 0.0, 0.0));
		CapsulePrimitive* capsule = new CapsulePrimitive("tc", 0.1, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d capTF = Affine3d::Identity();
		capTF.translation() << 0.0, 0.0, 0.5;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(plane, planeTF, capsule, capTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between z-plane and capsule, with plane rotation about Y
		cout << "--- TEST 4 ---" << endl; 
		PlanePrimitive* plane = new PlanePrimitive("tp", Vector3d(cos(M_PI/4), 0.0, sin(M_PI/4)), Vector3d(0.0, 0.0, 0.0));
		CapsulePrimitive* capsule = new CapsulePrimitive("tc", 0.1, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d capTF = Affine3d::Identity();
		capTF.translation() << 0.0, 0.0, 0.5;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(plane, planeTF, capsule, capTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between z-plane and capsule, with capsule rotation
		cout << "--- TEST 5 ---" << endl; 
		PlanePrimitive* plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		CapsulePrimitive* capsule = new CapsulePrimitive("tc", 0.1, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d capTF = Affine3d::Identity();
		capTF.translation() << 0.5, 0.5, 0.5;
		capTF.linear() << 0.0, 1.0, 0.0,
							-1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(plane, planeTF, capsule, capTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between z-plane and capsule, with capsule Y rotation
		cout << "--- TEST 6 ---" << endl; 
		PlanePrimitive* plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		CapsulePrimitive* capsule = new CapsulePrimitive("tc", 0.1, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d capTF = Affine3d::Identity();
		capTF.translation() << 0.5, 0.5, 0.5;
		capTF.linear() << 0.0, 0.0, -1.0,
							0.0, 1.0, 0.0,
							1.0, 0.0, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(plane, planeTF, capsule, capTF);
		printPrimPrimInfo(info);
	}
}

void testCapsuleCapsuleDistance() {
	cout << "000000000000 ----- CAPSULE - CAPSULE ----- 000000000000" << endl;
	{
		// test distance between parallel cylinders, aligned exactly
		cout << "--- TEST 1 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.5, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, aligned exactly
		cout << "--- TEST 2 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.2, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, -0.5, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, aligned exactly
		cout << "--- TEST 3 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.0, 1.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, aligned exactly
		cout << "--- TEST 4 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.4);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.0, -1.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, aligned exactly
		cout << "--- TEST 5 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.4);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.0, -1.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, aligned exactly, intersecting
		cout << "--- TEST 6 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.15, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, aligned exactly, displaced
		cout << "--- TEST 7 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.15, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.15, 0.2);
		Affine3d capATF = Affine3d::Identity();
		capATF.translation() << 0.0, 0.1, 0.0;
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.4, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, aligned exactly, rotated
		cout << "--- TEST 8 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		capATF.linear() << 0.0, 1.0, 0.0,
							-1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		Affine3d capBTF = Affine3d::Identity();
		capBTF.linear() << 0.0, 1.0, 0.0,
							-1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		capBTF.translation() << 0.3, 0.0, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, aligned inexactly
		cout << "--- TEST 9 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.1, 0.3, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, aligned inexactly far away
		cout << "--- TEST 10 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.5, 0.3, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, aligned end to end
		cout << "--- TEST 11 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.5, 0.0, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, aligned end to end, intersecting
		cout << "--- TEST 12 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.3, 0.0, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, capB rotated 90deg, in plane
		cout << "--- TEST 13 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.6, 0.0;
		capBTF.linear() << 0.0, -1.0, 0.0,
							1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, capA rotated 90deg, in plane
		cout << "--- TEST 14 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		capATF.linear() << 0.0, -1.0, 0.0,
							1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.6, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, capB rotated 45deg, in plane
		cout << "--- TEST 15 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.linear() << cos(M_PI/4.0), -sin(M_PI/4.0), 0.0,
							sin(M_PI/4.0), cos(M_PI/4.0), 0.0,
							0.0, 0.0, 1.0;
		capBTF.translation() << 0.0, 0.6, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, capA rotated 45deg, in plane
		cout << "--- TEST 16 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		capATF.linear() << cos(M_PI/4.0), -sin(M_PI/4.0), 0.0,
							sin(M_PI/4.0), cos(M_PI/4.0), 0.0,
							0.0, 0.0, 1.0;
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.6, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, capB rotated 90deg, in plane, end to end
		cout << "--- TEST 17 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.6, 0.6, 0.0;
		capBTF.linear() << 0.0, -1.0, 0.0,
							1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, capB rotated 90deg, out of plane
		cout << "--- TEST 18 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.6, 0.3;
		capBTF.linear() << 0.0, -1.0, 0.0,
							1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, capB rotated 90deg, out of plane, intersecting
		cout << "--- TEST 19 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.1, 0.15;
		capBTF.linear() << 0.0, -1.0, 0.0,
							1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel cylinders, aligned end to end, ambiguous case
		cout << "--- TEST LAST. THIS ASSERTS IN DEBUG ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		auto info = PrimPrimDistance::distancePrimitivePrimitive(capsuleA, capATF, capsuleB, capBTF);	
		printPrimPrimInfo(info);
	}
}

void testPlaneCylinderDistance() {
	cout << "000000000000 ----- PLANE - CYLINDER ----- 000000000000" << endl;
	{
		// test distance between z-plane and cylinder
		cout << "--- TEST 1 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto cylinder = new CylinderPrimitive("tc", 0.1, 0.2, 3);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d cylTF = Affine3d::Identity();
		cylTF.translation() << 0.0, 0.0, 0.5;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(plane, planeTF, cylinder, cylTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between z-plane and cylinder
		cout << "--- TEST 2 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto cylinder = new CylinderPrimitive("tc", 0.1, 0.2, 3);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d cylTF = Affine3d::Identity();
		cylTF.translation() << 0.0, 0.0, 0.5;
		cylTF.linear() << 0.0, -1.0, 0.0,
							1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(plane, planeTF, cylinder, cylTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between z-plane and cylinder: line contact
		cout << "--- TEST 3 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto cylinder = new CylinderPrimitive("tc", 0.1, 0.2, 3);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d cylTF = Affine3d::Identity();
		cylTF.translation() << 0.0, 0.0, 0.5;
		cylTF.linear() << 1.0, 0.0, 0.0,
							0.0, 0.0, -1.0,
							0.0, 1.0, 0.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(plane, planeTF, cylinder, cylTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between z-plane and cylinder: inverted plane contact
		cout << "--- TEST 4 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto cylinder = new CylinderPrimitive("tc", 0.1, 0.2, 3);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d cylTF = Affine3d::Identity();
		cylTF.translation() << 0.0, 0.0, 0.5;
		cylTF.linear() << 1.0, 0.0, 0.0,
							0.0, -1.0, 0.0,
							0.0, 0.0, -1.0;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(plane, planeTF, cylinder, cylTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between z-plane and cylinder: point contact, no intersection
		cout << "--- TEST 5 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto cylinder = new CylinderPrimitive("tc", 0.1, 0.2, 3);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d cylTF = Affine3d::Identity();
		cylTF.translation() << 0.0, 0.0, 0.5;
		cylTF.linear() << 1.0, 0.0, 0.0,
							0.0, cos(M_PI/4), -sin(M_PI/4),
							0.0, sin(M_PI/4), cos(M_PI/4);
		auto info = PrimPrimDistance::distancePrimitivePrimitive(plane, planeTF, cylinder, cylTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between z-plane and cylinder
		cout << "--- TEST 6 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto cylinder = new CylinderPrimitive("tc", 0.1, 0.2, 6);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d cylTF = Affine3d::Identity();
		cylTF.translation() << 0.0, 0.0, -0.02;
		auto info = PrimPrimDistance::distancePrimitivePrimitive(plane, planeTF, cylinder, cylTF);
		printPrimPrimInfo(info);
	}
}