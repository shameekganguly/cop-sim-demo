// COPGeometryTest.cpp

#include <iostream>
#include <Eigen/Dense>
#include "cop_simulator/geometry/Primitive.h"
#include "cop_simulator/geometry/Composite1PkN.h"
#include "cop_simulator/geometry/PrimPrimDistance.h"
#include "cop_simulator/geometry/IntersectionEdge.h"

using namespace std;
using namespace Eigen;
using namespace Sai2COPSim;

void printPrimPrimInfo(const PrimPrimContactInfo& info) {
	cout << "Min dist: " << info.min_distance << "\n";
	if(info.type == ContactType::POINT) cout << "Type: Point" << "\n";
	if(info.type == ContactType::LINE) cout << "Type: Line" << "\n";
	if(info.type == ContactType::SURFACE) cout << "Type: Surface" << "\n";
	if(info.type == ContactType::CONCAVE) cout << "Type: Concave" << "\n";
	if(info.type == ContactType::UNDEFINED) cout << "Type: Undefined" << "\n";
	if(info.type != ContactType::CONCAVE) {
		cout << "Normal: " << info.normal_dir.transpose() << "\n";
		cout << "Constraint dir 1: " << info.constraint_dir1.transpose() << "\n";
		cout << "Constraint dir 2: " << info.constraint_dir2.transpose() << "\n";
		for(const auto& pt: info.contact_points) {
			cout << "Point: " << pt.transpose() << "\n";
		}
	} else {
		for(uint i = 0; i < info.contact_points.size(); i++) {
			cout << "Point " << i << ": " << info.contact_points[i].transpose()
				 << "\n Distance: " << info.distances[i]
				 << "\n Normal: " << info.normal_dirs[i].transpose()
				 << "\n Constraint dir 1: " << info.constraint_dir1s[i].transpose()
				 << "\n Constraint dir 2: " << info.constraint_dir2s[i].transpose() << "\n";
		}
	}

	cout << endl;
}

void printContactPatch(const ContactPatch& patch) {
	cout << "Contact patch:" << endl;
	cout << "Interior point: " << patch._interior_point.transpose() << endl;
	cout << "Num curves: " << patch._intersection_curves.size() << endl;
	cout << "Num line segments: " << patch._line_segments.size() << endl;
	for(const auto& lineseg: patch._line_segments) {
		cout << "Line seg: (";
		cout << lineseg.point1.transpose() << "), (";
		cout << lineseg.point2.transpose() << ")" << endl;
	}
	cout << "Max extent: " << patch.max_extent << endl;
	cout << endl;
}

void printPointTestResult(const PointTestResult& result) {
	cout << "Test point result:" << endl;
	cout << "Is in patch: " << result.f_is_in_patch << endl;
	cout << "Is on line seg: " << result.f_is_on_line_seg << endl;
	cout << "Line seg id: " << result.line_seg_id << endl;
	cout << "Is on vertex: " << result.f_is_on_vertex << endl;
	cout << "Distance to boundary: " << result.min_dist_to_boundary << endl;
	cout << endl;
}

void testPlaneCapsuleDistance();
void testCapsuleCapsuleDistance();
void testPlaneCylinderDistance();
void testNegCapsuleCapsuleDistance();
void testCircleCapsuleDistance();

int main (int argc, char** argv) {
	// testPlaneCapsuleDistance();
	// testCapsuleCapsuleDistance();
	// testPlaneCylinderDistance();
	// testNegCapsuleCapsuleDistance();
	testCircleCapsuleDistance();
	return 0;
}

void testPlaneCapsuleDistance() {
	PrimPrimContactInfo info;
	cout << "000000000000 ----- PLANE - CAPSULE ----- 000000000000" << endl;
	{
		// test distance between z-plane and capsule
		cout << "--- TEST 1 ---" << endl; 
		PlanePrimitive* plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		CapsulePrimitive* capsule = new CapsulePrimitive("tc", 0.1, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d capTF = Affine3d::Identity();
		capTF.translation() << 0.0, 0.0, 0.5;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, capsule, capTF);
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
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, capsule, capTF);
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
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, capsule, capTF);
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
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, capsule, capTF);
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
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, capsule, capTF);
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
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, capsule, capTF);
		printPrimPrimInfo(info);
	}
}

void testCapsuleCapsuleDistance() {
	PrimPrimContactInfo info;
	cout << "000000000000 ----- CAPSULE - CAPSULE ----- 000000000000" << endl;
	{
		// test distance between parallel capsules, aligned exactly
		cout << "--- TEST 1 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.5, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, aligned exactly
		cout << "--- TEST 2 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.2, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, -0.5, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, aligned exactly
		cout << "--- TEST 3 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.0, 1.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, aligned exactly
		cout << "--- TEST 4 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.4);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.0, -1.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, aligned exactly
		cout << "--- TEST 5 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.4);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.0, -1.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, aligned exactly, intersecting
		cout << "--- TEST 6 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.15, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, aligned exactly, displaced
		cout << "--- TEST 7 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.15, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.15, 0.2);
		Affine3d capATF = Affine3d::Identity();
		capATF.translation() << 0.0, 0.1, 0.0;
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.4, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, aligned exactly, rotated
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
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, aligned inexactly
		cout << "--- TEST 9 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.1, 0.3, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, aligned inexactly far away
		cout << "--- TEST 10 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.5, 0.3, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, aligned end to end
		cout << "--- TEST 11 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.5, 0.0, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, aligned end to end, intersecting
		cout << "--- TEST 12 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.3, 0.0, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, capB rotated 90deg, in plane
		cout << "--- TEST 13 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.6, 0.0;
		capBTF.linear() << 0.0, -1.0, 0.0,
							1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, capA rotated 90deg, in plane
		cout << "--- TEST 14 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		capATF.linear() << 0.0, -1.0, 0.0,
							1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.6, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, capB rotated 45deg, in plane
		cout << "--- TEST 15 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.linear() << cos(M_PI/4.0), -sin(M_PI/4.0), 0.0,
							sin(M_PI/4.0), cos(M_PI/4.0), 0.0,
							0.0, 0.0, 1.0;
		capBTF.translation() << 0.0, 0.6, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, capA rotated 45deg, in plane
		cout << "--- TEST 16 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		capATF.linear() << cos(M_PI/4.0), -sin(M_PI/4.0), 0.0,
							sin(M_PI/4.0), cos(M_PI/4.0), 0.0,
							0.0, 0.0, 1.0;
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.6, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, capB rotated 90deg, in plane, end to end
		cout << "--- TEST 17 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.6, 0.6, 0.0;
		capBTF.linear() << 0.0, -1.0, 0.0,
							1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, capB rotated 90deg, out of plane
		cout << "--- TEST 18 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.6, 0.3;
		capBTF.linear() << 0.0, -1.0, 0.0,
							1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, capB rotated 90deg, out of plane, intersecting
		cout << "--- TEST 19 ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.1, 0.15;
		capBTF.linear() << 0.0, -1.0, 0.0,
							1.0, 0.0, 0.0,
							0.0, 0.0, 1.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, capB rotated 45deg, in plane, long capA
		cout << "--- TEST 20 ---" << endl;
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 2.0);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.linear() << cos(M_PI/4.0), -sin(M_PI/4.0), 0.0,
							sin(M_PI/4.0), cos(M_PI/4.0), 0.0,
							0.0, 0.0, 1.0;
		capBTF.translation() << 0.0, 0.6, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, aligned end to end, ambiguous case
		cout << "--- TEST LAST. THIS ASSERTS IN DEBUG ---" << endl; 
		CapsulePrimitive* capsuleA = new CapsulePrimitive("tc1", 0.1, 0.2);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tc2", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);	
		printPrimPrimInfo(info);
	}
}

void testPlaneCylinderDistance() {
	PrimPrimContactInfo info;
	cout << "000000000000 ----- PLANE - CYLINDER ----- 000000000000" << endl;
	{
		// test distance between z-plane and cylinder
		cout << "--- TEST 1 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto cylinder = new CylinderPrimitive("tc", 0.1, 0.2, 3);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d cylTF = Affine3d::Identity();
		cylTF.translation() << 0.0, 0.0, 0.5;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, cylinder, cylTF);
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
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, cylinder, cylTF);
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
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, cylinder, cylTF);
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
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, cylinder, cylTF);
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
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, cylinder, cylTF);
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
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, cylinder, cylTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between z-plane and cylinder
		cout << "--- TEST 7 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto cylinder = new CylinderPrimitive("tc", 0.1, 0.2, 3);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d cylTF = Affine3d::Identity();
		cylTF.translation() << 0.0, 0.0, 0.05;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, cylinder, cylTF);
		printPrimPrimInfo(info);
		// check contact patch info
		cout << "Contact patch interior point: " << info.contact_patch._interior_point.transpose() << endl;
		cout << "Contact patch max extent: " << info.contact_patch.max_extent << endl;
		auto result = info.contact_patch.testPoint(Vector2d(0.05, 0.0));
		cout << "Contact patch point test: is in patch: " << result.f_is_in_patch << endl;
		cout << "Contact patch point test: min dist to boundary: " << result.min_dist_to_boundary << endl;
		auto result2 = info.contact_patch.testPoint(Vector2d(0.05, 0.2));
		cout << "Contact patch point test: is in patch: " << result2.f_is_in_patch << endl;
		cout << "Contact patch point test: min dist to boundary: " << result2.min_dist_to_boundary << endl;
		auto result3 = info.contact_patch.distanceFromBoundaryAlongRay(Vector2d(0.0, 0.05), Vector2d(0.0, -1.0));
		cout << "Contact patch point test: dist along ray: " << result3 << endl;
		auto result4 = info.contact_patch.distanceFromBoundaryAlongRay(Vector2d(0.0, 0.05), Vector2d(1.0, 0.0));
		cout << "Contact patch point test: dist along ray: " << result4 << endl;
	}
	{
		// test distance between z-plane and box
		cout << "--- TEST 8 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto box = new BoxPrimitive("tb", 0.2, 0.2, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d boxTF = Affine3d::Identity();
		boxTF.translation() << 0.0, 0.0, 0.2;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, box, boxTF);
		printPrimPrimInfo(info);
		printContactPatch(info.contact_patch);
		Vector2d test(0, 0);
		cout << "distance to patch boundary from test (" << test.transpose() << "): " << endl;
		printPointTestResult(info.contact_patch.testPoint(test));
		cout << "distance from line seg 1 to test (" << test.transpose() << "): "
							<< info.contact_patch._line_segments[0].distanceToPoint(test)
							<< endl;
		Vector2d test2(0, 0.1);
		cout << "distance to patch boundary from test (" << test2.transpose() << "): " << endl;
		printPointTestResult(info.contact_patch.testPoint(test2));

		Vector2d test3(-0.1, -0.1);
		cout << "distance to patch boundary from test (" << test3.transpose() << "): " << endl;
		printPointTestResult(info.contact_patch.testPoint(test3));

	}
	{
		// test distance between -ve-x-plane and box
		cout << "--- TEST 9 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(-1.0, 0.0, 0.0), Vector3d(0.0, 0.0, 0.0));
		auto box = new BoxPrimitive("tb", 0.1, 0.2, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d boxTF = Affine3d::Identity();
		boxTF.translation() << -5.0, 0.0, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, box, boxTF);
		printPrimPrimInfo(info);
		printContactPatch(info.contact_patch);
		Vector2d test(0, 0);
		cout << "distance to patch boundary from test (" << test.transpose() << "): " << endl;
		printPointTestResult(info.contact_patch.testPoint(test));
		cout << "distance from line seg 1 to test (" << test.transpose() << "): "
							<< info.contact_patch._line_segments[0].distanceToPoint(test)
							<< endl;
		Vector2d test2(0, 0.1);
		cout << "distance to patch boundary from test (" << test2.transpose() << "): " << endl;
		printPointTestResult(info.contact_patch.testPoint(test2));
		Vector2d test3(4, 4);
		cout << "distance to patch boundary from test (" << test3.transpose() << "): " << endl;
		printPointTestResult(info.contact_patch.testPoint(test3));
	}
	{
		// test distance between z-plane and box
		cout << "--- TEST 10 ---" << endl;
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto box = new BoxPrimitive("tb", 0.2, 0.2, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d boxTF = Affine3d::Identity();
		boxTF.translation() << 0.0, 0.0, 0.4;
		boxTF.linear() << 0.707, 0, 0.707,
							0,   1,   0,
						  -0.707, 0, 0.707;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, box, boxTF);
		printPrimPrimInfo(info);
	}
	{
		// test line contact distance between z-plane and box
		cout << "--- TEST 11 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto box = new BoxPrimitive("tb", 0.2, 0.2, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d boxTF = Affine3d::Identity();
		boxTF.translation() << 0.0, 0.0, 0.4;
		boxTF.linear() << 0.707, 0, 0.707,
							0,   1,   0,
						  -0.707, 0, 0.707;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, box, boxTF);
		printPrimPrimInfo(info);
	}
	{
		// test point contact distance between z-plane and box
		cout << "--- TEST 12 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto box = new BoxPrimitive("tb", 0.2, 0.2, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d boxTF = Affine3d::Identity();
		boxTF.translation() << 0.0, 0.0, 0.4;
		boxTF.linear() << 0.707, 	-0.707, 	0,
						  0.499849, 0.499849, -0.707,
						  0.499849, 0.499849,  0.707;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, box, boxTF);
		printPrimPrimInfo(info);
	}
	{
		// test surface contact between z-plane and 4-sided pyramid
		cout << "--- TEST 13 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto pyramid = new PyramidPrimitive("tpyr", 4, 0.1, 0.2);
		cout << "Incircle radius: " << pyramid->incircleRadius() << endl;
		cout << "Circumcircle radius: " << pyramid->circumRadius() << endl;
		cout << "Side edge angle: " << pyramid->sideEdgeAngle() << endl;
		cout << "Side face angle: " << pyramid->sideFaceAngle() << endl;
		cout << "Included angle: " << pyramid->includedAngle() << endl;
		Affine3d planeTF = Affine3d::Identity();
		Affine3d pyramidTF = Affine3d::Identity();
		pyramidTF.translation() << 0.0, 0.0, 0.2;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, pyramid, pyramidTF);
		printPrimPrimInfo(info);
		printContactPatch(info.contact_patch);
	}
	{
		// test point contact between z-plane and 4-sided pyramid apex
		cout << "--- TEST 14 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto pyramid = new PyramidPrimitive("tpyr", 4, 0.1, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d pyramidTF = Affine3d::Identity();
		pyramidTF.translation() << 0.0, 0.0, 0.4;
		pyramidTF.linear() << 1.0, 0.0, 0.0,
							0.0, -1.0, 0.0,
							0.0, 0.0, -1.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, pyramid, pyramidTF);
		printPrimPrimInfo(info);
	}
	{
		// test point contact between z-plane and 4-sided pyramid tilted
		cout << "--- TEST 15 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto pyramid = new PyramidPrimitive("tpyr", 4, 0.1, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d pyramidTF = Affine3d::Identity();
		pyramidTF.translation() << 0.0, 0.0, 0.4;
		pyramidTF.linear() << 1.0, 0.0, 0.0,
							0.0, 0.707, -0.707,
							0.0, 0.707, 0.707;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, pyramid, pyramidTF);
		printPrimPrimInfo(info);
	}
	{
		// test line contact between z-plane and 4-sided pyramid tilted, base
		cout << "--- TEST 16 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto pyramid = new PyramidPrimitive("tpyr", 4, 0.1, 0.2);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d pyramidTF = Affine3d::Identity();
		pyramidTF.translation() << 0.0, 0.0, 0.4;
		pyramidTF.linear() = (Eigen::AngleAxisd(M_PI/4, Vector3d::UnitX()) *
								Eigen::AngleAxisd(M_PI/4, Vector3d::UnitZ())).toRotationMatrix();
		cout << pyramidTF.linear() << endl;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, pyramid, pyramidTF);
		printPrimPrimInfo(info);
	}
	{
		// test face contact between z-plane and 4-sided pyramid tilted, base
		cout << "--- TEST 17 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto pyramid = new PyramidPrimitive("tpyr", 4, 0.1, 0.05);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d pyramidTF = Affine3d::Identity();
		pyramidTF.translation() << 0.0, 0.0, 0.4;
		pyramidTF.linear() = (Eigen::AngleAxisd(3*M_PI/4, Vector3d::UnitX()) *
								Eigen::AngleAxisd(M_PI/4, Vector3d::UnitZ())).toRotationMatrix();
		// cout << pyramidTF.linear() << endl;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, pyramid, pyramidTF);
		printPrimPrimInfo(info);
		printContactPatch(info.contact_patch);
	}
	{
		// test point contact between z-plane and 4-sided pyramid tilted, base
		cout << "--- TEST 18 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto pyramid = new PyramidPrimitive("tpyr", 4, 0.1, 0.05);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d pyramidTF = Affine3d::Identity();
		pyramidTF.translation() << 0.0, 0.0, 0.4;
		pyramidTF.linear() = (Eigen::AngleAxisd(-M_PI/8, Vector3d::UnitX()) *
								Eigen::AngleAxisd(M_PI/8, Vector3d::UnitZ())).toRotationMatrix();
		// cout << pyramidTF.linear() << endl;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, pyramid, pyramidTF);
		printPrimPrimInfo(info);
	}
	{
		// test edge contact between z-plane and 4-sided pyramid tilted, base
		cout << "--- TEST 19 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto pyramid = new PyramidPrimitive("tpyr", 4, 0.1, 0.1);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d pyramidTF = Affine3d::Identity();
		pyramidTF.translation() << 0.0, 0.0, 0.4;
		pyramidTF.linear() = (Eigen::AngleAxisd((M_PI/2+pyramid->sideEdgeAngle()), Vector3d::UnitX())).toRotationMatrix();
		// cout << pyramidTF.linear() << endl;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, pyramid, pyramidTF);
		printPrimPrimInfo(info);
	}
	{
		// test surface contact between z-plane and 6-sided pyramid
		cout << "--- TEST 20 ---" << endl; 
		auto plane = new PlanePrimitive("tp", Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0));
		auto pyramid = new PyramidPrimitive("tpyr", 6, 0.2, 0.3);
		Affine3d planeTF = Affine3d::Identity();
		Affine3d pyramidTF = Affine3d::Identity();
		pyramidTF.translation() << 0.0, 0.0, 0.2;
		PrimPrimDistance::distancePrimitivePrimitive(info, plane, planeTF, pyramid, pyramidTF);
		printPrimPrimInfo(info);
		printContactPatch(info.contact_patch);
		Vector2d test(0, 0);
		cout << "distance to patch boundary from test (" << test.transpose() << "): " << endl;
		printPointTestResult(info.contact_patch.testPoint(test));
	}
}

void testNegCapsuleCapsuleDistance() {
	PrimPrimContactInfo info;
	cout << "000000000000 ----- NEG CAPSULE - CAPSULE ----- 000000000000" << endl;
	{
		// test distance between parallel capsules, not axisymmetric, completely contained
		cout << "--- TEST 1 ---" << endl;
		CapsulePrimitive* capsuleA = new NegCapsulePrimitive("tnc", 0.2, 0.4);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tpc", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.1, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, not axisymmetric, completely contained
		// capA translated
		cout << "--- TEST 2 ---" << endl;
		CapsulePrimitive* capsuleA = new NegCapsulePrimitive("tnc", 0.2, 0.4);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tpc", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capATF.translation() << 0.0, 0.1, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, not axisymmetric, completely contained
		// capA rotated
		cout << "--- TEST 3 ---" << endl;
		CapsulePrimitive* capsuleA = new NegCapsulePrimitive("tnc", 0.2, 0.4);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tpc", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capATF.linear() << 0.0, -1.0, 0.0,
						   1.0, 0.0, 0.0,
						   0.0, 0.0, 1.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between parallel capsules, axisymmetric, completely contained
		cout << "--- TEST 4 ---" << endl;
		CapsulePrimitive* capsuleA = new NegCapsulePrimitive("tnc", 0.2, 0.4);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tpc", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, positive larger than negative
		cout << "--- TEST 5 ---" << endl;
		CapsulePrimitive* capsuleA = new NegCapsulePrimitive("tnc", 0.2, 0.4);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tpc", 0.3, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, complete penetration
		cout << "--- TEST 6 ---" << endl;
		CapsulePrimitive* capsuleA = new NegCapsulePrimitive("tnc", 0.2, 0.4);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tpc", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.0, 0.4, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, partial penetration
		cout << "--- TEST 7 ---" << endl;
		CapsulePrimitive* capsuleA = new NegCapsulePrimitive("tnc", 0.2, 0.4);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tpc", 0.1, 0.3);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.1, 0.0, 0.0;
		capBTF.linear() << 0.0, -1.0, 0.0,
						   1.0, 0.0, 0.0,
						   0.0, 0.0, 1.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, partial penetration, asymmetric
		cout << "--- TEST 8 ---" << endl;
		CapsulePrimitive* capsuleA = new NegCapsulePrimitive("tnc", 0.2, 0.4);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tpc", 0.1, 0.3);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.1, 0.1, 0.0;
		capBTF.linear() << 0.0, -1.0, 0.0,
						   1.0, 0.0, 0.0,
						   0.0, 0.0, 1.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, partial penetration, end to end
		cout << "--- TEST 9 ---" << endl;
		CapsulePrimitive* capsuleA = new NegCapsulePrimitive("tnc", 0.2, 0.4);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tpc", 0.1, 0.4);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << 0.4, 0.1, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
	}
	{
		// test distance between capsules, no penetration, filter points
		cout << "--- TEST 10 ---" << endl;
		CapsulePrimitive* capsuleA = new NegCapsulePrimitive("tnc", 0.2, 0.4);
		CapsulePrimitive* capsuleB = new CapsulePrimitive("tpc", 0.1, 0.2);
		Affine3d capATF = Affine3d::Identity();
		Affine3d capBTF = Affine3d::Identity();
		capBTF.translation() << -0.2, 0.0, 0.0;
		PrimPrimDistance::distancePrimitivePrimitive(info, capsuleA, capATF, capsuleB, capBTF);
		printPrimPrimInfo(info);
		cout << "After filtering by max_distance (0.05)" << endl;
		info.filterContactPoints(0.05);
		printPrimPrimInfo(info);
	}
}

void testCircleCapsuleDistance() {
	cout << "000000000000 ----- CIRCLE - CAPSULE ----- 000000000000" << endl;
	{
		// capsule parallel to circle plane, far from circle
		cout << "--- TEST 1 ---" << endl;
		Circle3D circle(/*radius*/ 0.1,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.4);
		Affine3d capTF = Affine3d::Identity();
		capTF.translation() << 0.0, 0.0, 0.2;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule parallel to circle plane, far from circle, end point closest
		cout << "--- TEST 2 ---" << endl;
		Circle3D circle(/*radius*/ 0.1,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.2);
		Affine3d capTF = Affine3d::Identity();
		capTF.translation() << 0.2, 0.0, 0.2;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule parallel to circle plane, far from circle, partial overlap
		// case 1: end pt very close to circle edge
		cout << "--- TEST 3 ---" << endl;
		Circle3D circle(/*radius*/ 0.1,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.2);
		Affine3d capTF = Affine3d::Identity();
		capTF.translation() << 0.199, 0.0, 0.2;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule parallel to circle plane, far from circle, partial overlap
		// case 2: end pt very close to opp circle edge
		cout << "--- TEST 4 ---" << endl;
		Circle3D circle(/*radius*/ 0.1,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.2);
		Affine3d capTF = Affine3d::Identity();
		capTF.translation() << 0.01, 0.0, 0.2;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule parallel to circle plane, far from circle, line pt closest
		cout << "--- TEST 5 ---" << endl;
		Circle3D circle(/*radius*/ 0.1,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.2);
		Affine3d capTF = Affine3d::Identity();
		capTF.translation() << 0.0, 0.1, 0.2;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule perp to circle plane, line pt closest
		cout << "--- TEST 6 ---" << endl;
		Circle3D circle(/*radius*/ 0.2,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.2);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << 0, 0, -1,
						  0, 1, 0,
						  1, 0, 0;
		capTF.translation() << 0.1, 0.0, 0.0;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule perp to circle plane, line end closest, above circle
		cout << "--- TEST 7 ---" << endl;
		Circle3D circle(/*radius*/ 0.2,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.2);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << 0, 0, -1,
						  0, 1, 0,
						  1, 0, 0;
		capTF.translation() << 0.1, 0.0, 0.5;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule perp to circle plane, line end closest, below circle
		cout << "--- TEST 8 ---" << endl;
		Circle3D circle(/*radius*/ 0.2,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.2);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << 0, 0, -1,
						  0, 1, 0,
						  1, 0, 0;
		capTF.translation() << 0.1, 0.0, -0.5;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule perp to circle plane, line pt closest, outside circle
		cout << "--- TEST 9 ---" << endl;
		Circle3D circle(/*radius*/ 0.2,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.2);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << 0, 0, -1,
						  0, 1, 0,
						  1, 0, 0;
		capTF.translation() << 0.3, 0.0, 0.0;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule perp to circle plane passing through origin, symmetric case, line pt
		// closest
		cout << "--- TEST 10 ---" << endl;
		Circle3D circle(/*radius*/ 0.2,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.2);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << 0, 0, -1,
						  0, 1, 0,
						  1, 0, 0;
		capTF.translation() << 0.0, 0.0, 0.0;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule perp to circle plane passing through origin, symmetric case, line end
		// closest
		cout << "--- TEST 11 ---" << endl;
		Circle3D circle(/*radius*/ 0.2,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.2);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << 0, 0, -1,
						  0, 1, 0,
						  1, 0, 0;
		capTF.translation() << 0.0, 0.0, 0.2;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule at angle to circle plane passing through origin, symmetric case, line
		// pts closest
		cout << "--- TEST 12 ---" << endl;
		Circle3D circle(/*radius*/ 0.2,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.4);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << cos(M_PI/4), 0, -sin(M_PI/4),
						  0, 1, 0,
						  sin(M_PI/4), 0, cos(M_PI/4);
		capTF.translation() << 0.0, 0.0, 0.0;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule at angle to circle plane passing through origin, symmetric case, line
		// ends closest
		cout << "--- TEST 13 ---" << endl;
		Circle3D circle(/*radius*/ 0.3,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.2);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << cos(M_PI/4), 0, -sin(M_PI/4),
						  0, 1, 0,
						  sin(M_PI/4), 0, cos(M_PI/4);
		capTF.translation() << 0.0, 0.0, 0.0;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule at angle to circle plane passing through origin, symmetric case, only
		// one line pt closest
		cout << "--- TEST 14 ---" << endl;
		Circle3D circle(/*radius*/ 0.3,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.2);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << cos(M_PI/4), 0, -sin(M_PI/4),
						  0, 1, 0,
						  sin(M_PI/4), 0, cos(M_PI/4);
		capTF.translation() << 0.2, 0.0, 0.2;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule at angle to circle plane passing through origin, asymmetric case, one
		// line end and one line pt closest
		cout << "--- TEST 15 ---" << endl;
		Circle3D circle(/*radius*/ 0.3,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.3);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << cos(M_PI/4), 0, -sin(M_PI/4),
						  0, 1, 0,
						  sin(M_PI/4), 0, cos(M_PI/4);
		capTF.translation() << 0.05, 0.0, 0.05;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule at angle to circle plane not passing through origin, axis perp to
		// vector from origin to point in x-y plane where axis intersects
		// case 1: plane point closest
		cout << "--- TEST 16 ---" << endl;
		Circle3D circle(/*radius*/ 0.3,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.3);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << cos(3*M_PI/8), 0, -sin(3*M_PI/8),
						  0, 1, 0,
						  sin(3*M_PI/8), 0, cos(3*M_PI/8);
		capTF.translation() << 0.0, 0.2, 0.0;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule at angle to circle plane not passing through origin, axis perp to
		// vector from origin to point in x-y plane where axis intersects
		// case 2: line points closest
		cout << "--- TEST 17 ---" << endl;
		Circle3D circle(/*radius*/ 0.3,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 0.6);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << cos(M_PI/8), 0, -sin(M_PI/8),
						  0, 1, 0,
						  sin(M_PI/8), 0, cos(M_PI/8);
		capTF.translation() << 0.0, 0.2, 0.0;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule at angle to circle plane not passing through origin, axis in X-Z plane
		cout << "--- TEST 18 ---" << endl;
		Circle3D circle(/*radius*/ 0.3,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 1.0);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() << cos(M_PI/4), 0, -sin(M_PI/4),
						  0, 1, 0,
						  sin(M_PI/4), 0, cos(M_PI/4);
		capTF.translation() << 0.2, 0.0, 0.0;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}
	{
		// capsule at complex angle to circle plane not passing through origin
		cout << "--- TEST 19 ---" << endl;
		Circle3D circle(/*radius*/ 0.3,
						/*center*/ Vector3d(0.0, 0.0, 0.0),
						/*normal*/ Vector3d(0.0, 0.0, 1.0));
		CapsulePrimitive capsule("tc", 0.1, 1.0);
		Affine3d capTF = Affine3d::Identity();
		capTF.linear() = Matrix3d(AngleAxisd(M_PI/4, Vector3d::UnitZ()) *
						  AngleAxisd(-M_PI/4,  Vector3d::UnitY()));
		capTF.translation() << 0.2, 0.0, 0.0;
		auto info = Circle3DDistance::capsuleDist(capsule, circle, capTF);
		printPrimPrimInfo(info);
	}

}