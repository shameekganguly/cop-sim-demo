// WorldContactMap.h

#ifndef WORLDCONTACTMAP_H
#define WORLDCONTACTMAP_H

#include <string>
#include <list>
#include <unordered_set>

#include "geometry/PrimPrimContactInfo.h"

namespace Sai2COPSim {

class Primitive;

struct ContactPrimitivePair {
	Primitive* primA;
	Eigen::Affine3d primA_in_world;
	Primitive* primB;
	Eigen::Affine3d primB_in_world;
	PrimPrimContactInfo* info;
};

class ContactIsland {
public:
	typedef std::list<ContactPrimitivePair> ContactList;

public:
	ContactIsland(): _is_merged(false) {
		// nothing to do right now
	}

	// merge another contact island into this one
	void merge(ContactIsland& other);

public:
	//TODO: think about a more efficient data structure. 
	// Maybe a graph, where nodes are bodies and the edges contain the contact data?
	// But computing the island would still require a BFS. 
	// Also, this data structure is torn down and recomputed every time step right now
	// So, its not particularly efficient. But then again, if a particular contact breaks,
	// computing whether or not the rest of the graph is divided into two is not an easy
	// problem either

	// list of ContactPrimitivePairs in this island
	ContactList _contact_prim_pairs;

	// Set of bodies in this contact island.
	// Used to find if two ContactIslands share a body, in which case they should be treated
	// as one.
	std::unordered_set<std::string> _articulated_bodies;

	// has this contact island already been merged into another one?
	bool _is_merged;
};

class WorldContactMap {
public:
	WorldContactMap() {
		// nothing to do right now
	}

	// clear this world map
	void clear() { _islands.clear(); }

	// does world have any primitive pairs penetrating where they should not
	// input param: threshold, a (non-positive) number giving the allowed minimum
	// distance between two primitives before they are considered to be penetrating
	ContactIsland::ContactList penetratingPrimPairs(double threshold = 0.0) {
		// TODO: implement
		return ContactIsland::ContactList();
	}

	// TODO: provide const interfaces to retrieve islands

	// TODO: some more introspection functions. 
	// e.g. getting the contact force on a particular link


	// TODO: some visualization

public:
	std::vector<ContactIsland> _islands;
};

}

#endif//