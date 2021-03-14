// ContactDetection.cpp

#include <iostream>
#include "COPSimulator.h"

using namespace Eigen;

namespace Sai2COPSim {

// internal functions TODO: considering moving to a helper class
Eigen::Affine3d getPrimitiveTransformInWorld(
	const Primitive* prim,
	ArticulatedRigidBodyManager& arb_manager
) {
	auto linkTransform = Eigen::Affine3d::Identity();
	if(!prim->_is_static) {
		assert(!prim->_link_name.empty());
		assert(!prim->_articulated_body_name.empty());
		auto arb = arb_manager.getBody(prim->_articulated_body_name);
		arb->_model->transformInWorld(linkTransform, prim->_link_name);
	}
	return linkTransform*prim->_transform_in_link;
}

// returns the index of the island in the map if prim->_arb_name exists in an island
// else returns -1
int getPrimitiveARBContactIsland(const Primitive* prim, WorldContactMap& map) {
	if(prim->_is_static) return -1;
	assert(!prim->_articulated_body_name.empty());
	std::string arb_name = prim->_articulated_body_name;
	for(uint i = 0; i < map._islands.size(); i++) {
		if(map._islands[i]._is_merged) continue; // skip merged islands
		if(map._islands[i]._articulated_bodies.find(arb_name) != map._islands[i]._articulated_bodies.end()) {
			return i;
		}
	}
	return -1;
}

// works with the geometry manager to compute the WorldContactMap
void COPSimulator::computeWorldContactMap() {
	// reset existing world contact map
	// TODO: implement sleep on objects, so that we don't have to reset the full world
	// contact map
	_contact_map.clear();
	_max_penetration_current = -1;

	//TODO: Implement broad phase collision detection with sphere hierarchy in GeometryManager

	// narrow phase collision check primitive to primitive
	// TODO: parallelize prim - prim distance check
	for(uint i = 0; i < _geom_manager._primitives.size(); i++) {
		auto prim_i = _geom_manager._primitives[i];
		// get transform for prim_i
		// TODO: use a link transform to world cache
		auto prim_i_tf = getPrimitiveTransformInWorld(prim_i, _arb_manager);
		for(uint j = i+1; j < _geom_manager._primitives.size(); j++) {
			auto prim_j = _geom_manager._primitives[j];

			// TODO: flag for disallowing collision between parent and child links
			// within a sphere around the joint
			// e.g. collisions between two links on the same ARB,
			//  collision with an ARB for which collision is disabled, etc.

			// check for collision between two static meshes
			if(prim_i->_is_static && prim_j->_is_static) continue;

			// check for collision between two meshes on the same link
			if(!prim_i->_is_static && !prim_j->_is_static &&
				prim_i->_articulated_body_name.compare(prim_j->_articulated_body_name) == 0
			) {
				// TODO: flag for disallowing self-collisions
				if(prim_i->_link_name.compare(prim_j->_link_name) == 0) continue;
			}

			// get link transform
			auto prim_j_tf = getPrimitiveTransformInWorld(prim_j, _arb_manager);
			PrimPrimContactInfo* ppinfo = _geom_manager._prim_prim_distances[j][i];
			PrimPrimDistance::distancePrimitivePrimitive(
				*ppinfo,
				prim_i, prim_i_tf, prim_j, prim_j_tf);
			// std::cout << prim_i->_name << " " << prim_j->_name << std::endl;
			// std::cout << ppinfo.min_distance << std::endl;

			// normal info returned is from prim_i to prim_j
			// however, if either is static, then the direction must be from the static
			// primitive towards the ARB
			if(prim_j->_is_static) {
				ppinfo->flipNormal();
			}

			// check if in contact
			if(ppinfo->min_distance < COPAlgorithmicConstants::GEOMETRIC_CONTACT_DISTANCE_THRESHOLD) {
				_max_penetration_current = max(_max_penetration_current, -ppinfo->min_distance);
				// create ContactPrimitivePair
				ContactPrimitivePair contact_pair;
				contact_pair.primA = prim_i;
				contact_pair.primA_in_world = prim_i_tf;
				contact_pair.primB = prim_j;
				contact_pair.primB_in_world = prim_j_tf;
				contact_pair.info = ppinfo;

				// TODO: in a parallel implementation, contact island assignment has to be serialized
				// so remove to a separate loop over all generated ContactPrimitivePairs

				// check if the associated articulated bodies are already in a contact island
				// if they are, we should just add it to that island.
				// there is also a chance that they are in two separate contact islands, in which
				// case the islands need to be merged.

				int contact_island_prim_i = getPrimitiveARBContactIsland(prim_i, _contact_map);
				int contact_island_prim_j = getPrimitiveARBContactIsland(prim_j, _contact_map);

				// TODO: cache contact island info for prim_i
				// case 1: neither primitive exists in an island already, so create one
				if(contact_island_prim_i == -1 && contact_island_prim_j == -1) {
					ContactIsland island;
					island._contact_prim_pairs.push_back(contact_pair);
					if(!prim_i->_is_static) island._articulated_bodies.insert(prim_i->_articulated_body_name);
					if(!prim_j->_is_static) island._articulated_bodies.insert(prim_j->_articulated_body_name);
					_contact_map._islands.push_back(island);
				} else if(contact_island_prim_i > -1 && contact_island_prim_j == -1) {
					// case 2: add to contact island where prim_i->_arb_name already exists
					_contact_map._islands[contact_island_prim_i]._contact_prim_pairs.push_back(contact_pair);
					if(!prim_j->_is_static) {
						_contact_map._islands[contact_island_prim_i]._articulated_bodies.insert(prim_j->_articulated_body_name);
					}
				} else if(contact_island_prim_i == -1 && contact_island_prim_j > -1) {
					// case 3: add to contact island where prim_j->_arb_name already exists
					_contact_map._islands[contact_island_prim_j]._contact_prim_pairs.push_back(contact_pair);
					if(!prim_i->_is_static) {
						_contact_map._islands[contact_island_prim_j]._articulated_bodies.insert(prim_i->_articulated_body_name);
					}
				} else {
					// case 4: both primitive bodies are in separate islands. Merge and add in
					// island for prim_i
					_contact_map._islands[contact_island_prim_i].merge(_contact_map._islands[contact_island_prim_j]);
					_contact_map._islands[contact_island_prim_j]._is_merged = true;
					_contact_map._islands[contact_island_prim_i]._contact_prim_pairs.push_back(contact_pair);
				}
			}
		}
	}
}

}
