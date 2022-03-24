// Composite1PkN.cpp

#include "Composite1PkN.h"

namespace Sai2COPSim {

Composite1PkN::Composite1PkN(const std::string& name, Primitive* positivePrimitive) {
	if(name.empty()) throw(std::runtime_error("Name cannot be empty"));
	_name = name;
	_type = GeometryType::Composite1PkN;
	_positivePrimitive = positivePrimitive;
}

}