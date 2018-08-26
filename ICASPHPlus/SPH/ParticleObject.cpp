#include "ParticleObject.h"

RigidBodyParticleObject & RigidBodyParticleObject::operator=(RigidBodyParticleObject & a)
{
	rd = a.rd;
	staticRigidParticleList = a.staticRigidParticleList;
	return *this;
}
