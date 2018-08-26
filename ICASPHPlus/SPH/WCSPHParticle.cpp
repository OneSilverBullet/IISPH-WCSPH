#include "WCSPHParticle.h"

void WCSPHParticle::Initialization(Vector3f p, int v)
{
	density = DENSITY0;
	P = 0;
	position = p;
	velocity = new Vector3f(0, 0, 0);
	acceleration = new Vector3f(0, 0, 0);
	viscosity.Zero();
	pressure.Zero();
	hashIndex = v;
}

WCSPHParticle & WCSPHParticle::operator=(WCSPHParticle & a)
{
	this->position = a.position;
	this->velocity = a.velocity;
	this->acceleration = a.acceleration;
	this->density = a.density;
	this->hashIndex = a.hashIndex;
	this->P = a.P;
	return *this;
}
