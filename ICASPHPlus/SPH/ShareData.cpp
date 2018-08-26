#include "ShareData.h"

Vector3f SharedData::boundaryOrigin = Vector3f(0.0, 0.0, 0.0);
Vector3f SharedData::boundaryWLH = Vector3f(2, 4, 2);
double SharedData::offset = 0.6f;
double SharedData::timeStep = 0.001f;
double SharedData::particleRad = 0.025f;
double SharedData::particleSupportRad = 0.1f;
double SharedData::particleMass = 0.1f;
double SharedData::restDensity = 1000.0f;
int    SharedData::particleNum;
Vector3f SharedData::fluidWLH = Vector3f(1.0, 1.0, 1.0);
Vector3f SharedData::fluidOrigin = Vector3f(0.25, 0.25, 0.0);
double SharedData::cellLength = 0.1f;
int SharedData::hashCellNum = 0.0f;
double SharedData::stiffnesss = 50000;
int SharedData::gamma = 7;
double SharedData::viscosityConstant = 0.03f;
double SharedData::gravity = -9.81;

void SharedData::SetBoundaryOrigin(const Vector3f & a)
{
	boundaryOrigin = a;
}

Vector3f SharedData::GetBoundaryOrigin()
{
	return boundaryOrigin;
}

void SharedData::SetBoundaryWLH(const Vector3f & a)
{
	boundaryWLH = a;
}

Vector3f SharedData::GetBoundaryWLH()
{
	return boundaryWLH;
}

void SharedData::SetOffset(const double a)
{
	offset = a;
}

double SharedData::GetOffset()
{
	return offset;
}

void SharedData::SetTimeStep(const double a)
{
	timeStep = a;
}

double SharedData::GetTimeStep()
{
	return timeStep;
}

void SharedData::SetParticleRad(const double a)
{
	particleRad = a;
}

double SharedData::GetParticleRad()
{
	return particleRad;
}

void SharedData::SetParticleSupportRad(const double a)
{
	particleSupportRad = a;
}

double SharedData::GetParticleSupportRad()
{
	return particleRad;
}

void SharedData::SetParticleMass(const double a)
{
	particleMass = a;
}

double SharedData::GetParticleMass()
{
	return particleMass;
}

void SharedData::SetRestDensity(const double a)
{
	restDensity = a;
}

double SharedData::GetRestDensity()
{
	return restDensity;
}

void SharedData::SetParticleNum(const int a)
{
	particleNum = a;
}

int SharedData::GetParticleNum()
{
	return particleNum;
}

void SharedData::SetFluidWLH(const Vector3f & a)
{
	fluidWLH = a;
}

Vector3f SharedData::GetFluidWLH()
{
	return fluidWLH;
}

void SharedData::SetFluidOrigin(const Vector3f & a)
{
	fluidOrigin = a;
}

Vector3f SharedData::GetFluidOrigin()
{
	return fluidOrigin;
}

void SharedData::SetCellLength(const double a)
{
	cellLength = a;
}

double SharedData::GetCellLength()
{
	return cellLength;
}

void SharedData::SetHashCellNum(const int a)
{
	hashCellNum = a;
}

int SharedData::GetHashCellNum()
{
	return hashCellNum;
}

void SharedData::SetStiffnesss(const double a)
{
	stiffnesss = a;
}

double SharedData::GetStiffnesss()
{
	return stiffnesss;
}

void SharedData::SetGamma(const int a)
{
	gamma = a;
}

int SharedData::GetGamma()
{
	return gamma;
}

void SharedData::SetViscosityConstant(const double a)
{
	viscosityConstant = a;
}

double SharedData::GetViscosityConstant()
{
	return viscosityConstant;
}

void SharedData::SetGravity(const double a)
{
	gravity = a;
}

double SharedData::GetGravity()
{
	return gravity;
}
