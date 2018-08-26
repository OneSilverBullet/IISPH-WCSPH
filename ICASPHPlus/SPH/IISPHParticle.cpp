#include "IISPHParticle.h"

//对位置以及哈希值进行初始化
void IISPHParticle::Initialization(Vector3f pos, int v, double radius)
{
	position = pos;
	hashIndex = v;
	//刷新半径和支撑半径
    particleRad = radius;
	particleMass = 0.8f * pow(radius*2, 3) * 1000; //这里的1000是静止密度
	particleSupportRad = radius * 4;
	gridCellNum = FunctionKit::CalculateGridCellNum(radius);
}

void IISPHParticle::InitializationPVS(Vector3f pos, int v, double mass, double density)
{
	position = pos;
	hashIndex = v;
	//刷新半径和支撑半径
	particleRad = pow(mass*1.25 / density, double(1.0 / 3.0))*0.5;
	particleMass = mass; //这里的1000是静止密度
	particleSupportRad = particleRad * 4;
	gridCellNum = FunctionKit::CalculateGridCellNum(particleRad);
}

IISPHParticle & IISPHParticle::operator=(IISPHParticle & a)
{
	position = a.position;
	velocity = a.velocity;
	acceleration = a.acceleration;
	particleRad = a.particleRad;
	particleSupportRad = a.particleSupportRad;
	particleMass = a.particleMass;
	density = a.density;
	viscosity = a.viscosity;
	pressureAcceleration = a.pressureAcceleration;
	hashIndex = a.hashIndex;

	//IISPH属性初始化
	aii = a.aii;
	dii = a.dii;
	dij_pj = a.dij_pj;
	density_adv = a.density_adv;
	pressure = a.pressure;
	lastPressure = a.lastPressure;
	pressureAcceleration = a.pressureAcceleration;
	gridCellNum = a.gridCellNum;

	//无限自适应的扩展属性初始化
	distance = a.density;
	mark = a.mark;
	mopt = a.mopt;
	distance = a.distance;
	mark = a.mark;
	O = a.O;
	pt = a.pt;
	lastDensity = a.lastDensity;
	lastVelocity = a.lastVelocity;
	β = a.β;


	//属性修正
	//拷贝该粒子的边界粒子邻居
	fluidNeighbors.resize(a.fluidNeighbors.size());
	for (int i=0; i<a.fluidNeighbors.size(); i++)
	{
		//注意：这里应该不是深拷贝，而是拷贝的指针
		fluidNeighbors[i] = a.fluidNeighbors[i];
	}

	//拷贝该粒子的边界粒子邻居
	boundaryNeighbors.resize(a.boundaryNeighbors.size());
	for (int i = 0; i < a.boundaryNeighbors.size(); i++)
	{
		boundaryNeighbors[i] = a.boundaryNeighbors[i];
	}

	//拷贝该粒子的动态刚体粒子邻居
	dynamicRigidNeighbors.resize(a.dynamicRigidNeighbors.size());
	for (int i=0; i< a.dynamicRigidNeighbors.size(); i++)
	{
		dynamicRigidNeighbors[i] = a.dynamicRigidNeighbors[i];
	}

	return *this;
}

//在该函数调用的时候，关于粒子的所有的力都已经计算完毕
void IISPHParticle::Intigration(Boundary bound, double timeStep)
{
	//在重力的基础上增加压力加速度(黏度加速度已经在前面计算好了)
	velocity += pressureAcceleration * timeStep;
	velocity = (1 - β)*velocity + β*lastVelocity;
	position += velocity*timeStep;
	lastVelocity = velocity;
	//将粒子的位置进行限制
	ParticlePositionCorrection(bound);
}

void IISPHParticle::ParticlePositionCorrection(Boundary boundary)
{
	int Flag = 0;
	if (position.x < boundary.boundaryOriginX - boundary.offset / 2.0)
	{
		Flag = 1;
	}
	if (position.y < boundary.boundaryOriginY - boundary.offset / 2.0)
	{
		Flag = 1;
	}
	if (position.z < boundary.boundaryOriginZ - boundary.offset / 2.0)
	{
		Flag = 1;
	}
	if (position.x > boundary.boundaryWidth + boundary.boundaryOriginX + boundary.offset / 2.0)
	{
		Flag = 1;
	}
	if (position.y > boundary.boundaryLength + boundary.boundaryOriginY + boundary.offset / 2.0)
	{
		Flag = 1;
	}
	if (position.z > boundary.boundaryHeight + boundary.boundaryOriginZ + boundary.offset / 2.0)
	{
		Flag = 1;
	}
	//如果当前粒子超过相应的范围，那么就对其速度归0
	if (Flag)
	{
		velocity.Zero();
	}
}

