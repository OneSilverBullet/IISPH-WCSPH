#include "IISPHParticle.h"

//��λ���Լ���ϣֵ���г�ʼ��
void IISPHParticle::Initialization(Vector3f pos, int v, double radius)
{
	position = pos;
	hashIndex = v;
	//ˢ�°뾶��֧�Ű뾶
    particleRad = radius;
	particleMass = 0.8f * pow(radius*2, 3) * 1000; //�����1000�Ǿ�ֹ�ܶ�
	particleSupportRad = radius * 4;
	gridCellNum = FunctionKit::CalculateGridCellNum(radius);
}

void IISPHParticle::InitializationPVS(Vector3f pos, int v, double mass, double density)
{
	position = pos;
	hashIndex = v;
	//ˢ�°뾶��֧�Ű뾶
	particleRad = pow(mass*1.25 / density, double(1.0 / 3.0))*0.5;
	particleMass = mass; //�����1000�Ǿ�ֹ�ܶ�
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

	//IISPH���Գ�ʼ��
	aii = a.aii;
	dii = a.dii;
	dij_pj = a.dij_pj;
	density_adv = a.density_adv;
	pressure = a.pressure;
	lastPressure = a.lastPressure;
	pressureAcceleration = a.pressureAcceleration;
	gridCellNum = a.gridCellNum;

	//��������Ӧ����չ���Գ�ʼ��
	distance = a.density;
	mark = a.mark;
	mopt = a.mopt;
	distance = a.distance;
	mark = a.mark;
	O = a.O;
	pt = a.pt;
	lastDensity = a.lastDensity;
	lastVelocity = a.lastVelocity;
	�� = a.��;


	//��������
	//���������ӵı߽������ھ�
	fluidNeighbors.resize(a.fluidNeighbors.size());
	for (int i=0; i<a.fluidNeighbors.size(); i++)
	{
		//ע�⣺����Ӧ�ò�����������ǿ�����ָ��
		fluidNeighbors[i] = a.fluidNeighbors[i];
	}

	//���������ӵı߽������ھ�
	boundaryNeighbors.resize(a.boundaryNeighbors.size());
	for (int i = 0; i < a.boundaryNeighbors.size(); i++)
	{
		boundaryNeighbors[i] = a.boundaryNeighbors[i];
	}

	//���������ӵĶ�̬���������ھ�
	dynamicRigidNeighbors.resize(a.dynamicRigidNeighbors.size());
	for (int i=0; i< a.dynamicRigidNeighbors.size(); i++)
	{
		dynamicRigidNeighbors[i] = a.dynamicRigidNeighbors[i];
	}

	return *this;
}

//�ڸú������õ�ʱ�򣬹������ӵ����е������Ѿ��������
void IISPHParticle::Intigration(Boundary bound, double timeStep)
{
	//�������Ļ���������ѹ�����ٶ�(�ȼ��ٶ��Ѿ���ǰ��������)
	velocity += pressureAcceleration * timeStep;
	velocity = (1 - ��)*velocity + ��*lastVelocity;
	position += velocity*timeStep;
	lastVelocity = velocity;
	//�����ӵ�λ�ý�������
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
	//�����ǰ���ӳ�����Ӧ�ķ�Χ����ô�Ͷ����ٶȹ�0
	if (Flag)
	{
		velocity.Zero();
	}
}

