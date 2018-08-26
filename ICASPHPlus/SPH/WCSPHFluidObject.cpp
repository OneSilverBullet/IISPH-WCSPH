#include "WCSPHFluidObject.h"

void WCSPHFluidObject::Initialise()
{
	SetParticleRad();  //���µ�ǰ�⻬�˵������Ϣ
	ComputeFluidWLH(); //����Ĭ�ϲ�������õ�WLH
	ComputeParticleNum(fluidWLH.GetX(), fluidWLH.GetY(), fluidWLH.GetZ()); //���ݼ���õ������峤��߼���õ������������ҽ�������
	InitialiseParticlePosition(); //����������Ϣ��������λ�ó�ʼ��
	InitialiseParticleDensity();
	InitialiseParticleP();
	neighborList.resize(particleNum); //���ھӱ������չ
	fluidBoundaryNeighborList.resize(particleNum); //���������ӵı߽������ھӽ�����չ
												   //�Ա߽����ӽ��г�ʼ��
	InitializeBoundary();
	//��ȫ�ֱ������г�ʼ��
	Vector3f origin(originX, originY, originZ);
	SharedData::SetFluidOrigin(origin);
	Vector3f wlh(fluidWitdth, fluidLength, fluidHeight);
	SharedData::SetBoundaryWLH(wlh);
	SharedData::SetRestDensity(restDensity);
	SharedData::SetParticleRad(particleRad);
	SharedData::SetParticleSupportRad(particleSupportRad);
	SharedData::SetTimeStep(timeStep);
}

void WCSPHFluidObject::Initialise(Vector3f WLH, Vector3f fluidOrigin, double rD, double pR, double time, Boundary bound)
{
	//���ݳ�ʼ��
	fluidWitdth = WLH.GetX();
	fluidLength = WLH.GetY();
	fluidHeight = WLH.GetZ();
	originX = fluidOrigin.GetX();
	originY = fluidOrigin.GetY();
	originZ = fluidOrigin.GetZ();
	restDensity = rD;
	SetParticleRad(rD);
	timeStep = time;
	boundary = bound;
	//���ݻ�����Ϣ�����ٴμ���
	ComputeFluidWLH(); //����Ĭ�ϲ�������õ�WLH
	ComputeParticleNum(fluidWLH.GetX(), fluidWLH.GetY(), fluidWLH.GetZ()); //���ݼ���õ������峤��߼���õ������������ҽ�������
	InitialiseParticlePosition(); //����������Ϣ��������λ�ó�ʼ��
	InitialiseParticleDensity();
	InitialiseParticleP();
	neighborList.resize(particleNum); //��ʼ��ˮ����ھӱ�
									  //�Ա߽����ӽ��г�ʼ��
	InitializeBoundary();
	//�����SharedData���г�ʼ������
	SharedData::SetFluidOrigin(fluidOrigin);
	SharedData::SetBoundaryWLH(WLH);
	SharedData::SetRestDensity(rD);
	SharedData::SetParticleRad(pR);
	SharedData::SetParticleSupportRad(pR * 4);
	SharedData::SetTimeStep(time);
}

//���������Ҫ�����в����������֮�����
void WCSPHFluidObject::ComputeFluidWLH()
{
	int w = floor(double(fluidWitdth / (2 * particleRad))); //�õ�x�����ϵĳ���
	int l = floor(double(fluidLength / (2 * particleRad))); //�õ�y�����ϵĳ���
	int h = floor(double(fluidHeight / (2 * particleRad))); //�õ�z�����ϵĳ���
	fluidWLH.SetValue(w, l, h);                             //��ˮ��ĳ���߷�����г�ʼ��
}

void WCSPHFluidObject::ComputeFluidWLH(double fluidW, double fluidL, double fluidH, double r)
{
	int w = floor(double(fluidW / (2 * r))); //�õ�x�����ϵĳ���
	int l = floor(double(fluidL / (2 * r))); //�õ�y�����ϵĳ���
	int h = floor(double(fluidH / (2 * r))); //�õ�z�����ϵĳ���
	fluidWLH.SetValue(w, l, h);
}

void WCSPHFluidObject::InitialiseParticlePosition()
{
	//��ʼ�����ӵ�λ��
	//�����б��е����
	int particleIndex = 0;
	//���ӵ�ֱ��
	double particleLength = particleRad * 2;
	for (int i = 0; i < fluidWLH.GetX(); i++)
	{
		for (int j = 0; j < fluidWLH.GetY(); j++)
		{
			for (int k = 0; k < fluidWLH.GetZ(); k++)
			{
				Vector3f position(0, 0, 0); //��ʼ����ǰ���������ӵ�λ��
				position.SetX(originX + particleRad + particleLength*i); //���õ�ǰ���ӵ�x����
				position.SetY(originY + particleRad + particleLength*j); //���õ�ǰ���ӵ�y����
				position.SetZ(originZ + particleRad + particleLength*k); //���õ�ǰ���ӵ�z����
				int hashValue = FunctionKit::PositionMapHash(position,
					SharedData::GetCellLength(),
					SharedData::GetHashCellNum(),
					boundary.offset); //�õ���ǰposition����Ӧ�Ĺ�ϣֵ
				WCSPHParticle* newParticle = new WCSPHParticle(); //��ʼ��һ���µ�����
				newParticle->Initialization(position, hashValue); //�������ӽ��г�ʼ��
																  //����ʼ���õ�������װ�������б���
				particleList[particleIndex] = newParticle;
				particleIndex++;
			}
		}
	}
}

void WCSPHFluidObject::InitialiseParticleDensity()
{
	for (int i = 0; i < particleList.size(); i++)
	{
		particleList[i]->density = 1000.0;
	}
}

void WCSPHFluidObject::InitialiseParticleP()
{
	for (int i = 0; i < particleList.size(); i++)
	{
		particleList[i]->P = 0.0f;
	}
}

void WCSPHFluidObject::UpdateParticlePosition()
{
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static) 
		for (int i = 0; i < particleList.size(); i++)
		{
			//ˢ��ÿ�����ӵ�λ�ú��ٶ�
			Vector3f acc = particleList[i]->acceleration;
			// v = v + a * t 
			particleList[i]->velocity = particleList[i]->velocity + acc*timeStep;
			//p = p + v * t
			particleList[i]->position = particleList[i]->position + particleList[i]->velocity*timeStep;
			//�����ӵ��ٶȺ�λ�ý��о���				
			ParticleCorrection(particleList[i]);																						 																																									 //ParticleCorrection(particleList[i]);
		}
	}
}

void WCSPHFluidObject::UpdateParticleHashValue()
{
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static) 
		for (int i = 0; i < particleList.size(); i++)
		{
			double cellLength = SharedData::GetCellLength(); //�õ�����Ĵ�С
			int hashCellNum = SharedData::GetHashCellNum();  //�õ����й�ϣ��������
			int hashValue = FunctionKit::PositionMapHash(particleList[i]->position, cellLength, hashCellNum, boundary.offset);
			particleList[i]->hashIndex = hashValue;
		}
	}
}

//����boundary�ķ�Χ�������ӵ���ײ���
//����������Ϊ����߽磬�ᵼ�µ���������������ڴ������Ϊ�˱����������
//ʹ�ñ߽��⣬һ�������������ô�ͽ����ӵ��ٶ�ǿ������Ϊ0
void WCSPHFluidObject::ParticleCorrection(WCSPHParticle* a)
{
	int Flag = 0;
	if (a->position.x < boundary.boundaryOriginX - boundary.offset / 2.0)
	{
		Flag = 1;
	}
	if (a->position.y < boundary.boundaryOriginY - boundary.offset / 2.0)
	{
		Flag = 1;
	}
	if (a->position.z < boundary.boundaryOriginZ - boundary.offset / 2.0)
	{
		Flag = 1;
	}
	if (a->position.x > boundary.boundaryWidth + boundary.boundaryOriginX + boundary.offset / 2.0)
	{
		Flag = 1;
	}
	if (a->position.y > boundary.boundaryLength + boundary.boundaryOriginY + boundary.offset / 2.0)
	{
		Flag = 1;
	}
	if (a->position.z > boundary.boundaryHeight + boundary.boundaryOriginZ + boundary.offset / 2.0)
	{
		Flag = 1;
	}
	//�����ǰ���ӳ�����Ӧ�ķ�Χ����ô�Ͷ����ٶȹ�0
	if (Flag)
	{
		a->velocity.Zero();
	}
}

void WCSPHFluidObject::ComputeParticleNum(int x, int y, int z)
{
	particleNum = x*y*z;               //����x,y,z�᲻ͬ�ķ�����������������
	particleList.resize(particleNum);  //����������������ˢ�������б�����
}

void WCSPHFluidObject::SetParticleRad()
{
	Poly6Kernel::setRadius(particleSupportRad);
	SpikyKernel::setRadius(particleSupportRad);
	CubicKernel::setRadius(particleSupportRad);
}

void WCSPHFluidObject::SetParticleRad(double r)
{
	particleRad = r;
	particleSupportRad = 4 * r;
	Poly6Kernel::setRadius(particleSupportRad);
	SpikyKernel::setRadius(particleSupportRad);
	CubicKernel::setRadius(particleSupportRad);
}

//����ھӱ������е���Ϣ
void WCSPHFluidObject::ClearNeighborList()
{
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static) 
		//ע����������ǰ�ھ��б�����б�
		for (int i = 0; i < neighborList.size(); i++)
		{
			neighborList[i].clear();
		}
	}

}

//��������߽���ھ��б�
void WCSPHFluidObject::ClearFluidBoundaryNeighborList()
{
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static) 
		for (int i = 0; i < fluidBoundaryNeighborList.size(); i++)
		{
			fluidBoundaryNeighborList[i].clear();
		}
	}
}

//���������ִ�У������ڣ�
//1 ����ģ�͵�boundary�Ѿ���SPH Computation��ֵ���
//2 ����ģ�͵İ뾶��ֵ��Kernel����
void WCSPHFluidObject::InitializeBoundary()
{
	//���ݱ߽����ӣ������Ӧ�ı߽����ӵ�λ������
	vector<Vector3f> boundaryParticlePosition = boundary.initBoundaryData();
	//���ݼ���õ��ı߽����ӵĸ����Ա߽��ھ������б��������
	boundaryNeighborList.resize(boundaryParticlePosition.size());
	//��ʼ����̬���������
	StaticRigidObject rb;
	//ʹ�þ�̬������߽����ӵ�λ�öԾ�̬���ӱ߽���г�ʼ��
	AddRigidBodyObject(rb, boundaryParticlePosition.size(), boundaryParticlePosition);

	//������������SPHComputation������
	//���ݾ�̬���ӵ�λ�ã�������ӳ�䵽������
	//����ÿһ����̬�߽����ӵ�BoundaryPSI
}

//����λ���������ݣ���ÿһ�����ӽ���λ�õĸ�ֵ
void WCSPHFluidObject::AddRigidBodyObject(StaticRigidObject& rbo, const unsigned int numBoundaryParticles, vector<Vector3f> boundaryParticles)
{
	//�ӹ�ϣ�����б��еõ������ؼ��ľ�ֵ̬,���õĶ��Ǿ�ֵ̬
	int cellNum = SharedData::GetHashCellNum();
	double cellLength = SharedData::GetCellLength();

	std::cout << "��Ӹ������" << std::endl;
	//RigidBodyParticleObject *rb = new RigidBodyParticleObject(); //�´�����Ӧ�ĸ������Ӷ���
	//boundaryObj = rb;                                            //�������Ӷ��󿽱�������ı߽����

	boundaryObj.staticRigidParticleList.resize(numBoundaryParticles);    //��ʼ����ǰ��̬���������б�Ĵ�С													
	for (int i = 0; i < (int)numBoundaryParticles; i++)
	{
		boundaryObj.staticRigidParticleList[i] = new StaticRigidParticle();
	}
	for (int i = 0; i < (int)numBoundaryParticles; i++)
	{
		boundaryObj.staticRigidParticleList[i]->x0 = boundaryParticles[i];
		boundaryObj.staticRigidParticleList[i]->position = boundaryParticles[i];
		boundaryObj.staticRigidParticleList[i]->velocity.Zero();
		boundaryObj.staticRigidParticleList[i]->m_f.Zero();
		//���ݵ�ǰλ�ü���õ���Ӧ���ӵĹ�ϣֵ
		boundaryObj.staticRigidParticleList[i]->hashValue =
			FunctionKit::PositionMapHash(boundaryParticles[i], cellLength,
				cellNum, boundary.offset);
	}
	//�Ա߽����ĸ������Խ��и�ֵ
	boundaryObj.rd = rbo;
	std::cout << "�������������" << std::endl;
}


//�ú���һ��Ҫ��
void WCSPHFluidObject::ComputeBoundaryPsi(RigidBodyParticleObject& bound)
{
	const double density0 = restDensity; //�õ������ʼ�ܶ�
										 //�õ���̬���Ӹ���
	const unsigned int numBoundaryParticles = boundaryObj.GetStaticParticleNum();
	//�õ��߽����ӵ�֧�Ű뾶
	double supportRad = boundary.pRad * 4;
	//������ǰ�ı߽�����
	for (int i = 0; i < (int)numBoundaryParticles; i++)
	{
		double delta = CubicKernelOne::W_zero(supportRad);
		//�õ���Ӧ������bound���еľ�̬��������
		StaticRigidParticle* a = boundaryObj.staticRigidParticleList[i];
		for (unsigned int j = 0; j < boundaryNeighborList[i].size(); j++)
		{
			//�õ���Ӧ�ٽ����ӵ�ǰ���ٽ���
			StaticRigidParticle* b = boundaryNeighborList[i][j];
			//����õ���ǰ������������õ��Ĳ�ֵ
			Vector3f positionSub = a->position - b->position;
			delta += CubicKernelOne::W(supportRad, positionSub);
		}
		const double volume = 1.0 / delta;
		//�Ե�ǰ��̬���ӵ�boundaryPsi���и�ֵ
		boundaryObj.staticRigidParticleList[i]->boundaryPsi = density0 * volume;
		//std::cout << "  BoundaryPSI : " << boundaryObj.staticRigidParticleList[i]->boundaryPsi << std::endl;
	}
}
