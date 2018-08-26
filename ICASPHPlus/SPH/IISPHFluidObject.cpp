#include "IISPHFluidObject.h"

void IISPHFluidObject::ClearNeighborList()
{
    #pragma omp parallel default(shared)
	{
       #pragma omp for schedule(static) 
		//ע����������ǰ�ھ��б�����б�
		for (int i = 0; i < IISPH_ParticleNum; i++)
		{
			particleList[i]->fluidNeighbors.clear();
		}
	}
}

void IISPHFluidObject::ClearFluidBoundaryNeighborList()
{
    #pragma omp parallel default(shared)
	{
       #pragma omp for schedule(static) 
		for (int i = 0; i < IISPH_ParticleNum; i++)
		{
			particleList[i]->boundaryNeighbors.clear();
		}
	}
}

void IISPHFluidObject::ClearDynamicRigidNeighborList()
{
     #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static) 
		for (int i = 0; i < IISPH_ParticleNum; i++)
		{
			particleList[i]->dynamicRigidNeighbors.clear();
		}
	}
}

void IISPHFluidObject::InitializeBoundary()
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

void IISPHFluidObject::AddRigidBodyObject(StaticRigidObject & rbo, const unsigned int numBoundaryParticles, vector<Vector3f> boundaryParticles)
{
	//�ӹ�ϣ�����б��еõ������ؼ��ľ�ֵ̬,���õĶ��Ǿ�ֵ̬
	int cellNum = SharedData::GetHashCellNum();
	double cellLength = SharedData::GetCellLength();

	std::cout << "��Ӹ������" << std::endl;
	boundaryObj.staticRigidParticleList.clear();
	boundaryObj.staticRigidParticleList.resize(numBoundaryParticles);    //��ʼ����ǰ��̬���������б�Ĵ�С
	//������г�ʼ��
	//ÿһ��ָ�붼ָ��һ����̬��������
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
}

void IISPHFluidObject::ComputeBoundaryPsi(RigidBodyParticleObject & bound)
{
	boundaryNumber = boundary.boundaryNumber;
	double radius = boundary.pRad * 4;
	const double density0 = IISPH_RestDensity; //�õ������ʼ�ܶ�
	//�õ���̬���Ӹ���
	const unsigned int numBoundaryParticles = boundaryObj.GetStaticParticleNum();
	//������ǰ�ı߽�����
	for (int i = 0; i < (int)numBoundaryParticles; i++)
	{
		double delta = CubicKernelOne::W_zero(radius);
		//�õ���Ӧ������bound���еľ�̬��������
		StaticRigidParticle* a = boundaryObj.staticRigidParticleList[i];
		for (unsigned int j = 0; j < boundaryNeighborList[i].size(); j++)
		{
			//�õ���Ӧ�ٽ����ӵ�ǰ���ٽ���
			StaticRigidParticle* b = boundaryNeighborList[i][j];
			//����õ���ǰ������������õ��Ĳ�ֵ
			Vector3f positionSub = a->position - b->position;
			delta += CubicKernelOne::W(radius, positionSub);
		}
		const double volume = 1.0 / delta;
		//�Ե�ǰ��̬���ӵ�boundaryPsi���и�ֵ
		boundaryObj.staticRigidParticleList[i]->boundaryPsi = density0 * volume;
		std::cout << "  BoundaryPSI : " << boundaryObj.staticRigidParticleList[i]->boundaryPsi << std::endl;
	}
}

//���㵱ǰ���ж�̬���ӵ�boundaryPsi
void IISPHFluidObject::ComputeBoundaryPsi(DynamicRigidParticleObject& cube)
{
	boundaryNumber = boundary.boundaryNumber;
	double radius = boundary.pRad * 4;
	const double density0 = IISPH_RestDensity; //�õ������ʼ�ܶ�
											   //�õ���̬���Ӹ���
	const unsigned int numDynamicParticles = cube.dynamicParticles.size();
	//������ǰ�ı߽�����
	for (int i = 0; i < (int)numDynamicParticles; i++)
	{
		double delta = CubicKernelOne::W_zero(radius);
		//�õ���Ӧ������bound���еľ�̬��������
		DynamicRigidParticle* a = cube.dynamicParticles[i];
		for (unsigned int j = 0; j < dynamicRigidNeighborList[i].size(); j++)
		{
			//�õ���Ӧ�ٽ����ӵ�ǰ���ٽ���
			DynamicRigidParticle* b = dynamicRigidNeighborList[i][j];
			//����õ���ǰ������������õ��Ĳ�ֵ
			Vector3f positionSub = a->m_x - b->m_x;
			delta += CubicKernelOne::W(radius, positionSub);
		}
		const double volume = 1.0 / delta;
		//�Ե�ǰ��̬���ӵ�boundaryPsi���и�ֵ
		cube.dynamicParticles[i]->m_boundaryPsi = density0 * volume;
		cout << cube.dynamicParticles[i]->m_boundaryPsi << endl;
	}
}

void IISPHFluidObject::Initialise()
{
	particleList.clear();
	boundaryNeighborList.clear();
	dynamicRigidNeighborList.clear();

	//�Ƿ�����������ģ��
	boundary.staticRigidFlag = staticPilarFlag;
	ComputeFluidWLH();
	ComputeParticleNum(FluidWLH.GetX(), FluidWLH.GetY(), FluidWLH.GetZ()); 
	InitialiseParticlePosition(); 
	InitialiseParticleDensity();
	InitializeBoundary();

	//���ݳ�ʼ�������Ӱ뾶ȷ����ǰ��֧�Ű뾶
	Poly6Kernel::setRadius(4 * particleInitialPad);
	SpikyKernel::setRadius(4 * particleInitialPad);
	CubicKernel::setRadius(4 * particleInitialPad);
}

void IISPHFluidObject::ComputeFluidWLH()
{
	int numX = floor(double(IISPH_FluidWitdth / (2 * particleInitialPad)));
	int numY = floor(double(IISPH_FluidLength / (2 * particleInitialPad)));
	int numZ = floor(double(IISPH_FluidHeight / (2 * particleInitialPad)));
 	FluidWLH.SetValue(numX, numY, numZ);
	if (doubleDamBreak == true)
	{
		int numX2 = floor(double(IISPH_FluidWitdth2 / (2 * particleInitialPad)));
		int numY2 = floor(double(IISPH_FluidLength2 / (2 * particleInitialPad)));
		int numZ2 = floor(double(IISPH_FluidHeight2 / (2 * particleInitialPad)));
		FluidWLH2.SetValue(numX, numY, numZ);
	}
}

void IISPHFluidObject::InitialiseParticlePosition()
{
	//��ʼ�����ӵ�λ��
	//�����б��е����
	int particleIndex = 0;
	//���ӵ�ֱ��
	double particleLength = particleInitialPad * 2;
	for (int i = 0; i < FluidWLH.GetX(); i++)
	{
		for (int j = 0; j < FluidWLH.GetY(); j++)
		{
			for (int k = 0; k < FluidWLH.GetZ(); k++)
			{
				Vector3f position(0, 0, 0); //��ʼ����ǰ���������ӵ�λ��
				position.SetX(IISPH_OriginX + particleInitialPad + particleLength*i); //���õ�ǰ���ӵ�x����
				position.SetY(IISPH_OriginY + particleInitialPad + particleLength*j); //���õ�ǰ���ӵ�y����
				position.SetZ(IISPH_OriginZ + particleInitialPad + particleLength*k); //���õ�ǰ���ӵ�z����
				int hashValue = FunctionKit::PositionMapHash(position,
					SharedData::GetCellLength(),
					SharedData::GetHashCellNum(),
					boundary.offset); //�õ���ǰposition����Ӧ�Ĺ�ϣֵ
				IISPHParticle* newParticle=new IISPHParticle(); //��ʼ��һ���µ�����
				newParticle->Initialization(position, hashValue, particleInitialPad); //�������ӽ��г�ʼ��
                //����ʼ���õ�������װ�������б���
				particleList[particleIndex] = newParticle;
				particleIndex++;
			}
		}
	}
	if (doubleDamBreak == true)
	{
		for (int i = 0; i < FluidWLH2.GetX(); i++)
		{
			for (int j = 0; j < FluidWLH2.GetY(); j++)
			{
				for (int k = 0; k < FluidWLH2.GetZ(); k++)
				{
					Vector3f position(0, 0, 0); //��ʼ����ǰ���������ӵ�λ��
					position.SetX(IISPH_OriginX2 + particleInitialPad + particleLength*i); //���õ�ǰ���ӵ�x����
					position.SetY(IISPH_OriginY2 + particleInitialPad + particleLength*j); //���õ�ǰ���ӵ�y����
					position.SetZ(IISPH_OriginZ2 + particleInitialPad + particleLength*k); //���õ�ǰ���ӵ�z����
					int hashValue = FunctionKit::PositionMapHash(position,
						SharedData::GetCellLength(),
						SharedData::GetHashCellNum(),
						boundary.offset); //�õ���ǰposition����Ӧ�Ĺ�ϣֵ
					IISPHParticle* newParticle = new IISPHParticle(); //��ʼ��һ���µ�����
					newParticle->Initialization(position, hashValue, particleInitialPad); //�������ӽ��г�ʼ��
																						  //����ʼ���õ�������װ�������б���
					//newParticle->particleMass *= 0.5f;
					particleList[particleIndex] = newParticle;
					particleIndex++;
				}
			}
		}
	}
}

void IISPHFluidObject::InitialiseParticleDensity()
{
	//��ˮ���������ӵ��ܶ�����Ϊˮ��ĳ�ʼ�ܶ�
	for (int i = 0; i < particleList.size(); i++)
	{
		particleList[i]->density = IISPH_RestDensity;
		//����������
		//particleList[i]->lastDensity = IISPH_RestDensity;
	}
}

void IISPHFluidObject::UpdateParticleHashValue()
{
	//ˢ�µ�ǰ���ӵĹ�ϣֵ
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

void IISPHFluidObject::ComputeParticleNum(int x, int y, int z)
{
	IISPH_ParticleNum = x*y*z;               //����x,y,z�᲻ͬ�ķ�����������������
	if (doubleDamBreak == true)
	{
		//���˫����ģ��Ϊ�棬��ôˢ�µ�ǰIISPH����������
		IISPH_ParticleNum += FluidWLH2.GetX()*FluidWLH2.GetY()*FluidWLH2.GetZ();
	}
	particleList.resize(IISPH_ParticleNum);  //����������������ˢ�����������ھ��б�����
	//neighborList.resize(IISPH_ParticleNum);
	//fluidBoundaryNeighborList.resize(IISPH_ParticleNum);
}
