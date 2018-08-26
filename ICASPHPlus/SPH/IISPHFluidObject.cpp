#include "IISPHFluidObject.h"

void IISPHFluidObject::ClearNeighborList()
{
    #pragma omp parallel default(shared)
	{
       #pragma omp for schedule(static) 
		//注意仅仅清除当前邻居列表里的列表
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
	//依据边界粒子，计算对应的边界粒子的位置数组
	vector<Vector3f> boundaryParticlePosition = boundary.initBoundaryData();
	//依据计算得到的边界粒子的个数对边界邻居粒子列表进行扩容
	boundaryNeighborList.resize(boundaryParticlePosition.size());
	//初始化静态物体的属性
	StaticRigidObject rb;
	//使用静态粒子与边界粒子的位置对静态粒子边界进行初始化
	AddRigidBodyObject(rb, boundaryParticlePosition.size(), boundaryParticlePosition);

	//以下两步放在SPHComputation当中做
	//根据静态粒子的位置，将粒子映射到网格当中
	//计算每一个静态边界粒子的BoundaryPSI
}

void IISPHFluidObject::AddRigidBodyObject(StaticRigidObject & rbo, const unsigned int numBoundaryParticles, vector<Vector3f> boundaryParticles)
{
	//从哈希网格列表当中得到两个关键的静态值,调用的都是静态值
	int cellNum = SharedData::GetHashCellNum();
	double cellLength = SharedData::GetCellLength();

	std::cout << "添加刚体对象" << std::endl;
	boundaryObj.staticRigidParticleList.clear();
	boundaryObj.staticRigidParticleList.resize(numBoundaryParticles);    //初始化当前静态刚体粒子列表的大小
	//对其进行初始化
	//每一个指针都指向一个静态刚体粒子
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
		//根据当前位置计算得到对应粒子的哈希值
		boundaryObj.staticRigidParticleList[i]->hashValue =
			FunctionKit::PositionMapHash(boundaryParticles[i], cellLength,
				cellNum, boundary.offset);
	}
	//对边界对象的刚体属性进行赋值
	boundaryObj.rd = rbo;
}

void IISPHFluidObject::ComputeBoundaryPsi(RigidBodyParticleObject & bound)
{
	boundaryNumber = boundary.boundaryNumber;
	double radius = boundary.pRad * 4;
	const double density0 = IISPH_RestDensity; //得到流体初始密度
	//得到静态粒子个数
	const unsigned int numBoundaryParticles = boundaryObj.GetStaticParticleNum();
	//遍历当前的边界粒子
	for (int i = 0; i < (int)numBoundaryParticles; i++)
	{
		double delta = CubicKernelOne::W_zero(radius);
		//得到对应粒子在bound当中的静态刚性粒子
		StaticRigidParticle* a = boundaryObj.staticRigidParticleList[i];
		for (unsigned int j = 0; j < boundaryNeighborList[i].size(); j++)
		{
			//得到对应临近粒子当前的临近点
			StaticRigidParticle* b = boundaryNeighborList[i][j];
			//计算得到当前两个粒子相减得到的差值
			Vector3f positionSub = a->position - b->position;
			delta += CubicKernelOne::W(radius, positionSub);
		}
		const double volume = 1.0 / delta;
		//对当前静态粒子的boundaryPsi进行赋值
		boundaryObj.staticRigidParticleList[i]->boundaryPsi = density0 * volume;
		std::cout << "  BoundaryPSI : " << boundaryObj.staticRigidParticleList[i]->boundaryPsi << std::endl;
	}
}

//计算当前所有动态粒子的boundaryPsi
void IISPHFluidObject::ComputeBoundaryPsi(DynamicRigidParticleObject& cube)
{
	boundaryNumber = boundary.boundaryNumber;
	double radius = boundary.pRad * 4;
	const double density0 = IISPH_RestDensity; //得到流体初始密度
											   //得到静态粒子个数
	const unsigned int numDynamicParticles = cube.dynamicParticles.size();
	//遍历当前的边界粒子
	for (int i = 0; i < (int)numDynamicParticles; i++)
	{
		double delta = CubicKernelOne::W_zero(radius);
		//得到对应粒子在bound当中的静态刚性粒子
		DynamicRigidParticle* a = cube.dynamicParticles[i];
		for (unsigned int j = 0; j < dynamicRigidNeighborList[i].size(); j++)
		{
			//得到对应临近粒子当前的临近点
			DynamicRigidParticle* b = dynamicRigidNeighborList[i][j];
			//计算得到当前两个粒子相减得到的差值
			Vector3f positionSub = a->m_x - b->m_x;
			delta += CubicKernelOne::W(radius, positionSub);
		}
		const double volume = 1.0 / delta;
		//对当前静态粒子的boundaryPsi进行赋值
		cube.dynamicParticles[i]->m_boundaryPsi = density0 * volume;
		cout << cube.dynamicParticles[i]->m_boundaryPsi << endl;
	}
}

void IISPHFluidObject::Initialise()
{
	particleList.clear();
	boundaryNeighborList.clear();
	dynamicRigidNeighborList.clear();

	//是否进行流固耦合模拟
	boundary.staticRigidFlag = staticPilarFlag;
	ComputeFluidWLH();
	ComputeParticleNum(FluidWLH.GetX(), FluidWLH.GetY(), FluidWLH.GetZ()); 
	InitialiseParticlePosition(); 
	InitialiseParticleDensity();
	InitializeBoundary();

	//依据初始化的粒子半径确定当前的支撑半径
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
	//初始化粒子的位置
	//粒子列表当中的序号
	int particleIndex = 0;
	//粒子的直径
	double particleLength = particleInitialPad * 2;
	for (int i = 0; i < FluidWLH.GetX(); i++)
	{
		for (int j = 0; j < FluidWLH.GetY(); j++)
		{
			for (int k = 0; k < FluidWLH.GetZ(); k++)
			{
				Vector3f position(0, 0, 0); //初始化当前的最新粒子的位置
				position.SetX(IISPH_OriginX + particleInitialPad + particleLength*i); //设置当前粒子的x坐标
				position.SetY(IISPH_OriginY + particleInitialPad + particleLength*j); //设置当前粒子的y坐标
				position.SetZ(IISPH_OriginZ + particleInitialPad + particleLength*k); //设置当前粒子的z坐标
				int hashValue = FunctionKit::PositionMapHash(position,
					SharedData::GetCellLength(),
					SharedData::GetHashCellNum(),
					boundary.offset); //得到当前position所对应的哈希值
				IISPHParticle* newParticle=new IISPHParticle(); //初始化一个新的粒子
				newParticle->Initialization(position, hashValue, particleInitialPad); //对新粒子进行初始化
                //将初始化好的新粒子装入粒子列表当中
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
					Vector3f position(0, 0, 0); //初始化当前的最新粒子的位置
					position.SetX(IISPH_OriginX2 + particleInitialPad + particleLength*i); //设置当前粒子的x坐标
					position.SetY(IISPH_OriginY2 + particleInitialPad + particleLength*j); //设置当前粒子的y坐标
					position.SetZ(IISPH_OriginZ2 + particleInitialPad + particleLength*k); //设置当前粒子的z坐标
					int hashValue = FunctionKit::PositionMapHash(position,
						SharedData::GetCellLength(),
						SharedData::GetHashCellNum(),
						boundary.offset); //得到当前position所对应的哈希值
					IISPHParticle* newParticle = new IISPHParticle(); //初始化一个新的粒子
					newParticle->Initialization(position, hashValue, particleInitialPad); //对新粒子进行初始化
																						  //将初始化好的新粒子装入粒子列表当中
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
	//将水体所有粒子的密度设置为水体的初始密度
	for (int i = 0; i < particleList.size(); i++)
	{
		particleList[i]->density = IISPH_RestDensity;
		//！！！！！
		//particleList[i]->lastDensity = IISPH_RestDensity;
	}
}

void IISPHFluidObject::UpdateParticleHashValue()
{
	//刷新当前粒子的哈希值
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static) 
		for (int i = 0; i < particleList.size(); i++)
		{
			double cellLength = SharedData::GetCellLength(); //得到网格的大小
			int hashCellNum = SharedData::GetHashCellNum();  //得到所有哈希网格数量
			int hashValue = FunctionKit::PositionMapHash(particleList[i]->position, cellLength, hashCellNum, boundary.offset);
			particleList[i]->hashIndex = hashValue;
		}
	}
}

void IISPHFluidObject::ComputeParticleNum(int x, int y, int z)
{
	IISPH_ParticleNum = x*y*z;               //根据x,y,z轴不同的方向来计算粒子数量
	if (doubleDamBreak == true)
	{
		//如果双溃坝模型为真，那么刷新当前IISPH的粒子数量
		IISPH_ParticleNum += FluidWLH2.GetX()*FluidWLH2.GetY()*FluidWLH2.GetZ();
	}
	particleList.resize(IISPH_ParticleNum);  //并且依据粒子数量刷新流体粒子邻居列表容量
	//neighborList.resize(IISPH_ParticleNum);
	//fluidBoundaryNeighborList.resize(IISPH_ParticleNum);
}
