#include "WCSPHFluidObject.h"

void WCSPHFluidObject::Initialise()
{
	SetParticleRad();  //更新当前光滑核的相关信息
	ComputeFluidWLH(); //根据默认参数计算得到WLH
	ComputeParticleNum(fluidWLH.GetX(), fluidWLH.GetY(), fluidWLH.GetZ()); //依据计算得到的流体长宽高计算得到粒子数量并且进行扩容
	InitialiseParticlePosition(); //依据以上信息进行粒子位置初始化
	InitialiseParticleDensity();
	InitialiseParticleP();
	neighborList.resize(particleNum); //对邻居表进行扩展
	fluidBoundaryNeighborList.resize(particleNum); //对流体粒子的边界粒子邻居进行扩展
												   //对边界粒子进行初始化
	InitializeBoundary();
	//对全局变量进行初始化
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
	//数据初始化
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
	//依据基本信息进行再次计算
	ComputeFluidWLH(); //根据默认参数计算得到WLH
	ComputeParticleNum(fluidWLH.GetX(), fluidWLH.GetY(), fluidWLH.GetZ()); //依据计算得到的流体长宽高计算得到粒子数量并且进行扩容
	InitialiseParticlePosition(); //依据以上信息进行粒子位置初始化
	InitialiseParticleDensity();
	InitialiseParticleP();
	neighborList.resize(particleNum); //初始化水体的邻居表
									  //对边界粒子进行初始化
	InitializeBoundary();
	//下面对SharedData进行初始化设置
	SharedData::SetFluidOrigin(fluidOrigin);
	SharedData::SetBoundaryWLH(WLH);
	SharedData::SetRestDensity(rD);
	SharedData::SetParticleRad(pR);
	SharedData::SetParticleSupportRad(pR * 4);
	SharedData::SetTimeStep(time);
}

//这个函数需要在所有参数设置完毕之后调用
void WCSPHFluidObject::ComputeFluidWLH()
{
	int w = floor(double(fluidWitdth / (2 * particleRad))); //得到x方向上的长度
	int l = floor(double(fluidLength / (2 * particleRad))); //得到y方向上的长度
	int h = floor(double(fluidHeight / (2 * particleRad))); //得到z方向上的长度
	fluidWLH.SetValue(w, l, h);                             //对水体的长宽高方向进行初始化
}

void WCSPHFluidObject::ComputeFluidWLH(double fluidW, double fluidL, double fluidH, double r)
{
	int w = floor(double(fluidW / (2 * r))); //得到x方向上的长度
	int l = floor(double(fluidL / (2 * r))); //得到y方向上的长度
	int h = floor(double(fluidH / (2 * r))); //得到z方向上的长度
	fluidWLH.SetValue(w, l, h);
}

void WCSPHFluidObject::InitialiseParticlePosition()
{
	//初始化粒子的位置
	//粒子列表当中的序号
	int particleIndex = 0;
	//粒子的直径
	double particleLength = particleRad * 2;
	for (int i = 0; i < fluidWLH.GetX(); i++)
	{
		for (int j = 0; j < fluidWLH.GetY(); j++)
		{
			for (int k = 0; k < fluidWLH.GetZ(); k++)
			{
				Vector3f position(0, 0, 0); //初始化当前的最新粒子的位置
				position.SetX(originX + particleRad + particleLength*i); //设置当前粒子的x坐标
				position.SetY(originY + particleRad + particleLength*j); //设置当前粒子的y坐标
				position.SetZ(originZ + particleRad + particleLength*k); //设置当前粒子的z坐标
				int hashValue = FunctionKit::PositionMapHash(position,
					SharedData::GetCellLength(),
					SharedData::GetHashCellNum(),
					boundary.offset); //得到当前position所对应的哈希值
				WCSPHParticle* newParticle = new WCSPHParticle(); //初始化一个新的粒子
				newParticle->Initialization(position, hashValue); //对新粒子进行初始化
																  //将初始化好的新粒子装入粒子列表当中
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
			//刷新每个粒子的位置和速度
			Vector3f acc = particleList[i]->acceleration;
			// v = v + a * t 
			particleList[i]->velocity = particleList[i]->velocity + acc*timeStep;
			//p = p + v * t
			particleList[i]->position = particleList[i]->position + particleList[i]->velocity*timeStep;
			//对粒子的速度和位置进行纠正				
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
			double cellLength = SharedData::GetCellLength(); //得到网格的大小
			int hashCellNum = SharedData::GetHashCellNum();  //得到所有哈希网格数量
			int hashValue = FunctionKit::PositionMapHash(particleList[i]->position, cellLength, hashCellNum, boundary.offset);
			particleList[i]->hashIndex = hashValue;
		}
	}
}

//依据boundary的范围进行粒子的碰撞检测
//由于粒子作为刚体边界，会导致单个粒子溢出导致内存溢出，为了避免这种情况
//使用边界检测，一旦粒子溢出，那么就将粒子的速度强行设置为0
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
	//如果当前粒子超过相应的范围，那么就对其速度归0
	if (Flag)
	{
		a->velocity.Zero();
	}
}

void WCSPHFluidObject::ComputeParticleNum(int x, int y, int z)
{
	particleNum = x*y*z;               //根据x,y,z轴不同的方向来计算粒子数量
	particleList.resize(particleNum);  //并且依据粒子数量刷新粒子列表容量
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

//清除邻居表当中所有的信息
void WCSPHFluidObject::ClearNeighborList()
{
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static) 
		//注意仅仅清除当前邻居列表里的列表
		for (int i = 0; i < neighborList.size(); i++)
		{
			neighborList[i].clear();
		}
	}

}

//清理流体边界的邻居列表
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

//这个函数的执行，必须在：
//1 流体模型的boundary已经从SPH Computation赋值完成
//2 流体模型的半径赋值给Kernel函数
void WCSPHFluidObject::InitializeBoundary()
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

//依据位置数组数据，对每一个粒子进行位置的赋值
void WCSPHFluidObject::AddRigidBodyObject(StaticRigidObject& rbo, const unsigned int numBoundaryParticles, vector<Vector3f> boundaryParticles)
{
	//从哈希网格列表当中得到两个关键的静态值,调用的都是静态值
	int cellNum = SharedData::GetHashCellNum();
	double cellLength = SharedData::GetCellLength();

	std::cout << "添加刚体对象" << std::endl;
	//RigidBodyParticleObject *rb = new RigidBodyParticleObject(); //新创建相应的刚体粒子对象
	//boundaryObj = rb;                                            //刚体粒子对象拷贝给本类的边界对象

	boundaryObj.staticRigidParticleList.resize(numBoundaryParticles);    //初始化当前静态刚体粒子列表的大小													
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
	std::cout << "刚体对象添加完成" << std::endl;
}


//该函数一定要在
void WCSPHFluidObject::ComputeBoundaryPsi(RigidBodyParticleObject& bound)
{
	const double density0 = restDensity; //得到流体初始密度
										 //得到静态粒子个数
	const unsigned int numBoundaryParticles = boundaryObj.GetStaticParticleNum();
	//得到边界粒子的支撑半径
	double supportRad = boundary.pRad * 4;
	//遍历当前的边界粒子
	for (int i = 0; i < (int)numBoundaryParticles; i++)
	{
		double delta = CubicKernelOne::W_zero(supportRad);
		//得到对应粒子在bound当中的静态刚性粒子
		StaticRigidParticle* a = boundaryObj.staticRigidParticleList[i];
		for (unsigned int j = 0; j < boundaryNeighborList[i].size(); j++)
		{
			//得到对应临近粒子当前的临近点
			StaticRigidParticle* b = boundaryNeighborList[i][j];
			//计算得到当前两个粒子相减得到的差值
			Vector3f positionSub = a->position - b->position;
			delta += CubicKernelOne::W(supportRad, positionSub);
		}
		const double volume = 1.0 / delta;
		//对当前静态粒子的boundaryPsi进行赋值
		boundaryObj.staticRigidParticleList[i]->boundaryPsi = density0 * volume;
		//std::cout << "  BoundaryPSI : " << boundaryObj.staticRigidParticleList[i]->boundaryPsi << std::endl;
	}
}
