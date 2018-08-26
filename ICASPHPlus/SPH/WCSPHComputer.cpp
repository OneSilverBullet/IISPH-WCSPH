#include "WCSPHComputer.h"
//目前正在修改整理的版本2018——3——30

//初始化DEBUG 完成！！！
void WCSPHComputation::Initialization()
{
	//对哈希网格以及流体模型都使用默认参数进行初始化
	hashGridList.boundary = boundary;
	hashGridList.Intialization(); 
	fluidModel.boundary = boundary;
	fluidModel.Initialise();
	//将所有的粒子都装入对应的网格当中
	MapParticleToGrid();
	//将所有的边界粒子装入对应的网格当中
	MapBoundaryParticleToGrid();
	//根据网格当中的邻居边界粒子，计算得到边界粒子的边界粒子邻居列表
	ComputeBoundaryNeighborList();
	//计算边界粒子的PSI值
	fluidModel.ComputeBoundaryPsi(fluidModel.boundaryObj);
	//需要再一次将边界粒子装入网格当中，来更新BoundaryPSI值
	MapBoundaryParticleToGrid();
	//在所有的流体粒子以及边界粒子完成网格装填之后
	//更新当前流体粒子的流体粒子邻居表
	UpdateFluidNeighborList();
	//更新当前流体粒自的边界粒子邻居表
	UpdateFluidBoundaryNeighborList();

	//对WCSPHComputation当中的参数装入全局变量当中
	SharedData::SetStiffnesss(stiffnesss);
	SharedData::SetGamma(gamma);
	SharedData::SetGravity(gravity);
	SharedData::SetViscosityConstant(viscosityConstant);
}

//将粒子装入到对应的哈希表容器当中
//注意：调用这个函数之前，需要把所有的粒子都进行初始化
void WCSPHComputation::MapParticleToGrid()
{
	//清除网格当中的粒子之后再向其中添加新的粒子
	//将流体粒子装入对应位置
	hashGridList.ClearParticle();
	for (int i = 0; i < fluidModel.particleList.size(); i++)
	{
		hashGridList.PushParticle(fluidModel.particleList[i]);
	}
}

//边界粒子的网格映射只需要执行一次就好，不需要每次随着流体的变化而变化
void WCSPHComputation::MapBoundaryParticleToGrid()
{
	for (int i = 0; i < fluidModel.boundaryObj.GetStaticParticleNum(); i++)
	{
		StaticRigidParticle* tempStaticParticle = fluidModel.boundaryObj.staticRigidParticleList[i];
		hashGridList.PushParticle(tempStaticParticle);
	}
	//把所有的边界粒子存入边界网格当中
}

//锁帧，在实际运行当中只需要调用这个函数就可以了
//这里先不这么做
void WCSPHComputation::Frame()
{
	//帧数
	int framePreSecond = FRAME_PER_SECOND;
	//进行各种类型转换，得到每一帧程序需要间隔的时间点数，从而实现锁帧
	int CalculateDotNum = (int)(((float)1.0 / framePreSecond) / (float)(fluidModel.timeStep));
	int i = CalculateDotNum;
	for (;;)
	{
		Computation();
		if (i < 0)
		{
			break;
		}
		i--;
	}
}

//注意这个函数当中的计算顺序，不可变
void WCSPHComputation::Computation()
{
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int i = 0; i < fluidModel.particleList.size(); i++)
		{
			//得到流体粒子的流体粒子临近列表
			vector<WCSPHParticle*> neighborList = fluidModel.neighborList[i];
			//得到流体粒子的刚体粒子临近列表
			vector<StaticRigidParticle*> boundaryNeighborList = fluidModel.fluidBoundaryNeighborList[i];
			//依据临近列表得到当前粒子的密度
			double densityT = ComputeDensity(neighborList, boundaryNeighborList, fluidModel.particleList[i]);
			//cout << densityT << endl;
			fluidModel.particleList[i]->density = densityT;
			//依据临近列表得到当前粒子的粘性力
			Vector3f viscosityA = ComputeViscostyAcceleration(neighborList, boundaryNeighborList, fluidModel.particleList[i]);
			fluidModel.particleList[i]->viscosity = viscosityA;
			//对当前粒子密度进行校正，并且完成对粒子P的计算
			double correctDensity = CorrectionParticleDensity(fluidModel.particleList[i]);
			double  P = CalculateParticleTaite(correctDensity);
			fluidModel.particleList[i]->P = P;
		}
	}

    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static) 
		//经过上述计算，每个粒子的P可以确定，从而可以继续进行压力加速度、以及粒子的变化的计算
		for (int j = 0; j < fluidModel.particleList.size(); j++)
		{
			//得到临近列表
			vector<WCSPHParticle*> neighbors = fluidModel.neighborList[j];
			//得到流体粒子的刚体粒子临近列表
			vector<StaticRigidParticle*> boundaryNeighborList = fluidModel.fluidBoundaryNeighborList[j];
			//依据临近列表计算得到当前粒子的压力加速度
			Vector3f pressureA = ComputePressureAcceleration(neighbors, boundaryNeighborList, fluidModel.particleList[j]);
			fluidModel.particleList[j]->pressure = pressureA;
			//依据上面的计算，得到当前粒子所得到的总加速度
			Vector3f accelerationT = ComputeAcceleration(fluidModel.particleList[j]->pressure, fluidModel.particleList[j]->viscosity);
			fluidModel.particleList[j]->acceleration = accelerationT;
		}
	}
	
	//刷新流体模型当中粒子的加速度
	fluidModel.UpdateParticlePosition();
	//刷新所有粒子的哈希值，便于最终存入网格
	//运行笔记：如果这里不存入所有的粒子，那么所有粒子的邻近搜索在一定时间全部变为0
	fluidModel.UpdateParticleHashValue();
	//在所有粒子的位置和速度都计算完成之后，重新对粒子进行网格划分
	MapParticleToGrid();
	//刷新之后，更新流体的流体粒子邻居列表
	UpdateFluidNeighborList();
	//更新流体粒子的边界粒子邻居列表
	UpdateFluidBoundaryNeighborList();
	//完成一次计算
}

//这个函数返回当前粒子的周围粒子
vector<WCSPHParticle*> WCSPHComputation::ComputeNeighborParticle(WCSPHParticle* origin)
{
	double h = SharedData::GetCellLength(); //得到粒子的支撑半径，也就是网格单元的长度
	Vector3f originPosition = origin->position;
	//得到角落里的粒子表格
	Vector3f tempPosition; //用于对27个网格进行遍历的工具
	tempPosition.SetX(originPosition.GetX() - h);
	tempPosition.SetY(originPosition.GetY() - h);
	tempPosition.SetZ(originPosition.GetZ() - h);
	double originX = tempPosition.GetX();
	double originY = tempPosition.GetY();
	double originZ = tempPosition.GetZ();
	vector<WCSPHParticle*> neighborList; //粒子的邻居粒子列表
	int cellNum = hashGridList.GetCellNum(); //得到当前网格的数量
	for (int i = 0; i < 3; i++, tempPosition.SetX(tempPosition.GetX() + h))
	{
		for (int j = 0; j < 3; j++, tempPosition.SetY(tempPosition.GetY() + h))
		{
			for (int k = 0; k < 3; k++, tempPosition.SetZ(tempPosition.GetZ() + h))
			{
				//运行笔记：由于找不到对应的网格，返回的网格是错误的，所以，这个地方的代码是错误的
				//解决原因：其中一个函数调用参数用的是网格数量而不是哈希表格的数量
				GridCell targetGridCell = hashGridList.GetGridCell(tempPosition);
				//根据目标粒子的位置，直接查找到对应的哈希列表当中的网格，并且获取其中的粒子列表
				//运行笔记：找到对应网格之后，无法找到网格当中的粒子
				//解决方式：粒子映射网格拷贝错误，导致无法找到网格当中的粒子
				vector<WCSPHParticle*> tempParticleList = targetGridCell.GetParticleList();
				//将目标粒子列表当中的粒子一个一个的塞进邻居列表
				for (int q = 0; q < tempParticleList.size(); q++)
				{
					//如果在邻域搜索当中出现自己，那么就忽略掉
					if (origin->position == tempParticleList[q]->position)
					{
						continue;
					}
					//只压入在h范围之内（在粒子的影响范围之内）的邻居粒子
					if (FunctionKit::DistanceComputation(originPosition,
						tempParticleList[q]->position) < h)
					{
						neighborList.push_back(tempParticleList[q]);
					}
					//其他的舍弃
				}
			}
			tempPosition.SetZ(originZ);
		}
		tempPosition.SetY(originY);
	}
	return neighborList;
}

//查询流体粒子周围的边界粒子
vector<StaticRigidParticle*> WCSPHComputation::ComputeNeighborBoundaryParticle(WCSPHParticle* origin)
{
	double h = SharedData::GetCellLength(); //得到粒子的支撑半径，也就是网格单元的长度
	Vector3f originPosition = origin->position;
	//得到角落里的粒子表格
	Vector3f tempPosition; //用于对27个网格进行遍历的工具
	tempPosition.SetX(originPosition.GetX() - h);
	tempPosition.SetY(originPosition.GetY() - h);
	tempPosition.SetZ(originPosition.GetZ() - h);
	double originX = tempPosition.GetX();
	double originY = tempPosition.GetY();
	double originZ = tempPosition.GetZ();
	vector<StaticRigidParticle*> neighborList; //粒子的邻居粒子列表
	int cellNum = hashGridList.GetCellNum(); //得到当前网格的数量

	for (int i = 0; i < 3; i++, tempPosition.SetX(tempPosition.GetX() + h))
	{
		for (int j = 0; j < 3; j++, tempPosition.SetY(tempPosition.GetY() + h))
		{
			for (int k = 0; k < 3; k++, tempPosition.SetZ(tempPosition.GetZ() + h))
			{
				//运行笔记：由于找不到对应的网格，返回的网格是错误的，所以，这个地方的代码是错误的
				//解决原因：其中一个函数调用参数用的是网格数量而不是哈希表格的数量
				GridCell targetGridCell = hashGridList.GetGridCell(tempPosition);
				//根据目标粒子的位置，直接查找到对应的哈希列表当中的网格，并且获取其中的粒子列表
				//运行笔记：找到对应网格之后，无法找到网格当中的粒子
				//解决方式：粒子映射网格拷贝错误，导致无法找到网格当中的粒子
				vector<StaticRigidParticle*> tempParticleList = targetGridCell.GetBoundaryParticleList();
				//将目标粒子列表当中的粒子一个一个的塞进邻居列表
				for (int q = 0; q < tempParticleList.size(); q++)
				{
					//如果在邻域搜索当中出现自己，那么就忽略掉
					if (origin->position == tempParticleList[q]->position)
					{
						continue;
					}
					//只压入在h范围之内（在粒子的影响范围之内）的邻居粒子
					if (FunctionKit::DistanceComputation(originPosition,
						tempParticleList[q]->position) < h)
					{
						neighborList.push_back(tempParticleList[q]);
					}
					//其他的舍弃
				}
			}
			tempPosition.SetZ(originZ);
		}
		tempPosition.SetY(originY);
	}
	return neighborList;
}

vector<StaticRigidParticle*> WCSPHComputation::ComputeNeighborBoundaryParticle(StaticRigidParticle* origin)
{
	double h = SharedData::GetCellLength(); //得到粒子的支撑半径，也就是网格单元的长度
	Vector3f originPosition = origin->position;
	//得到角落里的粒子表格
	Vector3f tempPosition; //用于对27个网格进行遍历的工具
	tempPosition.SetX(originPosition.GetX() - h);
	tempPosition.SetY(originPosition.GetY() - h);
	tempPosition.SetZ(originPosition.GetZ() - h);
	double originX = tempPosition.GetX();
	double originY = tempPosition.GetY();
	double originZ = tempPosition.GetZ();
	vector<StaticRigidParticle*> neighborList; //粒子的邻居粒子列表
	int cellNum = hashGridList.GetCellNum(); //得到当前网格的数量

	for (int i = 0; i < 3; i++, tempPosition.SetX(tempPosition.GetX() + h))
	{
		for (int j = 0; j < 3; j++, tempPosition.SetY(tempPosition.GetY() + h))
		{
			for (int k = 0; k < 3; k++, tempPosition.SetZ(tempPosition.GetZ() + h))
			{
				//运行笔记：由于找不到对应的网格，返回的网格是错误的，所以，这个地方的代码是错误的
				//解决原因：其中一个函数调用参数用的是网格数量而不是哈希表格的数量
				GridCell targetGridCell = hashGridList.GetGridCell(tempPosition);
				//根据目标粒子的位置，直接查找到对应的哈希列表当中的网格，并且获取其中的粒子列表
				//运行笔记：找到对应网格之后，无法找到网格当中的粒子
				//解决方式：粒子映射网格拷贝错误，导致无法找到网格当中的粒子
				vector<StaticRigidParticle*> tempParticleList = targetGridCell.GetBoundaryParticleList();
				//将目标粒子列表当中的粒子一个一个的塞进邻居列表
				for (int q = 0; q < tempParticleList.size(); q++)
				{
					//如果在邻域搜索当中出现自己，那么就忽略掉
					if (origin->position == tempParticleList[q]->position)
					{
						continue;
					}
					//只压入在h范围之内（在粒子的影响范围之内）的邻居粒子
					if (FunctionKit::DistanceComputation(originPosition,
						tempParticleList[q]->position) < h)
					{
						neighborList.push_back(tempParticleList[q]);
					}
					//其他的舍弃
				}
			}
			tempPosition.SetZ(originZ);
		}
		tempPosition.SetY(originY);
	}
	return neighborList;
}

double WCSPHComputation::ComputeDensity(vector<WCSPHParticle*>& neighborList,vector<StaticRigidParticle*>& boundaryNeighborList, WCSPHParticle* origin)
{
	int i;
	double mass = fluidModel.particleMass;
	double result = CubicKernel::W_zero()*mass;
	double h = hashGridList.GetCellLength(); //得到当前单元的长度
	Vector3f xi = origin->position;
	//计算流体粒子对流体粒子之间的影响

    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (i = 0; i < neighborList.size(); i++)
		{
			Vector3f xj = neighborList[i]->position;
			double distance = FunctionKit::DistanceComputation(origin->position, neighborList[i]->position);
			//在核函数的范围之内的话
			Vector3f temp123 = xi - xj;
			double temp = CubicKernel::W(xi - xj);
			result += CubicKernel::W(xi - xj)*mass;
			//result += SPHKernal::Poly6Kernel(h, distance);
		}
	}
	//计算边界粒子对流体粒子的影响，计算密度补正

    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int j = 0; j < boundaryNeighborList.size(); j++)
		{
			Vector3f xj = boundaryNeighborList[j]->position;
			double kernelResult = CubicKernel::W(xi - xj);
			double boundPSI = boundaryNeighborList[j]->boundaryPsi;
			//调用得到对应的边界邻居粒子
			result += boundaryNeighborList[j]->boundaryPsi * CubicKernel::W(xi - xj);
		}
	}
	return result;
}

//注意：在该函数调用前，必须计算完所有粒子的密度以及压强
Vector3f WCSPHComputation::ComputePressureAcceleration(vector<WCSPHParticle*>& neighborList, vector<StaticRigidParticle*>& boundaryNeighborList, WCSPHParticle* origin)
{
	//注意 该函数使用的基本的SPH的密度计算方式S
	int j;
	double m = fluidModel.particleMass;
	double h = fluidModel.particleSupportRad; //得到支撑半径
	double pi = origin->P; //得到压强
	double di = origin->density; //得到当前粒子的密度
	double mass = fluidModel.particleMass; //得到粒子的质量
	Vector3f xi = origin->position; //得到当前粒子的位置
	//初始计算
	double dpi = pi / (di*di);
	Vector3f result(0, 0, 0);

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static) 
		for (j = 0; j < neighborList.size(); j++)
		{
			double pj = neighborList[j]->P;
			double dj = neighborList[j]->density;
			double r = FunctionKit::DistanceComputation(neighborList[j]->position, origin->position);
			Vector3f xj = neighborList[j]->position;
			//这里没有使用核函数，使用的是基本的SPH方法
			//运行笔记：由该核函数求值为负数
			//解决方法：省略掉对自身的计算
			/*
			Vector3f vecA = rj - ri;
			double additionP = pi + pj;
			double multi = 2 * di*dj;
			double disSub = h - distance;
			Vector3f kl = vecA * (additionP / multi *disSub*disSub / distance);
			result = result + kl;
			*/
			double dpj = pj / (dj*dj);
			result -= mass * (dpi + dpj) *  CubicKernel::gradW(xi - xj);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
    //根据边界粒子计算粒子应该受到的力
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static) 
		for (int j = 0; j < boundaryNeighborList.size(); j++)
		{
			StaticRigidParticle* boundaryNeighborParticle = boundaryNeighborList[j];
			Vector3f xj = boundaryNeighborParticle->position;
			Vector3f a = boundaryNeighborParticle->boundaryPsi * (dpi)*CubicKernel::gradW(xi - xj);
			result -= a;
		}
	}
	return result;
}

//计算粘性力，这里使用的是基础SPH的计算方式
//注意在该函数调用前，必须计算完所有粒子的速度
Vector3f WCSPHComputation::ComputeViscostyAcceleration(vector<WCSPHParticle*>& neighborList, vector<StaticRigidParticle*>& boundaryNeighborList, WCSPHParticle* origin)
{
	//此处调用的是XSPH当中的Viscosity的计算方法
	Vector3f vi = origin->velocity;
	Vector3f xi = origin->position;
	double m = fluidModel.particleMass;
	double viscosityCons = viscosityConstant;
	double h = fluidModel.particleSupportRad;
	double di = origin->density;
	Vector3f result(0, 0, 0);
	double invH = 1.0 / SharedData::GetTimeStep();
	//////////////////////////////////////////////////////////////////////////
	///流体自身的黏度计算
	//////////////////////////////////////////////////////////////////////////
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		for (int j = 0; j < neighborList.size(); j++)
		{
			Vector3f vj = neighborList[j]->velocity;
			Vector3f xj = neighborList[j]->position;
			double dj = neighborList[j]->density;
			double r = FunctionKit::DistanceComputation(xi, xj); //得到两个点的长度
			result -= invH*viscosityCons*(m / dj)*(vi - vj)*CubicKernel::W(xi - xj);
		}
	}
	//处理边界对流体黏度的影响
	//////////////////////////////////////////////////////////////////////////
	///边界对流体黏度的影响
	//////////////////////////////////////////////////////////////////////////
	//增添下面的代码，模拟出类似于火山岩浆的效果，流体黏度很大
  /*  #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		for (unsigned int j = 0; j < boundaryNeighborList.size(); j++)
		{
			const Vector3f xj = boundaryNeighborList[j].position;
			const Vector3f vj = boundaryNeighborList[j].velocity;
			result -= invH * viscosityCons * (boundaryNeighborList[j].boundaryPsi / di) * (vi)* CubicKernel::W(xi - xj);
		}
	}*/
	return result;
}

//计算合力加速度，从而完成全部计算
Vector3f WCSPHComputation::ComputeAcceleration(Vector3f pressureA, Vector3f viscosityA)
{
	//重力的加速度,注意其Z轴的标记，以及这个值是向下的
	Vector3f gA(0, 0, gravity); 
	Vector3f result = pressureA + viscosityA + gA;
	return result;
}

double WCSPHComputation::CorrectionParticleDensity(WCSPHParticle* a)
{
	return max(a->density, fluidModel.restDensity); //矫正当前粒子密度，用于计算压强
}

//调用泰特方程计算得到粒子当前的压强
double WCSPHComputation::CalculateParticleTaite(double a)
{
	return stiffnesss * (pow(a / fluidModel.restDensity, gamma) - 1.0);
}

void WCSPHComputation::UpdateFluidNeighborList()
{
	//在刷新新的粒子列表的时候，清除当前的邻居列表
	fluidModel.ClearNeighborList();
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		//遍历当前所有粒子的列表
		for (int i = 0; i < fluidModel.particleList.size(); i++)
		{
			//计算得到所有的临近的粒子
			vector<WCSPHParticle*> temp = ComputeNeighborParticle(fluidModel.particleList[i]);
			for (int j = 0; j < temp.size(); j++)
			{
				//将临近的所有粒子装入对应的粒子邻居表当中
				fluidModel.neighborList[i].push_back(temp[j]);
			}
		}
	}
}

//这个函数需要粒子每一次都需要计算的
void WCSPHComputation::UpdateFluidBoundaryNeighborList()
{
	//清理流体周围的边界粒子
	fluidModel.ClearFluidBoundaryNeighborList();
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		for (int i = 0; i < fluidModel.particleList.size(); i++)
		{
			//得到对应的流体粒子中的临近静态粒子列表
			vector<StaticRigidParticle*> temp = ComputeNeighborBoundaryParticle(fluidModel.particleList[i]);
			for (int j = 0; j < temp.size(); j++)
			{
				//将对应的静态流体粒子装入到流体周围的静态粒子列表当中
				fluidModel.fluidBoundaryNeighborList[i].push_back(temp[j]);
			}
		}
	}
}

//这个函数只需要在初始BoundaryPSI初始化的时候进行
void WCSPHComputation::ComputeBoundaryNeighborList()
{
	for (int i = 0; i < fluidModel.boundaryObj.staticRigidParticleList.size(); i++)
	{
		StaticRigidParticle* test2 = fluidModel.boundaryObj.staticRigidParticleList[i];
		//得到对应的流体粒子中的临近静态粒子列表
		vector<StaticRigidParticle*> temp = ComputeNeighborBoundaryParticle(fluidModel.boundaryObj.staticRigidParticleList[i]);
		for (int j = 0; j < temp.size(); j++)
		{
			//将对应的静态流体粒子装入到静态粒子周围的静态粒子列表当中
			fluidModel.boundaryNeighborList[i].push_back(temp[j]);
		}
	}	
}
