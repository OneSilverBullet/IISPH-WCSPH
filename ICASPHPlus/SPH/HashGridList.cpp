#include "HashGridList.h"

//对整个GridCellList使用默认的参数
void HashGridList::Intialization()
{
	ClearBoundaryParticle();
	ComputeGridWLH(cellLength, boundary);    //计算得到WLH方向的网格数量
	ComputeCellNum();    //再计算cell的数量
	InitializeHashCellNun(); //初始化当前哈希网格列表的的网格数量
	EnsureTheGrid();     //在当前的基础上确定完成整个网格的设置
	//把网格当中的数据存入到对应的的全局参数当中
	SharedData::SetCellLength(cellLength);
	SharedData::SetHashCellNum(hashCellNum);
}

void HashGridList::Intialization(double cellL, double thred, Boundary bound)
{
	//初始化参数
	cellLength = cellL;
	boundary = bound;
	hashThred = thred;
	ComputeGridWLH(cellLength, boundary);    //计算得到WLH方向的网格数量
	ComputeCellNum();    //再计算cell的数量
	InitializeHashCellNun(); //初始化当前哈希网格列表的的网格数量
	EnsureTheGrid();     //在当前的基础上确定完成整个网格的设置
}

//依据边界以及新的cellLength计算得到在X,Y,Z方向的cell数量
void HashGridList::ComputeGridWLH(double cellL, Boundary bound)
{
	gridWLH.x = ceil(double(bound.gridBoundaryWidth / cellL));
	gridWLH.y = ceil (double(bound.gridBoundaryLength / cellL));
	gridWLH.z = ceil (double(bound.gridBoundaryHeight / cellL));
}

//经过该函数计算得到最终的空间网格数量
void HashGridList::ComputeCellNum()
{
	cellNum = gridWLH.x*gridWLH.y*gridWLH.z; //依据gridWLH的长宽高方向的cell数量相乘得到最终网格的数量
}

//哈希表的大小初始化为当前网格数量之后的第一个素数
void HashGridList::InitializeHashCellNun()
{
	int temp = cellNum;
	for (;;)
	{
		//如果当前值为素数，那么赋值，程序结束
		if (FunctionKit::CheckPrime(temp)) 
		{
			hashCellNum = temp;
			//运行注释：这里的动态容器第一次扩容为3200万，不能完成扩容
			//解决：尝试重新设计水体的长度以及空间的长度
			gridCellList.resize(hashCellNum);
			return;
		}
		temp++; //程序继续运行，直至找到一个素数
	}
}

//这个函数只需要在初始化的时候运行一遍
void HashGridList::EnsureTheGrid()
{
	GenerateGridCell(); //首先根据当前的网格对空间进行划分
	//程序无限循环，直至找到一个合适的哈希表
	for (;;)
	{
		//检测当前哈希表是否符合要求
		if (CheckThred()) 
		{
			return; //如果符合要求你那么就结束程序
		}
		else //当前哈希表不符合要求
		{
			Rehash(); 
			GenerateGridCell(); //rehash之后继续进行网格的生成
		}
	}
	return;
}

void HashGridList::GenerateGridCell()
{
	//根据上述进行初始化
	for (int i = 0; i < gridWLH.x; i++)
	{
		for (int j = 0; j < gridWLH.y; j++)
		{
			for (int k = 0; k < gridWLH.z; k++)
			{
				//取得当前位置的哈希表值
				int hashValue = FunctionKit::GridHash(i, j, k, hashCellNum);
				GridCell grid(i, j, k, hashValue);
				gridCellList[hashValue].push_back(grid);
			}
		}
	}
	
}

//这个函数的目的是，取得原来的hash值的二倍之后的第一个素数
void HashGridList::Rehash()
{
	int hashValue = hashCellNum * 2; 
	for (;;)
	{
		if (FunctionKit::CheckPrime(hashValue))
		{
			hashCellNum = hashValue; //如果当前hash值为素数，那么就更新
			gridCellList.resize(hashCellNum);
			return; //结束程序
		}
		hashValue++; //不断的增加hash值
	}
}

//用来检测此时的哈希表是否需要rehash
bool HashGridList::CheckThred()
{
	int storeNum = 0; //哈希表当中已经存储的数量
	for (int i = 0; i < gridCellList.size(); i++)
	{
		if (gridCellList[i].size() != 0)
		{
			storeNum++; //一旦当前的空位不为空，那么就将storeNum增加
		}
	}
	//运行笔记：此处遇到了中断的错误，原因是表的长度为0,从而导致除以0的错误
	//解决方式：通过解决哈希函数等问题解决了表长度显示不对的问题
	float ratio = float(float(storeNum) / gridCellList.size()); 
	if (ratio > hashThred) //装填因子超过阈值，那么重新设置位置
	{
		return false;  //返回false的时候，此时哈希表需要重新设置位置
	}
	return true; //返回true的时候，哈希表不用重新设置
}

void HashGridList::PushParticle(WCSPHParticle* a)
{
	int particleHash = a->hashIndex; //得到当前粒子的哈希值
	vector<GridCell> tempList = gridCellList[particleHash]; //得到当前粒子对应的列表
	//得到粒子所在的位置index
	Vector3i particleIndex = FunctionKit::PositionMapIndex(a->position, cellLength, boundary.offset);
	for (int i = 0; i < tempList.size(); i++)
	{
		//如果当前的粒子号和对应的网格匹配，那么确定网格，结束程序
		if (tempList[i].CheckGridCell(particleIndex))
		{
			//tempList[i].PushParticle(a);
			//这里深入拷贝错误！
			gridCellList[particleHash][i].PushParticle(a);
			return;
		}
	}
}

void HashGridList::PushParticle(StaticRigidParticle* a)
{
	int particleHash = a->hashValue; //得到当前粒子的哈希值
	vector<GridCell> tempList = gridCellList[particleHash]; //得到当前粒子对应的列表
   //得到粒子所在的位置index
	Vector3i particleIndex = FunctionKit::PositionMapIndex(a->position, cellLength, boundary.offset);
	for (int i = 0; i < tempList.size(); i++)
	{
		//如果当前的粒子号和对应的网格匹配，那么确定网格，结束程序
		if (tempList[i].CheckGridCell(particleIndex))
		{
			//tempList[i].PushParticle(a);
			//这里深入拷贝错误！
			gridCellList[particleHash][i].PushBoundaryParticle(a);
			//	cout << gridCellList[13651][0].GetBoundaryParticleList().size() << endl;
			return;
		}
	}
	return; //如果没有找到网格，需要有所设置！****待定****
}

void HashGridList::PushParticle(IISPHParticle* a)
{
	int particleHash = a->hashIndex; //得到当前粒子的哈希值
	vector<GridCell> tempList = gridCellList[particleHash]; //得到当前粒子对应的列表
															//得到粒子所在的位置index
	Vector3i particleIndex = FunctionKit::PositionMapIndex(a->position, cellLength, boundary.offset);
	for (int i = 0; i < tempList.size(); i++)
	{
		//如果当前的粒子号和对应的网格匹配，那么确定网格，结束程序
		if (tempList[i].CheckGridCell(particleIndex))
		{
			gridCellList[particleHash][i].PushParticle(a);
			return;
		}
	}
}

void HashGridList::PushParticle(DynamicRigidParticle * a)
{
	int particleHash = a->hashValue; //得到当前粒子的哈希值
	vector<GridCell> tempList = gridCellList[particleHash]; //得到当前粒子对应的列表
	//当前网格需要得到当前的动态刚体粒子的哈希值才可以
	for (int i = 0; i < tempList.size(); i++)
	{
		//如果当前的粒子号和对应的网格匹配，那么确定网格，结束程序
		if (tempList[i].CheckGridCell(a))
		{
			gridCellList[particleHash][i].PushDynamicParticle(a);
			return;
		}
	}
}

void HashGridList::ClearParticle()
{
	for (int i = 0; i < gridCellList.size(); i++)
	{
		if (gridCellList[i].size() != 0)
		{
			for (int j = 0; j < gridCellList[i].size(); j++)
			{
				gridCellList[i][j].ResetParticleList(); //清理当前的粒子列表
			}
		}
	}
}

void HashGridList::ClearBoundaryParticle()
{
	for (int i = 0; i < gridCellList.size(); i++)
	{
		if (gridCellList[i].size() != 0)
		{
			for (int j = 0; j < gridCellList[i].size(); j++)
			{
				gridCellList[i][j].ResetBoundaryParticleList(); //清理当前的粒子列表
			}
		}
	}
}

double HashGridList::GetCellLength()
{
	return cellLength;
}

int HashGridList::GetCellNum()
{
	return cellNum;
}

int HashGridList::GetHashCellNum()
{
	//得到当前哈希表的大小
	return hashCellNum;
}

vector<GridCell> HashGridList::GetGridCellVector(int hashIndex)
{
	return gridCellList[hashIndex];
}

vector<GridCell> HashGridList::GetGridCellVector(Vector3i positionIndex)
{
	int hashIndex = FunctionKit::GridHash(positionIndex.GetX(), positionIndex.GetY(), positionIndex.GetZ(), hashCellNum);
	return GetGridCellVector(hashIndex); //调用上一级重载函数
}

vector<GridCell> HashGridList::GetGridCellVector(Vector3f position)
{
	//根据位置映射到位置Index
	Vector3i positionIndex = FunctionKit::PositionMapIndex(position, cellLength, boundary.offset); 
	return GetGridCellVector(positionIndex); //调用上一级重载函数
}

GridCell HashGridList::GetGridCellWithSecVec(vector<GridCell>& gridVec, Vector3f parPosition)
{
	//运行笔记：这里键入对应值不返回对应的网格值***
	int i;
	for (i = 0; i < gridVec.size(); i++)
	{
		if (gridVec[i].CheckGridCell(parPosition)) //根据当前粒子定位对应的网格
		{
			return gridVec[i]; //返回对应的网格
		}
	}
	return GridCell();  //否则返回默认的网格
}

GridCell HashGridList::GetGridCellWithSecVec(vector<GridCell>& gridVec, WCSPHParticle a)
{
	return GetGridCellWithSecVec(gridVec, a.position); //调用上一级重载函数来完成网格的查找
}

GridCell HashGridList::GetGridCell(WCSPHParticle a)
{
	vector<GridCell> temp = GetGridCellVector(a.position); //得到当前粒子a所在哈希键位处的容器
	return GetGridCellWithSecVec(temp, a); //调用上面的函数完成依据粒子进行查询的过程
}

GridCell HashGridList::GetGridCell(Vector3f parPosition)
{
	vector<GridCell> temp = GetGridCellVector(parPosition); //根据当前粒子的位置找到对应的哈希键当中的容器
	return GetGridCellWithSecVec(temp, parPosition); //调用上面的函数完成查找
}


