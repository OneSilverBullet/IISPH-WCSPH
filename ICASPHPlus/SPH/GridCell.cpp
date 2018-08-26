#include "GridCell.h"

void GridCell::Initialization(int i, int j, int k, int hashV)
{
	ComputeCellLength(); //计算当前网格长度
	ResetParticleList(); //重置粒子列表
	SetIndexVec(i, j, k); //设置该网格的位置信
	SetHashIndex(hashV); //进而设置其hashIndex
}

void GridCell::ComputeCellLength()
{
	//读取全局变量，设置当前网格
	gridCellLength = SharedData::GetCellLength();
}

void GridCell::PushParticle(WCSPHParticle* a)
{
	particleList.push_back(a);
}

void GridCell::PushParticle(IISPHParticle* a)
{
	IISPHparticleList.push_back(a);
}

void GridCell::PushBoundaryParticle(StaticRigidParticle* a)
{
	boundaryList.push_back(a); //将边界粒子装入对应的边界粒子李彪当中
}

void GridCell::PushDynamicParticle(DynamicRigidParticle * a)
{
	dynamicParticleList.push_back(a);
}

void GridCell::ResetParticleList()
{
	particleList.clear();
	IISPHparticleList.clear();
}

void GridCell::ResetBoundaryParticleList()
{
	boundaryList.clear();
}

void GridCell::ResetDynamicParticleList()
{
	dynamicParticleList.clear();
}

vector<WCSPHParticle*> GridCell::GetParticleList()
{
	return particleList;
}

vector<IISPHParticle*> GridCell::GetIISPHParticleList()
{
	return IISPHparticleList;
}

vector<StaticRigidParticle*> GridCell::GetBoundaryParticleList()
{
	return boundaryList;
}

vector<DynamicRigidParticle*> GridCell::GetDynamicParticleList()
{
	return dynamicParticleList;
}

bool GridCell::CheckGridCell(Vector3f parPosition)
{
	//访问全局变量，确定边界的偏移
	double temp = SharedData::GetOffset(); 
	//得到当前粒子的位置
	Vector3i indexVecPos = FunctionKit::PositionMapIndex(parPosition, gridCellLength, temp);
	//检测当前粒子所在的位置和当前网格的位置是一致的
	return CheckGridCell(indexVecPos); //调用上一级的重载函数
}

//通过位置标号来查询
bool GridCell::CheckGridCell(Vector3i parIndex)
{
	if (indexVec.GetX() == parIndex.GetX() &&
		indexVec.GetY() == parIndex.GetY() &&
		indexVec.GetZ() == parIndex.GetZ())
	{
		return true;
	}
	return false;
}

bool GridCell::CheckGridCell(WCSPHParticle particle)
{
	return CheckGridCell(particle.position); //调用上一级的重载函数
}

bool GridCell::CheckGridCell(IISPHParticle* particle)
{
	return CheckGridCell(particle->position); 
}

bool GridCell::CheckGridCell(DynamicRigidParticle * particle)
{
	return CheckGridCell(particle->m_x);
}

GridCell& GridCell::operator=(const GridCell & a)
{
	indexVec = a.indexVec;
	hashIndex = a.hashIndex;
	gridCellLength = a.gridCellLength;
	ResetParticleList();
	for (int i = 0; i < a.particleList.size(); i++)
	{
		PushParticle(a.particleList[i]);
	}
	//运行笔记：最基础的赋值操作符在扩展之后没有进行赋值边界粒子，导致无法查询到周围的粒子
	ResetBoundaryParticleList();
	for (int i = 0; i < a.boundaryList.size(); i++)
	{
		PushBoundaryParticle(a.boundaryList[i]);
	}
	ResetDynamicParticleList();
	for (int i=0; i<a.dynamicParticleList.size();i++)
	{
		PushDynamicParticle(a.dynamicParticleList[i]);
	}

	return *this;
}

void GridCell::SetGridCellLength(double length)
{
	gridCellLength = length;
}

double GridCell::GetGridCellLength()
{
	return gridCellLength;
}

void GridCell::SetIndexVec(Vector3i value)
{
	indexVec = value;
}

void GridCell::SetIndexVec(int i, int j, int k)
{
	Vector3i temp(i, j, k);
	indexVec = temp;
}

Vector3i GridCell::GetIndexVec()
{
	return indexVec;
}

void GridCell::SetHashIndex(int value)
{
	hashIndex = value;
}

int GridCell::GetHashIndex()
{
	return hashIndex;
}

