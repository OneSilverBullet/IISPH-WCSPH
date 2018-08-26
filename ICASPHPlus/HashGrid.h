#pragma once
#include "HGrid.h"


/***************
*哈希链表，网格单元类.
*每个对象代表一个网格单元。
*/
class GridCellDistance
{
public:
	GridCellDistance(double ox, double oy, double oz, GridIndex &gi)
	{
		this->ox = ox;
		this->oy = oy;
		this->oz = oz;

		gridIndex = gi;
	}

	GridCellDistance(GridCellDistance *gc)
	{
		if (gc == NULL)
			return;

		this->ox = gc->ox;
		this->oy = gc->oy;
		this->oz = gc->oz;
		//this->idxParticle = gc->idxParticle;
		this->gridIndex = gc->gridIndex;
		this->next = gc->next;
	}

	bool operator=(GridCellDistance *gc)
	{
		if (gc == NULL)
			return false;

		this->ox = gc->ox;
		this->oy = gc->oy;
		this->oz = gc->oz;
		//this->idxParticle = gc->idxParticle;
		this->gridIndex = gc->gridIndex;
		this->next = gc->next;

		return true;
	}
public:
	//vector<int> idxParticle;	// 网格单元内粒子索引，n个
	double ox, oy, oz;				// 网格单元起点，唯一
	GridIndex gridIndex;		// 网格单元的网格索引，唯一
	GridCellDistance *next = NULL;
};

/***************
*哈希网格类
*/
class HashGrid : public HGrid
{
public:
	HashGrid();
	HashGrid(double len, uint hashSize, double xSize, double ySize, double zSize, double xOrigin = 0.0, double yOrigin = 0.0, double zOrigin = 0.0);
	~HashGrid();

	bool operator=(HashGrid &hashGrid);
	int a, b, c;
	int d, e, f;
public:
	int numParticle = 0;//网格空间内粒子总数
	int numCell = 0;	//网格单元总数

public:
	vector<GridCellDistance*> hash;//网格哈希表

public:
	//根据网格单元长度，和空间xyz方向长度、起点，划分空间网格
	bool Partition(double &len, double xSize, double ySize, double zSize, double ox, double oy, double oz);
	//删除哈希表中每个网格单元
	bool DeleteHash();
	//清空网格空间的粒子
	bool ClearGrid();

public:
	//获取网格索引坐标，输入坐标值
	GridIndex GetGridIndex(double x, double y, double z);
	//获取哈希值，输入坐标位置
	int GetHashIndex(double x, double y, double z);
	//获取哈希值，通过网格索引
	int GetHashIndex(GridIndex &idxGrid);
	//获取哈希表长
	int GetHshSize();
	//通过哈希表获取网格单元
	GridCellDistance* GetGridCell(GridIndex &idxGrid);

public:
	//将网格单元添加到哈希表
	//bool AddHashCell(double x, double y, double z, GridCell* gridCell);
	//bool AddHashCell(int hashIndex, GridCell* gridCell);
	bool AddHashCell(GridIndex idxGrid, GridCellDistance* gridCell);
	//将粒子索引添加到网格单元
	//bool AddParticle(double x, double y, double z, int idx);
	//搜索点(x,y,z)周围半径radii范围内的邻点
	//bool SearchAdjParticles(double x, double y, double z, double radii, const vector<Particle> &particles, vector<int> &adjParticles);
};

