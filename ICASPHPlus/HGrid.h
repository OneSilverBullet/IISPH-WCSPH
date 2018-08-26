#pragma once
#include <vector>

using namespace std;

/***************
*大质数和常量
*/
#define P1 73856093
#define P2 19349663
#define P3 83492791

#define uint unsigned

#define MAX_DIS 65536

/***************
*网格索引类，每个索引对应空间中的一个网格单元
*/
class GridIndex{
public:
	int x;
	int y;
	int z;

public:
	GridIndex(int x = 0, int y = 0, int z = 0)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	bool operator=(GridIndex &gi)
	{
		this->x = gi.x;
		this->y = gi.y;
		this->z = gi.z;

		return true;
	}

	bool operator==(GridIndex &gi)
	{
		if (this->x == gi.x &&this->y == gi.y&& this->z == gi.z)
		{
			return true;
		}
		return false;
	}

	bool operator!=(GridIndex &gi)
	{
		const double delt = 0.000000001;
		if (abs(this->x - gi.x) < delt && abs(this->y - gi.y) < delt && abs(this->z - gi.z) < delt)
		{
			return false;
		}

		return true;
	}
};

/***************
*哈希网格基类，抽象类，包含常用方法
*/
class HGrid
{
public:
	HGrid();
	~HGrid();

public:
	double lenght = 0.1;//节点间距（网格单元边长）
	double xLen;	//空间网格边长
	double yLen;
	double zLen;
	double ox;		//空间网格起点
	double oy;
	double oz;
	int numX = 0;		//网格单元体数量
	int numY = 0;
	int numZ = 0;

public:
	//根据网格单元长度，和空间xyz方向长度、起点，将距离场划分为空间网格，确定各个节点
	virtual bool Partition(double &len, double xSize, double ySize, double zSize, double ox, double oy, double oz) = 0;
	//删除哈希表中每个网格节点
	virtual bool DeleteHash() = 0;
	//清空网格空间中的元素
	virtual bool ClearGrid() = 0;

public:
	//获取网格索引坐标，输入坐标值
	virtual GridIndex GetGridIndex(double x, double y, double z) = 0;
	//获取哈希值，输入坐标位置
	virtual int GetHashIndex(double x, double y, double z) = 0;
	//获取哈希值，通过网格索引
	virtual int GetHashIndex(GridIndex &idxGrid) = 0;
	//获取哈希表长
	int GetHshSize();
};

