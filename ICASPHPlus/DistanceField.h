#pragma once
#include <vector>
#include "HGrid.h"
#include "HashGrid.h"
#include "Point3D.h"
#include "../ICASPHPlus/SPH/IISPHParticle.h"
//#include "Particle.h"
//using namespace std;

/***************
*哈希链表，网格节点类.
*每个对象代表一个网格节点。
*/
class DistanceNode
{
public:
	DistanceNode(double ox, double oy, double oz, GridIndex &gi)
	{
		this->x = ox;
		this->y = oy;
		this->z = oz;
		this->data = MAX_DIS;
		gridIndex = gi;
	}

	DistanceNode(double ox, double oy, double oz)
	{
		this->x = ox;
		this->y = oy;
		this->z = oz;
	}

	DistanceNode(DistanceNode *gc)
	{
		if (gc == NULL)
			return;

		this->x = gc->x;
		this->y = gc->y;
		this->z = gc->z;
		this->gridIndex = gc->gridIndex;
		this->next = gc->next;
		this->data = gc->data;
	}

	bool operator=(DistanceNode *gc)
	{
		if (gc == NULL)
			return false;

		this->x = gc->x;
		this->y = gc->y;
		this->z = gc->z;
		this->gridIndex = gc->gridIndex;
		this->next = gc->next;
		this->data = gc->data;
		return true;
	}

	bool operator!=(DistanceNode &gc)
	{
		const double delt = 0.00000001;
		if (abs(this->x - gc.x) < delt && abs(this->y - gc.y) < delt && abs(this->z - gc.z) < delt)
		{
			return false;
		}

		return true;
	}

public:
	double data;	// 节点距离值
	double x, y, z;			// 坐标
	GridIndex gridIndex;		// 节点索引
	DistanceNode *next = NULL;
};

/***************
*哈希网格类
*/
class DistanceField : public HGrid
{
public:
	DistanceField();
	DistanceField(double len, uint hashSize, double xSize, double ySize, double zSize, double xOrigin = 0.0, double yOrigin = 0.0, double zOrigin = 0.0);
	~DistanceField();

public:
	int numNode = 0;	//节点总数

public:
	vector<DistanceNode*> hash;//节点哈希表
	vector<Point3D> idxGridCatalog;//网格索引目录

	void partIdxGrid();

public:
	//根据网格单元长度，和空间xyz方向长度、起点，将距离场划分为空间网格，确定各个节点
	bool Partition(double &len, double xSize, double ySize, double zSize, double ox, double oy, double oz);
	//删除哈希表中每个网格节点
	bool DeleteHash();
	//清空网格空间距离场
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
	//通过哈希表获取网格节点
	DistanceNode* GetGridNode(double x, double y, double z);
	//获取节点数据，距离值
	double GetNodeData(double x, double y, double z);

public:
	//将网格节点添加到哈希表
	bool AddGridNode(GridIndex idxGrid, DistanceNode* gridNode);
	//将节点距离值添加到网格节点
	bool AddData(double x, double y, double z, double data);
	//获取坐标距离值
	double GetDisValue(double x, double y, double z, const vector<IISPHParticle*> &particles,  double radius = 1.0);
};

