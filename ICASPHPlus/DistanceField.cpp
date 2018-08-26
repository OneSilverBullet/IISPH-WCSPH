#include "DistanceField.h"
#include <omp.h>

DistanceField::DistanceField()
{
}

DistanceField::DistanceField(double len, uint hashSize, double xSize, double ySize, double zSize, double xOrigin, double yOrigin, double zOrigin)
{
	if (xSize == 0)
	{
		xSize = 1;
	}
	if (ySize == 0)
	{
		ySize = 1;
	}
	if (zSize == 0)
	{
		zSize = 1;
	}

	this->lenght = len;
	this->xLen = xSize;
	this->yLen = ySize;
	this->zLen = zSize;
	this->ox = xOrigin;
	this->oy = yOrigin;
	this->oz = zOrigin;
	hash.resize(hashSize, NULL);

	//记录节点坐标
	partIdxGrid();
	//划分空间网格
	Partition(len, xSize, ySize, zSize, xOrigin, yOrigin, zOrigin);
}

DistanceField::~DistanceField()
{
	DeleteHash();
	idxGridCatalog.clear();
}

void DistanceField::partIdxGrid()
{
	//三个方向单元的数量
	int xSum = xLen / lenght + 1;
	int ySum = yLen / lenght + 1;
	int zSum = zLen / lenght + 1;

	numX = xSum;
	numY = ySum;
	numZ = zSum;

	//当前六面体单元格起点，(OpenGL，左下后)
	double curx = ox;
	double cury = oy;
	double curz = oz;

	//初始化网格索引目录
	idxGridCatalog.clear();
	idxGridCatalog.resize(xSum*ySum*zSum);

	//获取每个节点的坐标
	int ic = 0;
	for (int i = 0; i < xSum; i++)
	{
		for (int j = 0; j < ySum; j++)
		{
			for (int k = 0; k < zSum; k++)
			{
				idxGridCatalog[ic] = Point3D(curx, cury, curz);
				ic++;
				curz += lenght;//当前z坐标
			}
			curz = oz;
			cury += lenght;//当前y坐标
		}
		cury = oy;
		curx += lenght;//当前x坐标
	}
}

bool DistanceField::Partition(double &len, double xSize, double ySize, double zSize, double ox, double oy, double oz)
{
	//三个方向单元的数量
	int xSum = xSize / len + 1;
	int ySum = ySize / len + 1;
	int zSum = zSize / len + 1;

	//当前六面体单元格起点，(OpenGL，左下后)
	double curx = ox;
	double cury = oy;
	double curz = oz;

	GridIndex idxGrid;
	DistanceNode* gridNode;
	//开始划分单元格，并添加到哈希表
	for (int i = 0; i < xSum; i++)
	{
		for (int j = 0; j < ySum; j++)
		{
			for (int k = 0; k < zSum; k++)
			{
				//获取当前网格索引
				idxGrid = this->GetGridIndex(curx, cury, curz);

				//添加网格单元到哈希表
				gridNode = new DistanceNode(curx, cury, curz, idxGrid);

				//#pragma omp critical(max_arx)
				this->AddGridNode(idxGrid, gridNode);

				curz += lenght;//当前z坐标
			}
			curz = oz;
			cury += lenght;//当前y坐标
		}
		cury = oy;
		curx += lenght;//当前x坐标
	}

	return true;
}

bool DistanceField::DeleteHash()
{
	DistanceNode* gridNode = NULL;
	DistanceNode* next = NULL;
	for (int i = 0; i < hash.size(); i++)
	{
		gridNode = hash[i];
		if (gridNode != NULL)
		{
			next = gridNode->next;
			delete gridNode;
			while (next != NULL)
			{
				gridNode = next;
				next = gridNode->next;
				if (gridNode != NULL)
					delete gridNode;
			}
		}
	}
	return true;
}

bool DistanceField::ClearGrid()
{
	//三个方向单元的数量
	int xSum = xLen / lenght + 1;
	int ySum = yLen / lenght + 1;
	int zSum = zLen / lenght + 1;

	//当前六面体单元格起点，(OpenGL，左下后)
	double curx = ox;
	double cury = oy;
	double curz = oz;

	GridIndex idxGrid;
	DistanceNode* gridNode;
	//根据每个节点的坐标，遍历每个节点
	for (int i = 0; i < xSum; i++)
	{
		for (int j = 0; j < ySum; j++)
		{
			for (int k = 0; k < zSum; k++)
			{
				//获取网格单元
				gridNode = GetGridNode(curx, cury, curz);

				gridNode->data = MAX_DIS;

				curz += lenght;//当前z坐标
			}
			curz = oz;
			cury += lenght;//当前y坐标
		}
		cury = oy;
		curx += lenght;//当前x坐标
	}
	return true;
}

/********************查找*********************/
int DistanceField::GetHshSize()
{
	return hash.size();
}

GridIndex DistanceField::GetGridIndex(double x, double y, double z)
{
	//返回网格索引，向下取整(x/len)..
	return GridIndex(floor(float(x / lenght)), floor(float(y / lenght)), floor(float(z / lenght)));
}

//获取哈希值，输入坐标位置
int DistanceField::GetHashIndex(double x, double y, double z)
{
	if (hash.empty())
	{
		return -1;
	}	//网格索引
	int ix = floor(float(x / lenght));
	int iy = floor(float(y / lenght));
	int iz = floor(float(z / lenght));

	//哈希函数
	return (ix*P1 ^ iy*P2 ^ iz*P3) % hash.size();
}

//获取哈希值，通过网格索引
int DistanceField::GetHashIndex(GridIndex &idxGrid)
{
	if (hash.empty())
	{
		return -1;
	}
	//网格索引
	int ix = idxGrid.x;
	int iy = idxGrid.y;
	int iz = idxGrid.z;

	return (ix*P1 ^ iy*P2 ^ iz*P3) % hash.size();
}

DistanceNode* DistanceField::GetGridNode(double x, double y, double z)
{
	//根据网格索引获取哈希索引
	int hashIndex = GetHashIndex(x, y, z);
	DistanceNode* gridNode = hash[hashIndex];

	DistanceNode goalNode(x, y, z);
	//直到找到对应的网格单元
	while (gridNode != NULL && *gridNode != goalNode)
	{
		gridNode = gridNode->next;
	}

	if (gridNode == NULL)
		return gridNode;

	return gridNode;
}

double DistanceField::GetNodeData(double x, double y, double z)
{
	return GetGridNode(x, y, z)->data;
}

/********************添加*********************/
bool DistanceField::AddGridNode(GridIndex idxGrid, DistanceNode* gridNode)
{
	if (gridNode == NULL)
		return false;
	int idxHash = GetHashIndex(idxGrid);
	if (idxHash > hash.size())
	{
		return false;
	}

	/*在当前哈希节点的链表中，找到表尾，插入此网格单元*/
	DistanceNode* curNode = NULL;
	curNode = hash[idxHash];
	if (curNode == NULL)
	{
		hash[idxHash] = gridNode;

		numNode++;
		return true;
	}
	//遍历节点链表，找到表尾
	while (curNode->next != NULL)
	{
		curNode = curNode->next;
	}
	curNode->next = gridNode;

	numNode++;

	return true;
}

bool DistanceField::AddData(double x, double y, double z, double data)
{
	DistanceNode* gridNode = this->GetGridNode(x, y, z);

	if (gridNode == NULL)
		return false;

	gridNode->data = data;

	return true;
}

#define max(a,b) (b<a?a:b)

//这里的形参修正为当前x, y, z 所对应的流体邻居
double DistanceField::GetDisValue(double x, double y, double z, const vector<IISPHParticle*> &fluidNeighbors, double radius)
{
	double ri = radius;//每个粒子的半径
	const double SR = 8 * ri;//作用域半径

	//当前粒子的邻居
	
	if (fluidNeighbors.empty())
	{
		return MAX_DIS;
	}

	double _x = 0.0;
	double _y = 0.0;
	double _z = 0.0;
	double _r = 0.0;
	double R = SR;//搜索邻域半径
	double wi = 0.0;//权值
	double s = 0.0;
	double distValue = MAX_DIS;

	IISPHParticle* pi;//粒子i
	IISPHParticle* pj;//粒子j
	//遍历邻域内的每个粒子，计算加权平均值
	//#pragma omp parallel for reduction(+:_x,_y,_z,_r) 
	for (int i = 0; i < fluidNeighbors.size(); i++)
	{
		pi = fluidNeighbors[i];

		double uws = abs(sqrt(pow((pi->position.x - x), 2) + pow((pi->position.y - y), 2) + pow((pi->position.z - z), 2))) / R;
		double uw = max(0, pow((1 - uws*uws), 3));

		double dw = 0.0;
		//#pragma omp parallel for reduction(+:dw) 
		for (int j = 0; j < fluidNeighbors.size(); j++)
		{
			pj = fluidNeighbors[j];

			double dws = abs(sqrt(pow((pj->position.x - x), 2) + pow((pj->position.y - y), 2) + pow((pj->position.z - z), 2))) / R;
			dw += max(0, pow((1 - dws*dws), 3));
		}

		wi = ((uw == 0 && dw == 0) ? 1 : (uw / dw));

		_x += wi*pi->position.x;
		_y += wi*pi->position.y;
		_z += wi*pi->position.z;

		_r += wi * ri;
	}

	distValue = abs(sqrt(pow((x - _x), 2) + pow((y - _y), 2) + pow((z - _z), 2))) - _r;

	return distValue;
}