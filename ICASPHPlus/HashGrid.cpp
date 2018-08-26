#include "HashGrid.h"
#include <QtDebug>

HashGrid::HashGrid()
{
}

HashGrid::HashGrid(double len, uint hashSize, double xSize, double ySize, double zSize, double xOrigin, double yOrigin, double zOrigin)
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

	//划分空间网格
	Partition(len, xSize, ySize, zSize, xOrigin, yOrigin, zOrigin);
}

HashGrid::~HashGrid()
{
	DeleteHash();
}

bool HashGrid::operator=(HashGrid &hashGrid)
{

	this->numCell = hashGrid.numCell;
	this->numParticle = hashGrid.numParticle;
	lenght = hashGrid.lenght;
	this->xLen = hashGrid.xLen;
	this->yLen = hashGrid.yLen;
	this->zLen = hashGrid.zLen;
	this->ox = hashGrid.ox;
	this->oy = hashGrid.oy;
	this->oz = hashGrid.oz;
	this->numX = hashGrid.numX;
	this->numY = hashGrid.numY;
	this->numZ = hashGrid.numZ;


	//复制哈希表的每一项，深拷贝
	this->hash.resize(hashGrid.hash.size(), NULL);
	GridCellDistance* next = NULL;
	GridCellDistance* thisNext = NULL;
	for (int i = 0; i < hashGrid.hash.size(); i++)
	{
		//拷贝表首
		if (hashGrid.hash[i] != NULL)
		{
			this->hash[i] = new GridCellDistance(hashGrid.hash[i]);
			next = hashGrid.hash[i]->next;
			thisNext = this->hash[i];
		}
		//拷贝后续结点
		while (next != NULL)
		{
			thisNext->next = new GridCellDistance(next);
			next = next->next;
			thisNext = thisNext->next;
		}
	}

	return true;
}

bool HashGrid::Partition(double &len, double xSize, double ySize, double zSize, double ox, double oy, double oz)
{
	//三个方向单元的数量
	int xSum = abs(floor(ox / len)) + abs(ceil(xSize + ox) / len);
	int ySum = abs(floor(oy / len)) + abs(ceil(ySize + oy) / len);
	int zSum = abs(floor(oz / len)) + abs(ceil(zSize + oz) / len);

	numX = xSum;
	numY = ySum;
	numZ = zSum;

	//当前六面体单元格起点，(OpenGL，左下后)
	double curx = ox;
	double cury = oy;
	double curz = oz;

	GridIndex idxGrid;
	GridCellDistance* gridCellDistance;
	//开始划分单元格，并添加到哈希表
	for (int i = 0; i <= xSum; i++)
	{
		for (int j = 0; j <= ySum; j++)
		{
			for (int k = 0; k <= zSum; k++)
			{
				//获取当前网格索引
				idxGrid = this->GetGridIndex(curx, cury, curz);

				//添加网格单元到哈希表
				gridCellDistance = new GridCellDistance(curx, cury, curz, idxGrid);
				this->AddHashCell(idxGrid, gridCellDistance);
				curz += len;//当前z坐标
			}
			curz = oz;
			cury += len;//当前y坐标
		}
		cury = oy;
		curx += len;//当前x坐标
	}

	return true;
}

bool HashGrid::DeleteHash()
{
	GridCellDistance* gridCell = NULL;
	GridCellDistance* next = NULL;
	for (int i = 0; i < hash.size(); i++)
	{
		gridCell = hash[i];
		if (gridCell != NULL)
		{
			next = gridCell->next;
			delete gridCell;
			while (next != NULL)
			{
				gridCell = next;
				next = gridCell->next;
				if (gridCell != NULL)
					delete gridCell;
			}
		}
	}
	return true;
}

bool HashGrid::ClearGrid()
{
	//三个方向单元的数量
	int xSum = abs(floor(ox / lenght)) + abs(ceil(xLen + ox) / lenght);
	int ySum = abs(floor(oy / lenght)) + abs(ceil(yLen + oy) / lenght);
	int zSum = abs(floor(oz / lenght)) + abs(ceil(zLen + oz) / lenght);

	//当前六面体单元格起点，(OpenGL，左下后)
	double curx = ox;
	double cury = oy;
	double curz = oz;

	GridIndex idxGrid;
	GridCellDistance* gridCell;
	//开始划分单元格，并添加到哈希表
	for (int i = 0; i <= xSum; i++)
	{
		for (int j = 0; j <= ySum; j++)
		{
			for (int k = 0; k <= zSum; k++)
			{
				//获取当前网格索引
				idxGrid = this->GetGridIndex(curx, cury, curz);
				//获取网格单元
				gridCell = GetGridCell(idxGrid);

				//numParticle -= gridCell->idxParticle.size();

				//gridCell->idxParticle.clear();

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
int HashGrid::GetHshSize()
{
	return hash.size();
}

GridIndex HashGrid::GetGridIndex(double x, double y, double z)
{
	//返回网格索引，向下取整(x/len)..
	return GridIndex(floor(float(x / lenght)), floor(float(y / lenght)), floor(float(z / lenght)));
}

//获取哈希值，输入坐标位置
int HashGrid::GetHashIndex(double x, double y, double z)
{
	if (hash.empty())
	{
		return -1;
	}
	//网格索引
	int ix = floor(float(x / lenght));
	int iy = floor(float(y / lenght));
	int iz = floor(float(z / lenght));

	//哈希函数
	return (ix*P1 ^ iy*P2 ^ iz*P3) % hash.size();
}

//获取哈希值，通过网格索引
int HashGrid::GetHashIndex(GridIndex &idxGrid)
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


GridCellDistance* HashGrid::GetGridCell(GridIndex &idxGrid)
{
	a++;
	//根据网格索引获取哈希索引
	int hashIndex = GetHashIndex(idxGrid);
	//if (hashIndex == -1)
	//	return NULL;
	//if (hashIndex >= hash.size())
	//	return NULL;
	GridCellDistance* gridCell = hash[hashIndex];
	b++;
	
	//直到找到对应的网格单元
	while (gridCell != NULL && gridCell->gridIndex != idxGrid)
	{
		gridCell = gridCell->next;
	}
	c++;
	if (gridCell == NULL)
		return gridCell;

	return gridCell;
}

/********************添加*********************/
//bool HashGrid::AddHashCell(double x, double y, double z, GridCell* gridCell)
//{
//	if (gridCell == NULL)
//		return false;
//
//	int idxHash = GetHashIndex(x, y, x);
//	if (idxHash > hash.size())
//	{
//		return false;
//	}
//
//	/*在当前哈希节点的链表中，找到表尾，插入此网格单元*/
//	GridCell* curCell = NULL;
//	curCell = hash[idxHash];
//	if (curCell == NULL)
//	{
//		curCell = gridCell;
//		return true;
//	}
//	//遍历节点链表，找到表尾
//	while (curCell->next != NULL)
//	{
//		curCell = curCell->next;
//	}
//	curCell->next = gridCell;
//
//	return true;
//}
//bool HashGrid::AddHashCell(int hashIndex, GridCell* gridCell)
//{
//	if (gridCell == NULL)
//		return false;
//	if (hashIndex > hash.size())
//	{
//		return false;
//	}
//
//	/*在当前哈希节点的链表中，找到表尾，插入此网格单元*/
//	GridCell* curCell = NULL;
//	curCell = hash[hashIndex];
//	if (curCell == NULL)
//	{
//		curCell = gridCell;
//		return true;
//	}
//	//遍历节点链表，找到表尾
//	while (curCell->next != NULL)
//	{
//		curCell = curCell->next;
//	}
//	curCell->next = gridCell;
//
//	return true;
//}
bool HashGrid::AddHashCell(GridIndex idxGrid, GridCellDistance* gridCell)
{
	if (gridCell == NULL)
		return false;
	int idxHash = GetHashIndex(idxGrid);
	if (idxHash > hash.size())
	{
		return false;
	}

	/*在当前哈希节点的链表中，找到表尾，插入此网格单元*/
	GridCellDistance* curCell = NULL;
	curCell = hash[idxHash];
	if (curCell == NULL)
	{
		hash[idxHash] = gridCell;

		numCell++;
		return true;
	}
	//遍历节点链表，找到表尾
	while (curCell->next != NULL)
	{
		curCell = curCell->next;
	}
	curCell->next = gridCell;

	numCell++;

	return true;
}

//bool HashGrid::AddParticle(double x, double y, double z, int idx)
//{
//	d++;
//	GridIndex idxGrid = this->GetGridIndex(x, y, z);
//	e++;
//	GridCellDistance* gridCell = this->GetGridCell(idxGrid);
//	f++;
//	if (gridCell == NULL)
//		return false;
//
//	gridCell->idxParticle.push_back(idx);
//
//
//	numParticle++;
//
//	return true;
//}

//bool HashGrid::SearchAdjParticles(double x, double y, double z, double radii, const vector<Particle> &particles, vector<int> &adjParticles)
//{
//	adjParticles.clear();
//	//根据粒子坐标，遍历周围网格单元,3*3*3
//	int idxParticles;
//	double radii2 = radii*radii;
//	double dist = 0.0;
//	GridIndex idxGrid = this->GetGridIndex(x, y, z);//原网格索引
//	for (int ixc = -1; ixc <= 1; ixc++)
//	{
//		for (int iyc = -1; iyc <= 1; iyc++)
//		{
//			for (int izc = -1; izc <= 1; izc++)
//			{
//				GridCell* gridCell = this->GetGridCell(GridIndex(idxGrid.x + ixc, idxGrid.y + iyc, idxGrid.z + izc));//当前网格
//				if (gridCell != NULL)
//				{
//					for (int icell = 0; icell < gridCell->idxParticle.size(); icell++)
//					{
//						idxParticles = gridCell->idxParticle[icell];
//						Particle pj = particles[idxParticles];//取内循环粒子
//						dist = pow(x - pj.pos.x, 2) + pow(y - pj.pos.y, 2) + pow(z - pj.pos.z, 2);//计算两个粒子的距离平方
//						//如果两个粒子距离在范围内 进入判断
//						if (dist < radii2 && dist != 0.0)
//						{
//							adjParticles.push_back(idxParticles);
//						}
//					}
//				}
//			}
//		}
//	}
//
//
//	return true;
//}