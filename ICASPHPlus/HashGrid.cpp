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

	//���ֿռ�����
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


	//���ƹ�ϣ���ÿһ����
	this->hash.resize(hashGrid.hash.size(), NULL);
	GridCellDistance* next = NULL;
	GridCellDistance* thisNext = NULL;
	for (int i = 0; i < hashGrid.hash.size(); i++)
	{
		//��������
		if (hashGrid.hash[i] != NULL)
		{
			this->hash[i] = new GridCellDistance(hashGrid.hash[i]);
			next = hashGrid.hash[i]->next;
			thisNext = this->hash[i];
		}
		//�����������
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
	//��������Ԫ������
	int xSum = abs(floor(ox / len)) + abs(ceil(xSize + ox) / len);
	int ySum = abs(floor(oy / len)) + abs(ceil(ySize + oy) / len);
	int zSum = abs(floor(oz / len)) + abs(ceil(zSize + oz) / len);

	numX = xSum;
	numY = ySum;
	numZ = zSum;

	//��ǰ�����嵥Ԫ����㣬(OpenGL�����º�)
	double curx = ox;
	double cury = oy;
	double curz = oz;

	GridIndex idxGrid;
	GridCellDistance* gridCellDistance;
	//��ʼ���ֵ�Ԫ�񣬲���ӵ���ϣ��
	for (int i = 0; i <= xSum; i++)
	{
		for (int j = 0; j <= ySum; j++)
		{
			for (int k = 0; k <= zSum; k++)
			{
				//��ȡ��ǰ��������
				idxGrid = this->GetGridIndex(curx, cury, curz);

				//�������Ԫ����ϣ��
				gridCellDistance = new GridCellDistance(curx, cury, curz, idxGrid);
				this->AddHashCell(idxGrid, gridCellDistance);
				curz += len;//��ǰz����
			}
			curz = oz;
			cury += len;//��ǰy����
		}
		cury = oy;
		curx += len;//��ǰx����
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
	//��������Ԫ������
	int xSum = abs(floor(ox / lenght)) + abs(ceil(xLen + ox) / lenght);
	int ySum = abs(floor(oy / lenght)) + abs(ceil(yLen + oy) / lenght);
	int zSum = abs(floor(oz / lenght)) + abs(ceil(zLen + oz) / lenght);

	//��ǰ�����嵥Ԫ����㣬(OpenGL�����º�)
	double curx = ox;
	double cury = oy;
	double curz = oz;

	GridIndex idxGrid;
	GridCellDistance* gridCell;
	//��ʼ���ֵ�Ԫ�񣬲���ӵ���ϣ��
	for (int i = 0; i <= xSum; i++)
	{
		for (int j = 0; j <= ySum; j++)
		{
			for (int k = 0; k <= zSum; k++)
			{
				//��ȡ��ǰ��������
				idxGrid = this->GetGridIndex(curx, cury, curz);
				//��ȡ����Ԫ
				gridCell = GetGridCell(idxGrid);

				//numParticle -= gridCell->idxParticle.size();

				//gridCell->idxParticle.clear();

				curz += lenght;//��ǰz����
			}
			curz = oz;
			cury += lenght;//��ǰy����
		}
		cury = oy;
		curx += lenght;//��ǰx����
	}


	return true;
}

/********************����*********************/
int HashGrid::GetHshSize()
{
	return hash.size();
}

GridIndex HashGrid::GetGridIndex(double x, double y, double z)
{
	//������������������ȡ��(x/len)..
	return GridIndex(floor(float(x / lenght)), floor(float(y / lenght)), floor(float(z / lenght)));
}

//��ȡ��ϣֵ����������λ��
int HashGrid::GetHashIndex(double x, double y, double z)
{
	if (hash.empty())
	{
		return -1;
	}
	//��������
	int ix = floor(float(x / lenght));
	int iy = floor(float(y / lenght));
	int iz = floor(float(z / lenght));

	//��ϣ����
	return (ix*P1 ^ iy*P2 ^ iz*P3) % hash.size();
}

//��ȡ��ϣֵ��ͨ����������
int HashGrid::GetHashIndex(GridIndex &idxGrid)
{
	if (hash.empty())
	{
		return -1;
	}
	//��������
	int ix = idxGrid.x;
	int iy = idxGrid.y;
	int iz = idxGrid.z;

	return (ix*P1 ^ iy*P2 ^ iz*P3) % hash.size();
}


GridCellDistance* HashGrid::GetGridCell(GridIndex &idxGrid)
{
	a++;
	//��������������ȡ��ϣ����
	int hashIndex = GetHashIndex(idxGrid);
	//if (hashIndex == -1)
	//	return NULL;
	//if (hashIndex >= hash.size())
	//	return NULL;
	GridCellDistance* gridCell = hash[hashIndex];
	b++;
	
	//ֱ���ҵ���Ӧ������Ԫ
	while (gridCell != NULL && gridCell->gridIndex != idxGrid)
	{
		gridCell = gridCell->next;
	}
	c++;
	if (gridCell == NULL)
		return gridCell;

	return gridCell;
}

/********************���*********************/
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
//	/*�ڵ�ǰ��ϣ�ڵ�������У��ҵ���β�����������Ԫ*/
//	GridCell* curCell = NULL;
//	curCell = hash[idxHash];
//	if (curCell == NULL)
//	{
//		curCell = gridCell;
//		return true;
//	}
//	//�����ڵ������ҵ���β
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
//	/*�ڵ�ǰ��ϣ�ڵ�������У��ҵ���β�����������Ԫ*/
//	GridCell* curCell = NULL;
//	curCell = hash[hashIndex];
//	if (curCell == NULL)
//	{
//		curCell = gridCell;
//		return true;
//	}
//	//�����ڵ������ҵ���β
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

	/*�ڵ�ǰ��ϣ�ڵ�������У��ҵ���β�����������Ԫ*/
	GridCellDistance* curCell = NULL;
	curCell = hash[idxHash];
	if (curCell == NULL)
	{
		hash[idxHash] = gridCell;

		numCell++;
		return true;
	}
	//�����ڵ������ҵ���β
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
//	//�����������꣬������Χ����Ԫ,3*3*3
//	int idxParticles;
//	double radii2 = radii*radii;
//	double dist = 0.0;
//	GridIndex idxGrid = this->GetGridIndex(x, y, z);//ԭ��������
//	for (int ixc = -1; ixc <= 1; ixc++)
//	{
//		for (int iyc = -1; iyc <= 1; iyc++)
//		{
//			for (int izc = -1; izc <= 1; izc++)
//			{
//				GridCell* gridCell = this->GetGridCell(GridIndex(idxGrid.x + ixc, idxGrid.y + iyc, idxGrid.z + izc));//��ǰ����
//				if (gridCell != NULL)
//				{
//					for (int icell = 0; icell < gridCell->idxParticle.size(); icell++)
//					{
//						idxParticles = gridCell->idxParticle[icell];
//						Particle pj = particles[idxParticles];//ȡ��ѭ������
//						dist = pow(x - pj.pos.x, 2) + pow(y - pj.pos.y, 2) + pow(z - pj.pos.z, 2);//�����������ӵľ���ƽ��
//						//����������Ӿ����ڷ�Χ�� �����ж�
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