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

	//��¼�ڵ�����
	partIdxGrid();
	//���ֿռ�����
	Partition(len, xSize, ySize, zSize, xOrigin, yOrigin, zOrigin);
}

DistanceField::~DistanceField()
{
	DeleteHash();
	idxGridCatalog.clear();
}

void DistanceField::partIdxGrid()
{
	//��������Ԫ������
	int xSum = xLen / lenght + 1;
	int ySum = yLen / lenght + 1;
	int zSum = zLen / lenght + 1;

	numX = xSum;
	numY = ySum;
	numZ = zSum;

	//��ǰ�����嵥Ԫ����㣬(OpenGL�����º�)
	double curx = ox;
	double cury = oy;
	double curz = oz;

	//��ʼ����������Ŀ¼
	idxGridCatalog.clear();
	idxGridCatalog.resize(xSum*ySum*zSum);

	//��ȡÿ���ڵ������
	int ic = 0;
	for (int i = 0; i < xSum; i++)
	{
		for (int j = 0; j < ySum; j++)
		{
			for (int k = 0; k < zSum; k++)
			{
				idxGridCatalog[ic] = Point3D(curx, cury, curz);
				ic++;
				curz += lenght;//��ǰz����
			}
			curz = oz;
			cury += lenght;//��ǰy����
		}
		cury = oy;
		curx += lenght;//��ǰx����
	}
}

bool DistanceField::Partition(double &len, double xSize, double ySize, double zSize, double ox, double oy, double oz)
{
	//��������Ԫ������
	int xSum = xSize / len + 1;
	int ySum = ySize / len + 1;
	int zSum = zSize / len + 1;

	//��ǰ�����嵥Ԫ����㣬(OpenGL�����º�)
	double curx = ox;
	double cury = oy;
	double curz = oz;

	GridIndex idxGrid;
	DistanceNode* gridNode;
	//��ʼ���ֵ�Ԫ�񣬲���ӵ���ϣ��
	for (int i = 0; i < xSum; i++)
	{
		for (int j = 0; j < ySum; j++)
		{
			for (int k = 0; k < zSum; k++)
			{
				//��ȡ��ǰ��������
				idxGrid = this->GetGridIndex(curx, cury, curz);

				//�������Ԫ����ϣ��
				gridNode = new DistanceNode(curx, cury, curz, idxGrid);

				//#pragma omp critical(max_arx)
				this->AddGridNode(idxGrid, gridNode);

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
	//��������Ԫ������
	int xSum = xLen / lenght + 1;
	int ySum = yLen / lenght + 1;
	int zSum = zLen / lenght + 1;

	//��ǰ�����嵥Ԫ����㣬(OpenGL�����º�)
	double curx = ox;
	double cury = oy;
	double curz = oz;

	GridIndex idxGrid;
	DistanceNode* gridNode;
	//����ÿ���ڵ�����꣬����ÿ���ڵ�
	for (int i = 0; i < xSum; i++)
	{
		for (int j = 0; j < ySum; j++)
		{
			for (int k = 0; k < zSum; k++)
			{
				//��ȡ����Ԫ
				gridNode = GetGridNode(curx, cury, curz);

				gridNode->data = MAX_DIS;

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
int DistanceField::GetHshSize()
{
	return hash.size();
}

GridIndex DistanceField::GetGridIndex(double x, double y, double z)
{
	//������������������ȡ��(x/len)..
	return GridIndex(floor(float(x / lenght)), floor(float(y / lenght)), floor(float(z / lenght)));
}

//��ȡ��ϣֵ����������λ��
int DistanceField::GetHashIndex(double x, double y, double z)
{
	if (hash.empty())
	{
		return -1;
	}	//��������
	int ix = floor(float(x / lenght));
	int iy = floor(float(y / lenght));
	int iz = floor(float(z / lenght));

	//��ϣ����
	return (ix*P1 ^ iy*P2 ^ iz*P3) % hash.size();
}

//��ȡ��ϣֵ��ͨ����������
int DistanceField::GetHashIndex(GridIndex &idxGrid)
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

DistanceNode* DistanceField::GetGridNode(double x, double y, double z)
{
	//��������������ȡ��ϣ����
	int hashIndex = GetHashIndex(x, y, z);
	DistanceNode* gridNode = hash[hashIndex];

	DistanceNode goalNode(x, y, z);
	//ֱ���ҵ���Ӧ������Ԫ
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

/********************���*********************/
bool DistanceField::AddGridNode(GridIndex idxGrid, DistanceNode* gridNode)
{
	if (gridNode == NULL)
		return false;
	int idxHash = GetHashIndex(idxGrid);
	if (idxHash > hash.size())
	{
		return false;
	}

	/*�ڵ�ǰ��ϣ�ڵ�������У��ҵ���β�����������Ԫ*/
	DistanceNode* curNode = NULL;
	curNode = hash[idxHash];
	if (curNode == NULL)
	{
		hash[idxHash] = gridNode;

		numNode++;
		return true;
	}
	//�����ڵ������ҵ���β
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

//������β�����Ϊ��ǰx, y, z ����Ӧ�������ھ�
double DistanceField::GetDisValue(double x, double y, double z, const vector<IISPHParticle*> &fluidNeighbors, double radius)
{
	double ri = radius;//ÿ�����ӵİ뾶
	const double SR = 8 * ri;//������뾶

	//��ǰ���ӵ��ھ�
	
	if (fluidNeighbors.empty())
	{
		return MAX_DIS;
	}

	double _x = 0.0;
	double _y = 0.0;
	double _z = 0.0;
	double _r = 0.0;
	double R = SR;//��������뾶
	double wi = 0.0;//Ȩֵ
	double s = 0.0;
	double distValue = MAX_DIS;

	IISPHParticle* pi;//����i
	IISPHParticle* pj;//����j
	//���������ڵ�ÿ�����ӣ������Ȩƽ��ֵ
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