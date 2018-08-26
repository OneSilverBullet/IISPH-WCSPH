#pragma once
#include <vector>
#include "HGrid.h"
#include "HashGrid.h"
#include "Point3D.h"
#include "../ICASPHPlus/SPH/IISPHParticle.h"
//#include "Particle.h"
//using namespace std;

/***************
*��ϣ��������ڵ���.
*ÿ���������һ������ڵ㡣
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
	double data;	// �ڵ����ֵ
	double x, y, z;			// ����
	GridIndex gridIndex;		// �ڵ�����
	DistanceNode *next = NULL;
};

/***************
*��ϣ������
*/
class DistanceField : public HGrid
{
public:
	DistanceField();
	DistanceField(double len, uint hashSize, double xSize, double ySize, double zSize, double xOrigin = 0.0, double yOrigin = 0.0, double zOrigin = 0.0);
	~DistanceField();

public:
	int numNode = 0;	//�ڵ�����

public:
	vector<DistanceNode*> hash;//�ڵ��ϣ��
	vector<Point3D> idxGridCatalog;//��������Ŀ¼

	void partIdxGrid();

public:
	//��������Ԫ���ȣ��Ϳռ�xyz���򳤶ȡ���㣬�����볡����Ϊ�ռ�����ȷ�������ڵ�
	bool Partition(double &len, double xSize, double ySize, double zSize, double ox, double oy, double oz);
	//ɾ����ϣ����ÿ������ڵ�
	bool DeleteHash();
	//�������ռ���볡
	bool ClearGrid();

public:
	//��ȡ�����������꣬��������ֵ
	GridIndex GetGridIndex(double x, double y, double z);
	//��ȡ��ϣֵ����������λ��
	int GetHashIndex(double x, double y, double z);
	//��ȡ��ϣֵ��ͨ����������
	int GetHashIndex(GridIndex &idxGrid);
	//��ȡ��ϣ��
	int GetHshSize();
	//ͨ����ϣ���ȡ����ڵ�
	DistanceNode* GetGridNode(double x, double y, double z);
	//��ȡ�ڵ����ݣ�����ֵ
	double GetNodeData(double x, double y, double z);

public:
	//������ڵ���ӵ���ϣ��
	bool AddGridNode(GridIndex idxGrid, DistanceNode* gridNode);
	//���ڵ����ֵ��ӵ�����ڵ�
	bool AddData(double x, double y, double z, double data);
	//��ȡ�������ֵ
	double GetDisValue(double x, double y, double z, const vector<IISPHParticle*> &particles,  double radius = 1.0);
};

