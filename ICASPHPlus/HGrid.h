#pragma once
#include <vector>

using namespace std;

/***************
*�������ͳ���
*/
#define P1 73856093
#define P2 19349663
#define P3 83492791

#define uint unsigned

#define MAX_DIS 65536

/***************
*���������࣬ÿ��������Ӧ�ռ��е�һ������Ԫ
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
*��ϣ������࣬�����࣬�������÷���
*/
class HGrid
{
public:
	HGrid();
	~HGrid();

public:
	double lenght = 0.1;//�ڵ��ࣨ����Ԫ�߳���
	double xLen;	//�ռ�����߳�
	double yLen;
	double zLen;
	double ox;		//�ռ��������
	double oy;
	double oz;
	int numX = 0;		//����Ԫ������
	int numY = 0;
	int numZ = 0;

public:
	//��������Ԫ���ȣ��Ϳռ�xyz���򳤶ȡ���㣬�����볡����Ϊ�ռ�����ȷ�������ڵ�
	virtual bool Partition(double &len, double xSize, double ySize, double zSize, double ox, double oy, double oz) = 0;
	//ɾ����ϣ����ÿ������ڵ�
	virtual bool DeleteHash() = 0;
	//�������ռ��е�Ԫ��
	virtual bool ClearGrid() = 0;

public:
	//��ȡ�����������꣬��������ֵ
	virtual GridIndex GetGridIndex(double x, double y, double z) = 0;
	//��ȡ��ϣֵ����������λ��
	virtual int GetHashIndex(double x, double y, double z) = 0;
	//��ȡ��ϣֵ��ͨ����������
	virtual int GetHashIndex(GridIndex &idxGrid) = 0;
	//��ȡ��ϣ��
	int GetHshSize();
};

