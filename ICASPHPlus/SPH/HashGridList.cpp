#include "HashGridList.h"

//������GridCellListʹ��Ĭ�ϵĲ���
void HashGridList::Intialization()
{
	ClearBoundaryParticle();
	ComputeGridWLH(cellLength, boundary);    //����õ�WLH�������������
	ComputeCellNum();    //�ټ���cell������
	InitializeHashCellNun(); //��ʼ����ǰ��ϣ�����б�ĵ���������
	EnsureTheGrid();     //�ڵ�ǰ�Ļ�����ȷ������������������
	//�������е����ݴ��뵽��Ӧ�ĵ�ȫ�ֲ�������
	SharedData::SetCellLength(cellLength);
	SharedData::SetHashCellNum(hashCellNum);
}

void HashGridList::Intialization(double cellL, double thred, Boundary bound)
{
	//��ʼ������
	cellLength = cellL;
	boundary = bound;
	hashThred = thred;
	ComputeGridWLH(cellLength, boundary);    //����õ�WLH�������������
	ComputeCellNum();    //�ټ���cell������
	InitializeHashCellNun(); //��ʼ����ǰ��ϣ�����б�ĵ���������
	EnsureTheGrid();     //�ڵ�ǰ�Ļ�����ȷ������������������
}

//���ݱ߽��Լ��µ�cellLength����õ���X,Y,Z�����cell����
void HashGridList::ComputeGridWLH(double cellL, Boundary bound)
{
	gridWLH.x = ceil(double(bound.gridBoundaryWidth / cellL));
	gridWLH.y = ceil (double(bound.gridBoundaryLength / cellL));
	gridWLH.z = ceil (double(bound.gridBoundaryHeight / cellL));
}

//�����ú�������õ����յĿռ���������
void HashGridList::ComputeCellNum()
{
	cellNum = gridWLH.x*gridWLH.y*gridWLH.z; //����gridWLH�ĳ���߷����cell������˵õ��������������
}

//��ϣ��Ĵ�С��ʼ��Ϊ��ǰ��������֮��ĵ�һ������
void HashGridList::InitializeHashCellNun()
{
	int temp = cellNum;
	for (;;)
	{
		//�����ǰֵΪ��������ô��ֵ���������
		if (FunctionKit::CheckPrime(temp)) 
		{
			hashCellNum = temp;
			//����ע�ͣ�����Ķ�̬������һ������Ϊ3200�򣬲����������
			//����������������ˮ��ĳ����Լ��ռ�ĳ���
			gridCellList.resize(hashCellNum);
			return;
		}
		temp++; //����������У�ֱ���ҵ�һ������
	}
}

//�������ֻ��Ҫ�ڳ�ʼ����ʱ������һ��
void HashGridList::EnsureTheGrid()
{
	GenerateGridCell(); //���ȸ��ݵ�ǰ������Կռ���л���
	//��������ѭ����ֱ���ҵ�һ�����ʵĹ�ϣ��
	for (;;)
	{
		//��⵱ǰ��ϣ���Ƿ����Ҫ��
		if (CheckThred()) 
		{
			return; //�������Ҫ������ô�ͽ�������
		}
		else //��ǰ��ϣ������Ҫ��
		{
			Rehash(); 
			GenerateGridCell(); //rehash֮������������������
		}
	}
	return;
}

void HashGridList::GenerateGridCell()
{
	//�����������г�ʼ��
	for (int i = 0; i < gridWLH.x; i++)
	{
		for (int j = 0; j < gridWLH.y; j++)
		{
			for (int k = 0; k < gridWLH.z; k++)
			{
				//ȡ�õ�ǰλ�õĹ�ϣ��ֵ
				int hashValue = FunctionKit::GridHash(i, j, k, hashCellNum);
				GridCell grid(i, j, k, hashValue);
				gridCellList[hashValue].push_back(grid);
			}
		}
	}
	
}

//���������Ŀ���ǣ�ȡ��ԭ����hashֵ�Ķ���֮��ĵ�һ������
void HashGridList::Rehash()
{
	int hashValue = hashCellNum * 2; 
	for (;;)
	{
		if (FunctionKit::CheckPrime(hashValue))
		{
			hashCellNum = hashValue; //�����ǰhashֵΪ��������ô�͸���
			gridCellList.resize(hashCellNum);
			return; //��������
		}
		hashValue++; //���ϵ�����hashֵ
	}
}

//��������ʱ�Ĺ�ϣ���Ƿ���Ҫrehash
bool HashGridList::CheckThred()
{
	int storeNum = 0; //��ϣ�����Ѿ��洢������
	for (int i = 0; i < gridCellList.size(); i++)
	{
		if (gridCellList[i].size() != 0)
		{
			storeNum++; //һ����ǰ�Ŀ�λ��Ϊ�գ���ô�ͽ�storeNum����
		}
	}
	//���бʼǣ��˴��������жϵĴ���ԭ���Ǳ�ĳ���Ϊ0,�Ӷ����³���0�Ĵ���
	//�����ʽ��ͨ�������ϣ�������������˱�����ʾ���Ե�����
	float ratio = float(float(storeNum) / gridCellList.size()); 
	if (ratio > hashThred) //װ�����ӳ�����ֵ����ô��������λ��
	{
		return false;  //����false��ʱ�򣬴�ʱ��ϣ����Ҫ��������λ��
	}
	return true; //����true��ʱ�򣬹�ϣ������������
}

void HashGridList::PushParticle(WCSPHParticle* a)
{
	int particleHash = a->hashIndex; //�õ���ǰ���ӵĹ�ϣֵ
	vector<GridCell> tempList = gridCellList[particleHash]; //�õ���ǰ���Ӷ�Ӧ���б�
	//�õ��������ڵ�λ��index
	Vector3i particleIndex = FunctionKit::PositionMapIndex(a->position, cellLength, boundary.offset);
	for (int i = 0; i < tempList.size(); i++)
	{
		//�����ǰ�����ӺźͶ�Ӧ������ƥ�䣬��ôȷ�����񣬽�������
		if (tempList[i].CheckGridCell(particleIndex))
		{
			//tempList[i].PushParticle(a);
			//�������뿽������
			gridCellList[particleHash][i].PushParticle(a);
			return;
		}
	}
}

void HashGridList::PushParticle(StaticRigidParticle* a)
{
	int particleHash = a->hashValue; //�õ���ǰ���ӵĹ�ϣֵ
	vector<GridCell> tempList = gridCellList[particleHash]; //�õ���ǰ���Ӷ�Ӧ���б�
   //�õ��������ڵ�λ��index
	Vector3i particleIndex = FunctionKit::PositionMapIndex(a->position, cellLength, boundary.offset);
	for (int i = 0; i < tempList.size(); i++)
	{
		//�����ǰ�����ӺźͶ�Ӧ������ƥ�䣬��ôȷ�����񣬽�������
		if (tempList[i].CheckGridCell(particleIndex))
		{
			//tempList[i].PushParticle(a);
			//�������뿽������
			gridCellList[particleHash][i].PushBoundaryParticle(a);
			//	cout << gridCellList[13651][0].GetBoundaryParticleList().size() << endl;
			return;
		}
	}
	return; //���û���ҵ�������Ҫ�������ã�****����****
}

void HashGridList::PushParticle(IISPHParticle* a)
{
	int particleHash = a->hashIndex; //�õ���ǰ���ӵĹ�ϣֵ
	vector<GridCell> tempList = gridCellList[particleHash]; //�õ���ǰ���Ӷ�Ӧ���б�
															//�õ��������ڵ�λ��index
	Vector3i particleIndex = FunctionKit::PositionMapIndex(a->position, cellLength, boundary.offset);
	for (int i = 0; i < tempList.size(); i++)
	{
		//�����ǰ�����ӺźͶ�Ӧ������ƥ�䣬��ôȷ�����񣬽�������
		if (tempList[i].CheckGridCell(particleIndex))
		{
			gridCellList[particleHash][i].PushParticle(a);
			return;
		}
	}
}

void HashGridList::PushParticle(DynamicRigidParticle * a)
{
	int particleHash = a->hashValue; //�õ���ǰ���ӵĹ�ϣֵ
	vector<GridCell> tempList = gridCellList[particleHash]; //�õ���ǰ���Ӷ�Ӧ���б�
	//��ǰ������Ҫ�õ���ǰ�Ķ�̬�������ӵĹ�ϣֵ�ſ���
	for (int i = 0; i < tempList.size(); i++)
	{
		//�����ǰ�����ӺźͶ�Ӧ������ƥ�䣬��ôȷ�����񣬽�������
		if (tempList[i].CheckGridCell(a))
		{
			gridCellList[particleHash][i].PushDynamicParticle(a);
			return;
		}
	}
}

void HashGridList::ClearParticle()
{
	for (int i = 0; i < gridCellList.size(); i++)
	{
		if (gridCellList[i].size() != 0)
		{
			for (int j = 0; j < gridCellList[i].size(); j++)
			{
				gridCellList[i][j].ResetParticleList(); //����ǰ�������б�
			}
		}
	}
}

void HashGridList::ClearBoundaryParticle()
{
	for (int i = 0; i < gridCellList.size(); i++)
	{
		if (gridCellList[i].size() != 0)
		{
			for (int j = 0; j < gridCellList[i].size(); j++)
			{
				gridCellList[i][j].ResetBoundaryParticleList(); //����ǰ�������б�
			}
		}
	}
}

double HashGridList::GetCellLength()
{
	return cellLength;
}

int HashGridList::GetCellNum()
{
	return cellNum;
}

int HashGridList::GetHashCellNum()
{
	//�õ���ǰ��ϣ��Ĵ�С
	return hashCellNum;
}

vector<GridCell> HashGridList::GetGridCellVector(int hashIndex)
{
	return gridCellList[hashIndex];
}

vector<GridCell> HashGridList::GetGridCellVector(Vector3i positionIndex)
{
	int hashIndex = FunctionKit::GridHash(positionIndex.GetX(), positionIndex.GetY(), positionIndex.GetZ(), hashCellNum);
	return GetGridCellVector(hashIndex); //������һ�����غ���
}

vector<GridCell> HashGridList::GetGridCellVector(Vector3f position)
{
	//����λ��ӳ�䵽λ��Index
	Vector3i positionIndex = FunctionKit::PositionMapIndex(position, cellLength, boundary.offset); 
	return GetGridCellVector(positionIndex); //������һ�����غ���
}

GridCell HashGridList::GetGridCellWithSecVec(vector<GridCell>& gridVec, Vector3f parPosition)
{
	//���бʼǣ���������Ӧֵ�����ض�Ӧ������ֵ***
	int i;
	for (i = 0; i < gridVec.size(); i++)
	{
		if (gridVec[i].CheckGridCell(parPosition)) //���ݵ�ǰ���Ӷ�λ��Ӧ������
		{
			return gridVec[i]; //���ض�Ӧ������
		}
	}
	return GridCell();  //���򷵻�Ĭ�ϵ�����
}

GridCell HashGridList::GetGridCellWithSecVec(vector<GridCell>& gridVec, WCSPHParticle a)
{
	return GetGridCellWithSecVec(gridVec, a.position); //������һ�����غ������������Ĳ���
}

GridCell HashGridList::GetGridCell(WCSPHParticle a)
{
	vector<GridCell> temp = GetGridCellVector(a.position); //�õ���ǰ����a���ڹ�ϣ��λ��������
	return GetGridCellWithSecVec(temp, a); //��������ĺ�������������ӽ��в�ѯ�Ĺ���
}

GridCell HashGridList::GetGridCell(Vector3f parPosition)
{
	vector<GridCell> temp = GetGridCellVector(parPosition); //���ݵ�ǰ���ӵ�λ���ҵ���Ӧ�Ĺ�ϣ�����е�����
	return GetGridCellWithSecVec(temp, parPosition); //��������ĺ�����ɲ���
}


