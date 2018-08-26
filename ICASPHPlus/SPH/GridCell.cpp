#include "GridCell.h"

void GridCell::Initialization(int i, int j, int k, int hashV)
{
	ComputeCellLength(); //���㵱ǰ���񳤶�
	ResetParticleList(); //���������б�
	SetIndexVec(i, j, k); //���ø������λ����
	SetHashIndex(hashV); //����������hashIndex
}

void GridCell::ComputeCellLength()
{
	//��ȡȫ�ֱ��������õ�ǰ����
	gridCellLength = SharedData::GetCellLength();
}

void GridCell::PushParticle(WCSPHParticle* a)
{
	particleList.push_back(a);
}

void GridCell::PushParticle(IISPHParticle* a)
{
	IISPHparticleList.push_back(a);
}

void GridCell::PushBoundaryParticle(StaticRigidParticle* a)
{
	boundaryList.push_back(a); //���߽�����װ���Ӧ�ı߽�������뵱��
}

void GridCell::PushDynamicParticle(DynamicRigidParticle * a)
{
	dynamicParticleList.push_back(a);
}

void GridCell::ResetParticleList()
{
	particleList.clear();
	IISPHparticleList.clear();
}

void GridCell::ResetBoundaryParticleList()
{
	boundaryList.clear();
}

void GridCell::ResetDynamicParticleList()
{
	dynamicParticleList.clear();
}

vector<WCSPHParticle*> GridCell::GetParticleList()
{
	return particleList;
}

vector<IISPHParticle*> GridCell::GetIISPHParticleList()
{
	return IISPHparticleList;
}

vector<StaticRigidParticle*> GridCell::GetBoundaryParticleList()
{
	return boundaryList;
}

vector<DynamicRigidParticle*> GridCell::GetDynamicParticleList()
{
	return dynamicParticleList;
}

bool GridCell::CheckGridCell(Vector3f parPosition)
{
	//����ȫ�ֱ�����ȷ���߽��ƫ��
	double temp = SharedData::GetOffset(); 
	//�õ���ǰ���ӵ�λ��
	Vector3i indexVecPos = FunctionKit::PositionMapIndex(parPosition, gridCellLength, temp);
	//��⵱ǰ�������ڵ�λ�ú͵�ǰ�����λ����һ�µ�
	return CheckGridCell(indexVecPos); //������һ�������غ���
}

//ͨ��λ�ñ������ѯ
bool GridCell::CheckGridCell(Vector3i parIndex)
{
	if (indexVec.GetX() == parIndex.GetX() &&
		indexVec.GetY() == parIndex.GetY() &&
		indexVec.GetZ() == parIndex.GetZ())
	{
		return true;
	}
	return false;
}

bool GridCell::CheckGridCell(WCSPHParticle particle)
{
	return CheckGridCell(particle.position); //������һ�������غ���
}

bool GridCell::CheckGridCell(IISPHParticle* particle)
{
	return CheckGridCell(particle->position); 
}

bool GridCell::CheckGridCell(DynamicRigidParticle * particle)
{
	return CheckGridCell(particle->m_x);
}

GridCell& GridCell::operator=(const GridCell & a)
{
	indexVec = a.indexVec;
	hashIndex = a.hashIndex;
	gridCellLength = a.gridCellLength;
	ResetParticleList();
	for (int i = 0; i < a.particleList.size(); i++)
	{
		PushParticle(a.particleList[i]);
	}
	//���бʼǣ�������ĸ�ֵ����������չ֮��û�н��и�ֵ�߽����ӣ������޷���ѯ����Χ������
	ResetBoundaryParticleList();
	for (int i = 0; i < a.boundaryList.size(); i++)
	{
		PushBoundaryParticle(a.boundaryList[i]);
	}
	ResetDynamicParticleList();
	for (int i=0; i<a.dynamicParticleList.size();i++)
	{
		PushDynamicParticle(a.dynamicParticleList[i]);
	}

	return *this;
}

void GridCell::SetGridCellLength(double length)
{
	gridCellLength = length;
}

double GridCell::GetGridCellLength()
{
	return gridCellLength;
}

void GridCell::SetIndexVec(Vector3i value)
{
	indexVec = value;
}

void GridCell::SetIndexVec(int i, int j, int k)
{
	Vector3i temp(i, j, k);
	indexVec = temp;
}

Vector3i GridCell::GetIndexVec()
{
	return indexVec;
}

void GridCell::SetHashIndex(int value)
{
	hashIndex = value;
}

int GridCell::GetHashIndex()
{
	return hashIndex;
}

