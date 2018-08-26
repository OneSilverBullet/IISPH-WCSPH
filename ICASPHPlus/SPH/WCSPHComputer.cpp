#include "WCSPHComputer.h"
//Ŀǰ�����޸�����İ汾2018����3����30

//��ʼ��DEBUG ��ɣ�����
void WCSPHComputation::Initialization()
{
	//�Թ�ϣ�����Լ�����ģ�Ͷ�ʹ��Ĭ�ϲ������г�ʼ��
	hashGridList.boundary = boundary;
	hashGridList.Intialization(); 
	fluidModel.boundary = boundary;
	fluidModel.Initialise();
	//�����е����Ӷ�װ���Ӧ��������
	MapParticleToGrid();
	//�����еı߽�����װ���Ӧ��������
	MapBoundaryParticleToGrid();
	//���������е��ھӱ߽����ӣ�����õ��߽����ӵı߽������ھ��б�
	ComputeBoundaryNeighborList();
	//����߽����ӵ�PSIֵ
	fluidModel.ComputeBoundaryPsi(fluidModel.boundaryObj);
	//��Ҫ��һ�ν��߽�����װ�������У�������BoundaryPSIֵ
	MapBoundaryParticleToGrid();
	//�����е����������Լ��߽������������װ��֮��
	//���µ�ǰ�������ӵ����������ھӱ�
	UpdateFluidNeighborList();
	//���µ�ǰ�������Եı߽������ھӱ�
	UpdateFluidBoundaryNeighborList();

	//��WCSPHComputation���еĲ���װ��ȫ�ֱ�������
	SharedData::SetStiffnesss(stiffnesss);
	SharedData::SetGamma(gamma);
	SharedData::SetGravity(gravity);
	SharedData::SetViscosityConstant(viscosityConstant);
}

//������װ�뵽��Ӧ�Ĺ�ϣ����������
//ע�⣺�����������֮ǰ����Ҫ�����е����Ӷ����г�ʼ��
void WCSPHComputation::MapParticleToGrid()
{
	//��������е�����֮��������������µ�����
	//����������װ���Ӧλ��
	hashGridList.ClearParticle();
	for (int i = 0; i < fluidModel.particleList.size(); i++)
	{
		hashGridList.PushParticle(fluidModel.particleList[i]);
	}
}

//�߽����ӵ�����ӳ��ֻ��Ҫִ��һ�ξͺã�����Ҫÿ����������ı仯���仯
void WCSPHComputation::MapBoundaryParticleToGrid()
{
	for (int i = 0; i < fluidModel.boundaryObj.GetStaticParticleNum(); i++)
	{
		StaticRigidParticle* tempStaticParticle = fluidModel.boundaryObj.staticRigidParticleList[i];
		hashGridList.PushParticle(tempStaticParticle);
	}
	//�����еı߽����Ӵ���߽�������
}

//��֡����ʵ�����е���ֻ��Ҫ������������Ϳ�����
//�����Ȳ���ô��
void WCSPHComputation::Frame()
{
	//֡��
	int framePreSecond = FRAME_PER_SECOND;
	//���и�������ת�����õ�ÿһ֡������Ҫ�����ʱ��������Ӷ�ʵ����֡
	int CalculateDotNum = (int)(((float)1.0 / framePreSecond) / (float)(fluidModel.timeStep));
	int i = CalculateDotNum;
	for (;;)
	{
		Computation();
		if (i < 0)
		{
			break;
		}
		i--;
	}
}

//ע������������еļ���˳�򣬲��ɱ�
void WCSPHComputation::Computation()
{
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int i = 0; i < fluidModel.particleList.size(); i++)
		{
			//�õ��������ӵ����������ٽ��б�
			vector<WCSPHParticle*> neighborList = fluidModel.neighborList[i];
			//�õ��������ӵĸ��������ٽ��б�
			vector<StaticRigidParticle*> boundaryNeighborList = fluidModel.fluidBoundaryNeighborList[i];
			//�����ٽ��б�õ���ǰ���ӵ��ܶ�
			double densityT = ComputeDensity(neighborList, boundaryNeighborList, fluidModel.particleList[i]);
			//cout << densityT << endl;
			fluidModel.particleList[i]->density = densityT;
			//�����ٽ��б�õ���ǰ���ӵ�ճ����
			Vector3f viscosityA = ComputeViscostyAcceleration(neighborList, boundaryNeighborList, fluidModel.particleList[i]);
			fluidModel.particleList[i]->viscosity = viscosityA;
			//�Ե�ǰ�����ܶȽ���У����������ɶ�����P�ļ���
			double correctDensity = CorrectionParticleDensity(fluidModel.particleList[i]);
			double  P = CalculateParticleTaite(correctDensity);
			fluidModel.particleList[i]->P = P;
		}
	}

    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static) 
		//�����������㣬ÿ�����ӵ�P����ȷ�����Ӷ����Լ�������ѹ�����ٶȡ��Լ����ӵı仯�ļ���
		for (int j = 0; j < fluidModel.particleList.size(); j++)
		{
			//�õ��ٽ��б�
			vector<WCSPHParticle*> neighbors = fluidModel.neighborList[j];
			//�õ��������ӵĸ��������ٽ��б�
			vector<StaticRigidParticle*> boundaryNeighborList = fluidModel.fluidBoundaryNeighborList[j];
			//�����ٽ��б����õ���ǰ���ӵ�ѹ�����ٶ�
			Vector3f pressureA = ComputePressureAcceleration(neighbors, boundaryNeighborList, fluidModel.particleList[j]);
			fluidModel.particleList[j]->pressure = pressureA;
			//��������ļ��㣬�õ���ǰ�������õ����ܼ��ٶ�
			Vector3f accelerationT = ComputeAcceleration(fluidModel.particleList[j]->pressure, fluidModel.particleList[j]->viscosity);
			fluidModel.particleList[j]->acceleration = accelerationT;
		}
	}
	
	//ˢ������ģ�͵������ӵļ��ٶ�
	fluidModel.UpdateParticlePosition();
	//ˢ���������ӵĹ�ϣֵ���������մ�������
	//���бʼǣ�������ﲻ�������е����ӣ���ô�������ӵ��ڽ�������һ��ʱ��ȫ����Ϊ0
	fluidModel.UpdateParticleHashValue();
	//���������ӵ�λ�ú��ٶȶ��������֮�����¶����ӽ������񻮷�
	MapParticleToGrid();
	//ˢ��֮�󣬸�����������������ھ��б�
	UpdateFluidNeighborList();
	//�����������ӵı߽������ھ��б�
	UpdateFluidBoundaryNeighborList();
	//���һ�μ���
}

//����������ص�ǰ���ӵ���Χ����
vector<WCSPHParticle*> WCSPHComputation::ComputeNeighborParticle(WCSPHParticle* origin)
{
	double h = SharedData::GetCellLength(); //�õ����ӵ�֧�Ű뾶��Ҳ��������Ԫ�ĳ���
	Vector3f originPosition = origin->position;
	//�õ�����������ӱ��
	Vector3f tempPosition; //���ڶ�27��������б����Ĺ���
	tempPosition.SetX(originPosition.GetX() - h);
	tempPosition.SetY(originPosition.GetY() - h);
	tempPosition.SetZ(originPosition.GetZ() - h);
	double originX = tempPosition.GetX();
	double originY = tempPosition.GetY();
	double originZ = tempPosition.GetZ();
	vector<WCSPHParticle*> neighborList; //���ӵ��ھ������б�
	int cellNum = hashGridList.GetCellNum(); //�õ���ǰ���������
	for (int i = 0; i < 3; i++, tempPosition.SetX(tempPosition.GetX() + h))
	{
		for (int j = 0; j < 3; j++, tempPosition.SetY(tempPosition.GetY() + h))
		{
			for (int k = 0; k < 3; k++, tempPosition.SetZ(tempPosition.GetZ() + h))
			{
				//���бʼǣ������Ҳ�����Ӧ�����񣬷��ص������Ǵ���ģ����ԣ�����ط��Ĵ����Ǵ����
				//���ԭ������һ���������ò����õ����������������ǹ�ϣ��������
				GridCell targetGridCell = hashGridList.GetGridCell(tempPosition);
				//����Ŀ�����ӵ�λ�ã�ֱ�Ӳ��ҵ���Ӧ�Ĺ�ϣ�б��е����񣬲��һ�ȡ���е������б�
				//���бʼǣ��ҵ���Ӧ����֮���޷��ҵ������е�����
				//�����ʽ������ӳ�����񿽱����󣬵����޷��ҵ������е�����
				vector<WCSPHParticle*> tempParticleList = targetGridCell.GetParticleList();
				//��Ŀ�������б��е�����һ��һ���������ھ��б�
				for (int q = 0; q < tempParticleList.size(); q++)
				{
					//����������������г����Լ�����ô�ͺ��Ե�
					if (origin->position == tempParticleList[q]->position)
					{
						continue;
					}
					//ֻѹ����h��Χ֮�ڣ������ӵ�Ӱ�췶Χ֮�ڣ����ھ�����
					if (FunctionKit::DistanceComputation(originPosition,
						tempParticleList[q]->position) < h)
					{
						neighborList.push_back(tempParticleList[q]);
					}
					//����������
				}
			}
			tempPosition.SetZ(originZ);
		}
		tempPosition.SetY(originY);
	}
	return neighborList;
}

//��ѯ����������Χ�ı߽�����
vector<StaticRigidParticle*> WCSPHComputation::ComputeNeighborBoundaryParticle(WCSPHParticle* origin)
{
	double h = SharedData::GetCellLength(); //�õ����ӵ�֧�Ű뾶��Ҳ��������Ԫ�ĳ���
	Vector3f originPosition = origin->position;
	//�õ�����������ӱ��
	Vector3f tempPosition; //���ڶ�27��������б����Ĺ���
	tempPosition.SetX(originPosition.GetX() - h);
	tempPosition.SetY(originPosition.GetY() - h);
	tempPosition.SetZ(originPosition.GetZ() - h);
	double originX = tempPosition.GetX();
	double originY = tempPosition.GetY();
	double originZ = tempPosition.GetZ();
	vector<StaticRigidParticle*> neighborList; //���ӵ��ھ������б�
	int cellNum = hashGridList.GetCellNum(); //�õ���ǰ���������

	for (int i = 0; i < 3; i++, tempPosition.SetX(tempPosition.GetX() + h))
	{
		for (int j = 0; j < 3; j++, tempPosition.SetY(tempPosition.GetY() + h))
		{
			for (int k = 0; k < 3; k++, tempPosition.SetZ(tempPosition.GetZ() + h))
			{
				//���бʼǣ������Ҳ�����Ӧ�����񣬷��ص������Ǵ���ģ����ԣ�����ط��Ĵ����Ǵ����
				//���ԭ������һ���������ò����õ����������������ǹ�ϣ��������
				GridCell targetGridCell = hashGridList.GetGridCell(tempPosition);
				//����Ŀ�����ӵ�λ�ã�ֱ�Ӳ��ҵ���Ӧ�Ĺ�ϣ�б��е����񣬲��һ�ȡ���е������б�
				//���бʼǣ��ҵ���Ӧ����֮���޷��ҵ������е�����
				//�����ʽ������ӳ�����񿽱����󣬵����޷��ҵ������е�����
				vector<StaticRigidParticle*> tempParticleList = targetGridCell.GetBoundaryParticleList();
				//��Ŀ�������б��е�����һ��һ���������ھ��б�
				for (int q = 0; q < tempParticleList.size(); q++)
				{
					//����������������г����Լ�����ô�ͺ��Ե�
					if (origin->position == tempParticleList[q]->position)
					{
						continue;
					}
					//ֻѹ����h��Χ֮�ڣ������ӵ�Ӱ�췶Χ֮�ڣ����ھ�����
					if (FunctionKit::DistanceComputation(originPosition,
						tempParticleList[q]->position) < h)
					{
						neighborList.push_back(tempParticleList[q]);
					}
					//����������
				}
			}
			tempPosition.SetZ(originZ);
		}
		tempPosition.SetY(originY);
	}
	return neighborList;
}

vector<StaticRigidParticle*> WCSPHComputation::ComputeNeighborBoundaryParticle(StaticRigidParticle* origin)
{
	double h = SharedData::GetCellLength(); //�õ����ӵ�֧�Ű뾶��Ҳ��������Ԫ�ĳ���
	Vector3f originPosition = origin->position;
	//�õ�����������ӱ��
	Vector3f tempPosition; //���ڶ�27��������б����Ĺ���
	tempPosition.SetX(originPosition.GetX() - h);
	tempPosition.SetY(originPosition.GetY() - h);
	tempPosition.SetZ(originPosition.GetZ() - h);
	double originX = tempPosition.GetX();
	double originY = tempPosition.GetY();
	double originZ = tempPosition.GetZ();
	vector<StaticRigidParticle*> neighborList; //���ӵ��ھ������б�
	int cellNum = hashGridList.GetCellNum(); //�õ���ǰ���������

	for (int i = 0; i < 3; i++, tempPosition.SetX(tempPosition.GetX() + h))
	{
		for (int j = 0; j < 3; j++, tempPosition.SetY(tempPosition.GetY() + h))
		{
			for (int k = 0; k < 3; k++, tempPosition.SetZ(tempPosition.GetZ() + h))
			{
				//���бʼǣ������Ҳ�����Ӧ�����񣬷��ص������Ǵ���ģ����ԣ�����ط��Ĵ����Ǵ����
				//���ԭ������һ���������ò����õ����������������ǹ�ϣ��������
				GridCell targetGridCell = hashGridList.GetGridCell(tempPosition);
				//����Ŀ�����ӵ�λ�ã�ֱ�Ӳ��ҵ���Ӧ�Ĺ�ϣ�б��е����񣬲��һ�ȡ���е������б�
				//���бʼǣ��ҵ���Ӧ����֮���޷��ҵ������е�����
				//�����ʽ������ӳ�����񿽱����󣬵����޷��ҵ������е�����
				vector<StaticRigidParticle*> tempParticleList = targetGridCell.GetBoundaryParticleList();
				//��Ŀ�������б��е�����һ��һ���������ھ��б�
				for (int q = 0; q < tempParticleList.size(); q++)
				{
					//����������������г����Լ�����ô�ͺ��Ե�
					if (origin->position == tempParticleList[q]->position)
					{
						continue;
					}
					//ֻѹ����h��Χ֮�ڣ������ӵ�Ӱ�췶Χ֮�ڣ����ھ�����
					if (FunctionKit::DistanceComputation(originPosition,
						tempParticleList[q]->position) < h)
					{
						neighborList.push_back(tempParticleList[q]);
					}
					//����������
				}
			}
			tempPosition.SetZ(originZ);
		}
		tempPosition.SetY(originY);
	}
	return neighborList;
}

double WCSPHComputation::ComputeDensity(vector<WCSPHParticle*>& neighborList,vector<StaticRigidParticle*>& boundaryNeighborList, WCSPHParticle* origin)
{
	int i;
	double mass = fluidModel.particleMass;
	double result = CubicKernel::W_zero()*mass;
	double h = hashGridList.GetCellLength(); //�õ���ǰ��Ԫ�ĳ���
	Vector3f xi = origin->position;
	//�����������Ӷ���������֮���Ӱ��

    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (i = 0; i < neighborList.size(); i++)
		{
			Vector3f xj = neighborList[i]->position;
			double distance = FunctionKit::DistanceComputation(origin->position, neighborList[i]->position);
			//�ں˺����ķ�Χ֮�ڵĻ�
			Vector3f temp123 = xi - xj;
			double temp = CubicKernel::W(xi - xj);
			result += CubicKernel::W(xi - xj)*mass;
			//result += SPHKernal::Poly6Kernel(h, distance);
		}
	}
	//����߽����Ӷ��������ӵ�Ӱ�죬�����ܶȲ���

    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int j = 0; j < boundaryNeighborList.size(); j++)
		{
			Vector3f xj = boundaryNeighborList[j]->position;
			double kernelResult = CubicKernel::W(xi - xj);
			double boundPSI = boundaryNeighborList[j]->boundaryPsi;
			//���õõ���Ӧ�ı߽��ھ�����
			result += boundaryNeighborList[j]->boundaryPsi * CubicKernel::W(xi - xj);
		}
	}
	return result;
}

//ע�⣺�ڸú�������ǰ������������������ӵ��ܶ��Լ�ѹǿ
Vector3f WCSPHComputation::ComputePressureAcceleration(vector<WCSPHParticle*>& neighborList, vector<StaticRigidParticle*>& boundaryNeighborList, WCSPHParticle* origin)
{
	//ע�� �ú���ʹ�õĻ�����SPH���ܶȼ��㷽ʽS
	int j;
	double m = fluidModel.particleMass;
	double h = fluidModel.particleSupportRad; //�õ�֧�Ű뾶
	double pi = origin->P; //�õ�ѹǿ
	double di = origin->density; //�õ���ǰ���ӵ��ܶ�
	double mass = fluidModel.particleMass; //�õ����ӵ�����
	Vector3f xi = origin->position; //�õ���ǰ���ӵ�λ��
	//��ʼ����
	double dpi = pi / (di*di);
	Vector3f result(0, 0, 0);

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static) 
		for (j = 0; j < neighborList.size(); j++)
		{
			double pj = neighborList[j]->P;
			double dj = neighborList[j]->density;
			double r = FunctionKit::DistanceComputation(neighborList[j]->position, origin->position);
			Vector3f xj = neighborList[j]->position;
			//����û��ʹ�ú˺�����ʹ�õ��ǻ�����SPH����
			//���бʼǣ��ɸú˺�����ֵΪ����
			//���������ʡ�Ե�������ļ���
			/*
			Vector3f vecA = rj - ri;
			double additionP = pi + pj;
			double multi = 2 * di*dj;
			double disSub = h - distance;
			Vector3f kl = vecA * (additionP / multi *disSub*disSub / distance);
			result = result + kl;
			*/
			double dpj = pj / (dj*dj);
			result -= mass * (dpi + dpj) *  CubicKernel::gradW(xi - xj);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
    //���ݱ߽����Ӽ�������Ӧ���ܵ�����
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static) 
		for (int j = 0; j < boundaryNeighborList.size(); j++)
		{
			StaticRigidParticle* boundaryNeighborParticle = boundaryNeighborList[j];
			Vector3f xj = boundaryNeighborParticle->position;
			Vector3f a = boundaryNeighborParticle->boundaryPsi * (dpi)*CubicKernel::gradW(xi - xj);
			result -= a;
		}
	}
	return result;
}

//����ճ����������ʹ�õ��ǻ���SPH�ļ��㷽ʽ
//ע���ڸú�������ǰ������������������ӵ��ٶ�
Vector3f WCSPHComputation::ComputeViscostyAcceleration(vector<WCSPHParticle*>& neighborList, vector<StaticRigidParticle*>& boundaryNeighborList, WCSPHParticle* origin)
{
	//�˴����õ���XSPH���е�Viscosity�ļ��㷽��
	Vector3f vi = origin->velocity;
	Vector3f xi = origin->position;
	double m = fluidModel.particleMass;
	double viscosityCons = viscosityConstant;
	double h = fluidModel.particleSupportRad;
	double di = origin->density;
	Vector3f result(0, 0, 0);
	double invH = 1.0 / SharedData::GetTimeStep();
	//////////////////////////////////////////////////////////////////////////
	///����������ȼ���
	//////////////////////////////////////////////////////////////////////////
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		for (int j = 0; j < neighborList.size(); j++)
		{
			Vector3f vj = neighborList[j]->velocity;
			Vector3f xj = neighborList[j]->position;
			double dj = neighborList[j]->density;
			double r = FunctionKit::DistanceComputation(xi, xj); //�õ�������ĳ���
			result -= invH*viscosityCons*(m / dj)*(vi - vj)*CubicKernel::W(xi - xj);
		}
	}
	//����߽�������ȵ�Ӱ��
	//////////////////////////////////////////////////////////////////////////
	///�߽�������ȵ�Ӱ��
	//////////////////////////////////////////////////////////////////////////
	//��������Ĵ��룬ģ��������ڻ�ɽ�ҽ���Ч���������Ⱥܴ�
  /*  #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		for (unsigned int j = 0; j < boundaryNeighborList.size(); j++)
		{
			const Vector3f xj = boundaryNeighborList[j].position;
			const Vector3f vj = boundaryNeighborList[j].velocity;
			result -= invH * viscosityCons * (boundaryNeighborList[j].boundaryPsi / di) * (vi)* CubicKernel::W(xi - xj);
		}
	}*/
	return result;
}

//����������ٶȣ��Ӷ����ȫ������
Vector3f WCSPHComputation::ComputeAcceleration(Vector3f pressureA, Vector3f viscosityA)
{
	//�����ļ��ٶ�,ע����Z��ı�ǣ��Լ����ֵ�����µ�
	Vector3f gA(0, 0, gravity); 
	Vector3f result = pressureA + viscosityA + gA;
	return result;
}

double WCSPHComputation::CorrectionParticleDensity(WCSPHParticle* a)
{
	return max(a->density, fluidModel.restDensity); //������ǰ�����ܶȣ����ڼ���ѹǿ
}

//����̩�ط��̼���õ����ӵ�ǰ��ѹǿ
double WCSPHComputation::CalculateParticleTaite(double a)
{
	return stiffnesss * (pow(a / fluidModel.restDensity, gamma) - 1.0);
}

void WCSPHComputation::UpdateFluidNeighborList()
{
	//��ˢ���µ������б��ʱ�������ǰ���ھ��б�
	fluidModel.ClearNeighborList();
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		//������ǰ�������ӵ��б�
		for (int i = 0; i < fluidModel.particleList.size(); i++)
		{
			//����õ����е��ٽ�������
			vector<WCSPHParticle*> temp = ComputeNeighborParticle(fluidModel.particleList[i]);
			for (int j = 0; j < temp.size(); j++)
			{
				//���ٽ�����������װ���Ӧ�������ھӱ���
				fluidModel.neighborList[i].push_back(temp[j]);
			}
		}
	}
}

//���������Ҫ����ÿһ�ζ���Ҫ�����
void WCSPHComputation::UpdateFluidBoundaryNeighborList()
{
	//����������Χ�ı߽�����
	fluidModel.ClearFluidBoundaryNeighborList();
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		for (int i = 0; i < fluidModel.particleList.size(); i++)
		{
			//�õ���Ӧ�����������е��ٽ���̬�����б�
			vector<StaticRigidParticle*> temp = ComputeNeighborBoundaryParticle(fluidModel.particleList[i]);
			for (int j = 0; j < temp.size(); j++)
			{
				//����Ӧ�ľ�̬��������װ�뵽������Χ�ľ�̬�����б���
				fluidModel.fluidBoundaryNeighborList[i].push_back(temp[j]);
			}
		}
	}
}

//�������ֻ��Ҫ�ڳ�ʼBoundaryPSI��ʼ����ʱ�����
void WCSPHComputation::ComputeBoundaryNeighborList()
{
	for (int i = 0; i < fluidModel.boundaryObj.staticRigidParticleList.size(); i++)
	{
		StaticRigidParticle* test2 = fluidModel.boundaryObj.staticRigidParticleList[i];
		//�õ���Ӧ�����������е��ٽ���̬�����б�
		vector<StaticRigidParticle*> temp = ComputeNeighborBoundaryParticle(fluidModel.boundaryObj.staticRigidParticleList[i]);
		for (int j = 0; j < temp.size(); j++)
		{
			//����Ӧ�ľ�̬��������װ�뵽��̬������Χ�ľ�̬�����б���
			fluidModel.boundaryNeighborList[i].push_back(temp[j]);
		}
	}	
}
