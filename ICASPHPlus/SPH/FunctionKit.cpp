#include "FunctionKit.h"

double FunctionKit::DistanceComputation(Vector3f a, Vector3f b)
{
	double xTem = a.GetX() - b.GetX();
	double yTem = a.GetY() - b.GetY();
	double zTem = a.GetZ() - b.GetZ();
	double result = xTem*xTem + yTem*yTem + zTem*zTem;
	return double(sqrt(result));
}

//��ϣ����
int FunctionKit::GridHash(int i, int j, int k, int cellNum)
{
	//���бʼǣ������е�k=26��ʱ�򣬴�������ˣ�26*83492791�����˴���
	//���������ʹ����ǰȡģ�ķ�ʽ�Ӷ������˴������ĳ���
	long long int a = i*(P1%cellNum) ^ j*(P2%cellNum) ^ k*(P3%cellNum);
	int b = a % cellNum;
	return b;
}

//�����ӵ�λ��ӳ�䵽�����λ��
Vector3i FunctionKit::PositionMapIndex(Vector3f parPosition, double cellLength, double offset)
{
	int i = floor(float((parPosition.GetX() + offset) / cellLength));
	int j = floor(float((parPosition.GetY() + offset) / cellLength));
	int k = floor(float((parPosition.GetZ() + offset) / cellLength));
	Vector3i result(i, j, k);
	return result;
}

//ֱ�Ӵ����ӵ�λ��ֱ��ӳ�䵽��ϣ�����index
int FunctionKit::PositionMapHash(Vector3f parPosition, double cellLength, int cellNum, double offset)
{
	Vector3i particleCellIndex = PositionMapIndex(parPosition, cellLength, offset);
	int result = GridHash(particleCellIndex.GetX(),
		particleCellIndex.GetY(),
		particleCellIndex.GetZ(), cellNum);
	return result;
}

//����Ӧ�����Ƿ�Ϊ����
bool FunctionKit::CheckPrime(int n)
{
	int t = 1;
	if (n < 2){
		return false; //��ʱnΪ1�����䲻������
	}
	for (int i = 2; i < sqrt(n) + 1; i++)
	{
		if (n%i == 0){
			t = 0;
			break;
		}
	}
	if (t){
		return true;
	}
	return false;
}

double FunctionKit::GetVectorsCos(Vector3f a, Vector3f b)
{
	double result = a*b / (a.Length()*b.Length());
	return result;
}

double FunctionKit::PressureConstantB(double density0, double soundVelocity, double r)
{
	return double(density0*soundVelocity*soundVelocity / r);
}

double FunctionKit::TaitEquation(double density, double B ,int gamma, double restDensity)
{
	return double(B*(pow(density / restDensity, gamma) - 1));
}

int FunctionKit::CalculateGridCellNum(double radius)
{
	double d = 0.00001;
	//��ǰֻ��Ҫ���һ����λ����Χ����
	if (abs(radius - 0.025) < d)
	{
		return 1;
	}
	if ((radius > 0.025&&radius < 0.05) || (abs(radius - 0.05) < d))
	{
		return 2;
	}
	if ((radius > 0.05&&radius < 0.075) || (abs(radius - 0.075) < d))
	{
		return 3;
	}
	if ((radius > 0.075&&radius < 0.1) || (abs(radius - 0.1) < d))
	{
		return 4;
	}
	if ((radius > 0.1&&radius < 0.125) || (abs(radius - 0.125) < d))
	{
		return 5;
	}
	if ((radius > 0.125&&radius < 0.15) || (abs(radius - 0.15) < d))
	{
		return 6;
	}
	if ((radius > 0.15&&radius < 0.175) || (abs(radius - 0.175) < d))
	{
		return 7;
	}
	if ((radius > 0.175&&radius < 0.2) || (abs(radius - 0.2) < d))
	{
		return 8;
	}
	if ((radius > 0.2&&radius < 0.225) || (abs(radius - 0.225) < d))
	{
		return 9;
	}
	if ((radius > 0.2&&radius < 0.225) || (abs(radius - 0.25) < d))
	{
		return 10;
	}
}

Vector3f FunctionKit::Reflection(Vector3f normalize, Vector3f a)
{
	//��׼��
	Vector3f normal = normalize.Normalize(); 
	Vector3f aNormal = a.Normalize();
	double cos = GetVectorsCos(aNormal, -a);
	Vector3f tempVector = normal*a.Length()*cos; 
	Vector3f P = a + tempVector;  //�õ�ˮƽ����
	Vector3f result = P * 2 - a; 
	return result;
}
