#pragma once
#include <math.h>
#include "Vector3f.h"
#include "Vector3i.h"
#include "ShareData.h"


//////////////////////////////////////////////////////////////////////////
//FunctionKit��
//���ã�����һЩɢ�ҵĹ��ߺ���������������ȵ�
//////////////////////////////////////////////////////////////////////////
class FunctionKit
{
public:
	//�����������ĺ���
	static double DistanceComputation(Vector3f a, Vector3f b);
	//����Ĺ�ϣ���㺯��
	static int GridHash(int i, int j, int k, int cellNum); 
	//ͨ����ǰλ�ã�����ӳ�䵱ǰ��������
	static Vector3i PositionMapIndex(Vector3f parPosition, double cellLength, double offset);
	//ͨ������λ��ӳ�䵱ǰ��������
	static int PositionMapHash(Vector3f parPosition, double cellLength, int cellNum, double offset); 
	//����Ӧ�������Ƿ�Ϊ����
	static bool CheckPrime(int n);
	//������������
	//��������������cosֵ
	static double GetVectorsCos(Vector3f a, Vector3f b);
	//���ݷ����������ض�Ӧ�ķ�������
	static Vector3f Reflection(Vector3f normal, Vector3f a);
	//�����������ڼ���״̬���̵�ʱ��ʹ��
	static double PressureConstantB(double density0, double soundVelocity, double r);
	//̩�ط��̼���
	static double TaitEquation(double density, double B,int gamma, double restDensity);
	//�������ӵİ뾶������õ���ǰ������Ҫ���ٸ�����
	static int CalculateGridCellNum(double radius);

};