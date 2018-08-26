#pragma once
//����ű����м�¼һЩ���������
#include <math.h>
#include "Vector3f.h"
#define PI acos(-1)
//�������ӵĳ�ʼֵ
#define DENSITY0 1000.0           //����涨��ʼ�ܶ�Ϊ1000
#define MASS0    1.0             //���е����Ӷ�����Ϊ10��λ������
//ע�⣺�������������ˮ�嵱�о����ж��ٸ����ӣ�һ���ռ䣬����Ӧ��
//��ϣ�����
#define P1 73856093
#define P2 19349663
#define P3 83492791
#define MAX_DIS 65536
//�������������֡��
#define FRAME_PER_SECOND 25.0   //ÿ������ܶ���֡
//Ϊ�˺�������һЩoffset


const Vector3f POS_X_NORM(-1, 0, 0);
const Vector3f NEG_X_NORM(1, 0, 0);
const Vector3f POS_Y_NORM(0, -1, 0);
const Vector3f NEG_Y_NORM(0, 1, 0);
const Vector3f POS_Z_NORM(0, 0, -1);
const Vector3f NEG_Z_NORM(0, 0, 1);

enum PARTICLE_TYPE { S, s, o, l, L };


//////////////////////////////////////////////////////////////////////////
//SharedData��
//���ã������Ϊȫ�ֳ����ṩ��һ��ȫ�ֽӿڣ��Ӷ�ΪSPH�ĳ����ṩ��
//ȫ������͸�������ԡ�����ͨ�����������Խ����޸ġ�
//////////////////////////////////////////////////////////////////////////
class SharedData
{
public:
	static Vector3f boundaryOrigin;      //�߽����ʼ��
	static Vector3f boundaryWLH;         //�߽�ĳ����
	static double offset;                //�߽��ƫ��
	static double timeStep;              //WCSPH�����ʱ�䲽��
	static double particleRad;           //���ӵ����Ӱ뾶
	static double particleSupportRad;    //���ӵ�֧�Ű뾶
	static double particleMass;          //���ӵ�����
	static double restDensity;           //���ӵľ�ֹ�ܶ�
	static int    particleNum;           //���ӵ�����
	static Vector3f fluidWLH;            //�����WLH�ı߽糤��
	static Vector3f fluidOrigin;         //�������ɵ���ʼ��
	static double cellLength;            //����ĳ������ã������ӵ�֧�Ű뾶ϢϢ���
	static int hashCellNum;              //��ϣ��������ȷ���ĸ���
	static double stiffnesss;            //���嵱�еĸ���
	static int gamma;                    //���嵱�е�gammaֵ
	static double viscosityConstant;     //����������Ҫ�ĳ���
	static double gravity;               //���������

public:
	static void SetBoundaryOrigin(const Vector3f& a);
	static Vector3f GetBoundaryOrigin();
	static void SetBoundaryWLH(const Vector3f& a);
	static Vector3f GetBoundaryWLH();
	static void SetOffset(const double a);
	static double GetOffset();
	static void SetTimeStep(const double a);
	static double GetTimeStep();
	static void SetParticleRad(const double a);
	static double GetParticleRad();
	static void SetParticleSupportRad(const double a);
	static double GetParticleSupportRad();
	static void SetParticleMass(const double a);
	static double GetParticleMass();
	static void SetRestDensity(const double a);
	static double GetRestDensity();
	static void SetParticleNum(const int a);
	static int GetParticleNum();
	static void SetFluidWLH(const Vector3f& a);
	static Vector3f GetFluidWLH();
	static void SetFluidOrigin(const Vector3f& a);
	static Vector3f GetFluidOrigin();
	static void SetCellLength(const double a);
	static double GetCellLength();
	static void SetHashCellNum(const int a);
	static int GetHashCellNum();
	static void SetStiffnesss(const double a);
	static double GetStiffnesss();
	static void SetGamma(const int a);
	static int GetGamma();
	static void SetViscosityConstant(const double a);
	static double GetViscosityConstant();
	static void SetGravity(const double a);
	static double GetGravity();

};