#pragma once
//这个脚本当中记录一些共享的数据
#include <math.h>
#include "Vector3f.h"
#define PI acos(-1)
//关于粒子的初始值
#define DENSITY0 1000.0           //这里规定初始密度为1000
#define MASS0    1.0             //所有的粒子都设置为10单位的质量
//注意：质量决定了这个水体当中究竟有多少个粒子，一定空间，自适应的
//哈希表参数
#define P1 73856093
#define P2 19349663
#define P3 83492791
#define MAX_DIS 65536
//关于整个程序的帧数
#define FRAME_PER_SECOND 25.0   //每秒程序跑多少帧
//为核函数增大一些offset


const Vector3f POS_X_NORM(-1, 0, 0);
const Vector3f NEG_X_NORM(1, 0, 0);
const Vector3f POS_Y_NORM(0, -1, 0);
const Vector3f NEG_Y_NORM(0, 1, 0);
const Vector3f POS_Z_NORM(0, 0, -1);
const Vector3f NEG_Z_NORM(0, 0, 1);

enum PARTICLE_TYPE { S, s, o, l, L };


//////////////////////////////////////////////////////////////////////////
//SharedData类
//作用：这个类为全局程序提供了一个全局接口，从而为SPH的程序提供了
//全局属性透明的特性。不可通过这个类对属性进行修改。
//////////////////////////////////////////////////////////////////////////
class SharedData
{
public:
	static Vector3f boundaryOrigin;      //边界的起始点
	static Vector3f boundaryWLH;         //边界的长宽高
	static double offset;                //边界的偏移
	static double timeStep;              //WCSPH计算的时间步骤
	static double particleRad;           //粒子的粒子半径
	static double particleSupportRad;    //粒子的支撑半径
	static double particleMass;          //粒子的质量
	static double restDensity;           //粒子的静止密度
	static int    particleNum;           //粒子的数量
	static Vector3f fluidWLH;            //流体的WLH的边界长度
	static Vector3f fluidOrigin;         //流体生成的起始点
	static double cellLength;            //网格的长度设置，和粒子的支撑半径息息相关
	static int hashCellNum;              //哈希网格最终确定的个数
	static double stiffnesss;            //流体当中的刚性
	static int gamma;                    //流体当中的gamma值
	static double viscosityConstant;     //计算黏度所需要的常量
	static double gravity;               //流体的重力

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