#pragma once
#include "Vector3f.h"
#include "ShareData.h"
#include "FunctionKit.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////////
//Particle类
//作用：用于模拟流体的基本粒子，包含流体粒子的必须属性
//////////////////////////////////////////////////////////////////////////
class WCSPHParticle
{
public:
	Vector3f position;        //粒子的位置
	Vector3f velocity;        //粒子的速度
	Vector3f acceleration;    //粒子的加速度
	double P;                 //粒子的压强
	double density;           //粒子的密度
	Vector3f viscosity;       //粒子的黏度力
	Vector3f pressure;        //粒子的压力
	int hashIndex;            //当前粒子所属的网格哈希号

public:
	//构造函数直接调用初始化函数
	WCSPHParticle() {}
	~WCSPHParticle() {}
	void Initialization(Vector3f position, int v); //对粒子进行初始化更新

												   //对粒子进行赋值
	WCSPHParticle& operator=(WCSPHParticle& a);
};