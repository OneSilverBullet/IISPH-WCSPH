#pragma once
#include "Vector3f.h"
#include "Matrix3f.h"

//////////////////////////////////////////////////////////////////////////
//StaticRigidObject类
//作用：用于单个静态物体的属性说明，这里把由粒子组成的静态物体看成
//单个物体，这个类当中描述了当前单个静态物体的相关属性。
//////////////////////////////////////////////////////////////////////////
class StaticRigidObject
{
public:
	Vector3f m_postion;          //当前静态刚性体的位置
	Vector3f m_velocity;         //当前静态刚性体的速度，设置为0
	Vector3f m_angularVelocity;  //静态物体的角速度，为0
	Matrix3f m_rotation;         //当前静态刚性物体的旋转角度

public:
	//初始化，静态刚性体的速度为0
	StaticRigidObject() { m_velocity.Zero(); m_angularVelocity.Zero(); }
	double const GetMass()const { return 0.0; }
	Vector3f const & GetPosition() const { return m_postion; }
	Vector3f const &GetVelocity()const { return m_velocity; }
	Matrix3f const &GetRotation()const { return m_rotation; }
	Vector3f const &GetAngularVelocity() const { return m_angularVelocity; }
	void addForce(const Vector3f& f){}	                                 //受力
	void addTorque(const Vector3f& f){}                                  //力矩
	void SetPosition(const Vector3f& x) { m_postion = x; }               //位置
	void SetRotation(const Matrix3f& r) { m_rotation = r; }              //旋转
	//必要的操作符重载，用于包含该类的对象的物体进行赋值时，隐性调用该函数
	StaticRigidObject& operator=(StaticRigidObject& a)
	{
		m_postion = a.GetPosition();
		m_velocity = a.GetVelocity();
		m_angularVelocity = a.GetAngularVelocity();
		m_rotation = a.GetRotation();
		return *this;
	}
};