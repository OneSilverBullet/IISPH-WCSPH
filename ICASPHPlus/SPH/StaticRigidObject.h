#pragma once
#include "Vector3f.h"
#include "Matrix3f.h"

//////////////////////////////////////////////////////////////////////////
//StaticRigidObject��
//���ã����ڵ�����̬���������˵�����������������ɵľ�̬���忴��
//�������壬����൱�������˵�ǰ������̬�����������ԡ�
//////////////////////////////////////////////////////////////////////////
class StaticRigidObject
{
public:
	Vector3f m_postion;          //��ǰ��̬�������λ��
	Vector3f m_velocity;         //��ǰ��̬��������ٶȣ�����Ϊ0
	Vector3f m_angularVelocity;  //��̬����Ľ��ٶȣ�Ϊ0
	Matrix3f m_rotation;         //��ǰ��̬�����������ת�Ƕ�

public:
	//��ʼ������̬��������ٶ�Ϊ0
	StaticRigidObject() { m_velocity.Zero(); m_angularVelocity.Zero(); }
	double const GetMass()const { return 0.0; }
	Vector3f const & GetPosition() const { return m_postion; }
	Vector3f const &GetVelocity()const { return m_velocity; }
	Matrix3f const &GetRotation()const { return m_rotation; }
	Vector3f const &GetAngularVelocity() const { return m_angularVelocity; }
	void addForce(const Vector3f& f){}	                                 //����
	void addTorque(const Vector3f& f){}                                  //����
	void SetPosition(const Vector3f& x) { m_postion = x; }               //λ��
	void SetRotation(const Matrix3f& r) { m_rotation = r; }              //��ת
	//��Ҫ�Ĳ��������أ����ڰ�������Ķ����������и�ֵʱ�����Ե��øú���
	StaticRigidObject& operator=(StaticRigidObject& a)
	{
		m_postion = a.GetPosition();
		m_velocity = a.GetVelocity();
		m_angularVelocity = a.GetAngularVelocity();
		m_rotation = a.GetRotation();
		return *this;
	}
};