#pragma once

#include "../ICASPHPlus/SPH/Vector3f.h"
#include "../ICASPHPlus/SPH/Matrix3f.h"
#include "Quaternion.h"

#define FORCE_INLINE __forceinline

class DynamicRigidObject
{
public:
	//质量
	double m_mass;   
	//逆质量
	double m_invMass;
	//质点的位置
	Vector3f m_x;
	Vector3f m_lastX;
	Vector3f m_oldX;
	Vector3f m_x0;
	//当前刚体的速度
	Vector3f m_v;
	//当前收到外力影响的加速度
	Vector3f m_a;

	/** Inertia tensor in the principal axis system: \n
	* After the main axis transformation the inertia tensor is a diagonal matrix.
	* So only three values are required to store the inertia tensor. These values
	* are constant over time.
	*/
	Vector3f m_inertiaTensor;
	/** Inverse inertia tensor in body space */
	Vector3f m_inertiaTensorInverse;
	/** 3x3 matrix, inverse of the inertia tensor in world space */
	//惯性张量在世界坐标系当中的转置矩阵
	Matrix3f m_inertiaTensorInverseW;
	/** Quaternion that describes the rotation of the body in world space */
	//Quaternion m_q;
	//Quaternion m_lastQ;
	//Quaternion m_oldQ;
	//Quaternion m_q0;
	/** Quaternion representing the rotation of the main axis transformation
	that is performed to get a diagonal inertia tensor */
	//Quaternion m_q_mat;
	/** Quaternion representing the initial rotation of the geometry */
	//Quaternion m_q_initial;
	/** difference of the initial translation and the translation of the main axis transformation */
	//Vector3f m_x0_mat;
	/** rotationMatrix = 3x3 matrix.
	* Important for the transformation from world in body space and vice versa.
	* When using quaternions the rotation matrix is computed out of the quaternion.
	*/
	Matrix3f m_rot;
	/** Angular velocity, defines rotation axis and velocity (magnitude of the vector) */
	//角速度，定义旋转轴以及其速度
	Vector3f m_omega;
	/** external torque */
	Vector3f m_torque;

	//double m_restitutionCoeff;
	//double m_frictionCoeff;

	// RigidBodyGeometry m_geometry;

	// transformation required to transform a point to local space or vice vera
	//Matrix3f m_transformation_R;
	//Vector3f m_transformation_v1;
	//Vector3f m_transformation_v2;
	//Vector3f m_transformation_R_X_v1;

public:
	DynamicRigidObject(void)
	{
	}

	~DynamicRigidObject(void)
	{
	}

	//void initBody(const double mass, const Vector3f &x,
	//	const Vector3f &inertiaTensor, const Quaternion &rotation,
	//	//const VertexData &vertices, const IndexedFaceMesh &mesh,
	//	const Vector3f &scale = Vector3f(1.0, 1.0, 1.0))
	//{
	//	setMass(mass);
	//	m_x = x;
	//	m_x0 = x;
	//	m_lastX = x;
	//	m_oldX = x;
	//	m_v.Zero();
	//	m_a.Zero();

	//	setInertiaTensor(inertiaTensor);
	//	m_q = rotation;
	//	m_q0 = rotation;
	//	m_lastQ = rotation;
	//	m_oldQ = rotation;
	//	m_rot = m_q.matrix();
	//	//m_q_mat = Quaternionr(1.0, 0.0, 0.0, 0.0);
	//	//m_q_initial = Quaternionr(1.0, 0.0, 0.0, 0.0);
	//	m_x0_mat.Zero();
	//	//rotationUpdated();
	//	m_omega.Zero();
	//	m_torque.Zero();

	//	m_restitutionCoeff = 0.6;
	//	m_frictionCoeff = 0.2;

	//	//getGeometry().initMesh(vertices.size(), mesh.numFaces(), &vertices.getPosition(0), mesh.getFaces().data(), mesh.getUVIndices(), mesh.getUVs(), scale);
	//	//getGeometry().updateMeshTransformation(getPosition(), getRotationMatrix());
	//}

	/*void initBody(const Real density, const Vector3r &x, const Quaternionr &rotation,
		const VertexData &vertices, const IndexedFaceMesh &mesh, const Vector3r &scale = Vector3r(1.0, 1.0, 1.0))
	{
		m_mass = 1.0;
		m_inertiaTensor = Vector3r(1.0, 1.0, 1.0);
		m_x = x;
		m_x0 = x;
		m_lastX = x;
		m_oldX = x;
		m_v.setZero();
		m_a.setZero();

		m_q = rotation;
		m_q0 = rotation;
		m_lastQ = rotation;
		m_oldQ = rotation;
		m_rot = m_q.matrix();
		rotationUpdated();
		m_omega.setZero();
		m_torque.setZero();

		m_restitutionCoeff = 0.6;
		m_frictionCoeff = 0.2;

		getGeometry().initMesh(vertices.size(), mesh.numFaces(), &vertices.getPosition(0), mesh.getFaces().data(), mesh.getUVIndices(), mesh.getUVs(), scale);
		determineMassProperties(density);
		getGeometry().updateMeshTransformation(getPosition(), getRotationMatrix());
	}
*/
	void reset()
	{
		getPosition() = getPosition0();
		getOldPosition() = getPosition0();
		getLastPosition() = getPosition0();

		//getRotation() = getRotation0();
		//getOldRotation() = getRotation0();
		//getLastRotation() = getRotation0();

		getVelocity().Zero();
		getAngularVelocity().Zero();

		getAcceleration().Zero();
		getTorque().Zero();

		//rotationUpdated();
	}

	void updateInverseTransformation()
	{
		 //remove the rotation of the main axis transformation that is performed
		 //to get a diagonal inertia tensor since the distance function is 
		 //evaluated in local coordinates
		
		 //transformation world to local:
		 //p_local = R_initial^T ( R_MAT R^T (p_world - x) - x_initial + x_MAT)
		 
		 //transformation local to world:
		 //p_world = R R_MAT^T (R_initial p_local + x_initial - x_MAT) + x
		
	/*	m_transformation_R = (getRotationInitial().inverse() * getRotationMAT() * getRotation().inverse()).matrix();
		m_transformation_v1 = -getRotationInitial().inverse().matrix() * getPositionInitial_MAT();*/
		//m_transformation_v2 = (getRotation()*getRotationMAT().inverse()).matrix() * getPositionInitial_MAT() + getPosition();
		//m_transformation_R_X_v1 = -m_transformation_R * getPosition() + m_transformation_v1;
	}

	/*void rotationUpdated()
	{
		if (m_mass != 0.0)
		{
			m_rot = m_q.matrix();
			updateInverseInertiaW();
			updateInverseTransformation();
		}
	}*/

	/*void updateInverseInertiaW()
	{
		if (m_mass != 0.0)
		{
			m_inertiaTensorInverseW = m_rot * m_inertiaTensorInverse.asDiagonal() * m_rot.transpose();
		}
	}*/

	///** Determine mass and inertia tensor of the given geometry.
	//*/
	//void determineMassProperties(const Real density)
	//{
	//	// apply initial rotation
	//	VertexData &vd = m_geometry.getVertexDataLocal();
	//	for (unsigned int i = 0; i < vd.size(); i++)
	//		vd.getPosition(i) = m_rot * vd.getPosition(i) + m_x0;

	//	VolumeIntegration vi(m_geometry.getMesh(), m_geometry.getVertexDataLocal());
	//	vi.compute_inertia_tensor(density);

	//	// Diagonalize Inertia Tensor
	//	Eigen::SelfAdjointEigenSolver<Matrix3r> es(vi.getInertia());
	//	Vector3r inertiaTensor = es.eigenvalues();
	//	Matrix3r R = es.eigenvectors();

	//	setMass(vi.getMass());
	//	setInertiaTensor(inertiaTensor);

	//	if (R.determinant() < 0.0)
	//		R = -R;

	//	// rotate vertices back				
	//	for (unsigned int i = 0; i < vd.size(); i++)
	//		vd.getPosition(i) = R.transpose() * (vd.getPosition(i) - vi.getCenterOfMass());

	//	// set rotation
	//	Quaternionr qR = Quaternionr(R);
	//	qR.normalize();
	//	m_q_mat = qR;
	//	m_q_initial = m_q0;
	//	m_x0_mat = m_x0 - vi.getCenterOfMass();

	//	m_q0 = qR;
	//	m_q = m_q0;
	//	m_lastQ = m_q0;
	//	m_oldQ = m_q0;
	//	rotationUpdated();

	//	// set translation
	//	m_x0 = vi.getCenterOfMass();
	//	m_x = m_x0;
	//	m_lastX = m_x0;
	//	m_oldX = m_x0;
	//	updateInverseTransformation();
	//}

	//const Matrix3r &getTransformationR() { return m_transformation_R; }
	//const Vector3r &getTransformationV1() { return m_transformation_v1; }
	//const Vector3r &getTransformationV2() { return m_transformation_v2; }
	//const Vector3r &getTransformationRXV1() { return m_transformation_R_X_v1; }

	FORCE_INLINE double &getMass()
	{
		return m_mass;
	}

	FORCE_INLINE const double &getMass() const
	{
		return m_mass;
	}

	FORCE_INLINE void setMass(const double &value)
	{
		m_mass = value;
		if (m_mass != 0.0)
			m_invMass = 1.0 / m_mass;
		else
			m_invMass = 0.0;
	}

	FORCE_INLINE const double &getInvMass() const
	{
		return m_invMass;
	}

	FORCE_INLINE Vector3f &getPosition()
	{
		return m_x;
	}

	FORCE_INLINE const Vector3f &getPosition() const
	{
		return m_x;
	}

	FORCE_INLINE void setPosition(const Vector3f &pos)
	{
		m_x = pos;
	}

	FORCE_INLINE Vector3f &getLastPosition()
	{
		return m_lastX;
	}

	FORCE_INLINE const Vector3f &getLastPosition() const
	{
		return m_lastX;
	}

	FORCE_INLINE void setLastPosition(const Vector3f &pos)
	{
		m_lastX = pos;
	}

	FORCE_INLINE Vector3f &getOldPosition()
	{
		return m_oldX;
	}

	FORCE_INLINE const Vector3f &getOldPosition() const
	{
		return m_oldX;
	}

	FORCE_INLINE void setOldPosition(const Vector3f &pos)
	{
		m_oldX = pos;
	}

	FORCE_INLINE Vector3f &getPosition0()
	{
		return m_x0;
	}

	FORCE_INLINE const Vector3f &getPosition0() const
	{
		return m_x0;
	}

	FORCE_INLINE void setPosition0(const Vector3f &pos)
	{
		m_x0 = pos;
	}

	/*FORCE_INLINE Vector3f &getPositionInitial_MAT()
	{
		return m_x0_mat;
	}

	FORCE_INLINE const Vector3f &getPositionInitial_MAT() const
	{
		return m_x0_mat;
	}

	FORCE_INLINE void setPositionInitial_MAT(const Vector3f &pos)
	{
		m_x0_mat = pos;
	}
*/
	FORCE_INLINE Vector3f &getVelocity()
	{
		return m_v;
	}

	FORCE_INLINE const Vector3f &getVelocity() const
	{
		return m_v;
	}

	FORCE_INLINE void setVelocity(const Vector3f &value)
	{
		m_v = value;
	}

	FORCE_INLINE Vector3f &getAcceleration()
	{
		return m_a;
	}

	FORCE_INLINE const Vector3f &getAcceleration() const
	{
		return m_a;
	}

	FORCE_INLINE void setAcceleration(const Vector3f &accel)
	{
		m_a = accel;
	}

	FORCE_INLINE const Vector3f &getInertiaTensor() const
	{
		return m_inertiaTensor;
	}

	FORCE_INLINE void setInertiaTensor(const Vector3f &value)
	{
		m_inertiaTensor = value;
		m_inertiaTensorInverse = Vector3f(1.0 / value.x, 1.0 / value.y, 1.0 / value.z);
	}

	FORCE_INLINE const Vector3f &getInertiaTensorInverse() const
	{
		return m_inertiaTensorInverse;
	}

	FORCE_INLINE Matrix3f &getInertiaTensorInverseW()
	{
		return m_inertiaTensorInverseW;
	}

	FORCE_INLINE const Matrix3f &getInertiaTensorInverseW() const
	{
		return m_inertiaTensorInverseW;
	}

	FORCE_INLINE void setInertiaTensorInverseW(const Matrix3f &value)
	{
		m_inertiaTensorInverseW = value;
	}

	//FORCE_INLINE Quaternionr &getRotation()
	//{
	//	return m_q;
	//}

	//FORCE_INLINE const Quaternionr &getRotation() const
	//{
	//	return m_q;
	//}

	//FORCE_INLINE void setRotation(const Quaternionr &value)
	//{
	//	m_q = value;
	//}

	//FORCE_INLINE Quaternionr &getLastRotation()
	//{
	//	return m_lastQ;
	//}

	//FORCE_INLINE const Quaternionr &getLastRotation() const
	//{
	//	return m_lastQ;
	//}

	//FORCE_INLINE void setLastRotation(const Quaternionr &value)
	//{
	//	m_lastQ = value;
	//}

	//FORCE_INLINE Quaternionr &getOldRotation()
	//{
	//	return m_oldQ;
	//}

	//FORCE_INLINE const Quaternionr &getOldRotation() const
	//{
	//	return m_oldQ;
	//}

	//FORCE_INLINE void setOldRotation(const Quaternionr &value)
	//{
	//	m_oldQ = value;
	//}

	//FORCE_INLINE Quaternionr &getRotation0()
	//{
	//	return m_q0;
	//}

	//FORCE_INLINE const Quaternionr &getRotation0() const
	//{
	//	return m_q0;
	//}

	//FORCE_INLINE void setRotation0(const Quaternionr &value)
	//{
	//	m_q0 = value;
	//}

	//FORCE_INLINE Quaternionr &getRotationMAT()
	//{
	//	return m_q_mat;
	//}

	//FORCE_INLINE const Quaternionr &getRotationMAT() const
	//{
	//	return m_q_mat;
	//}

	//FORCE_INLINE void setRotationMAT(const Quaternionr &value)
	//{
	//	m_q_mat = value;
	//}

	//FORCE_INLINE Quaternionr &getRotationInitial()
	//{
	//	return m_q_initial;
	//}

	//FORCE_INLINE const Quaternionr &getRotationInitial() const
	//{
	//	return m_q_initial;
	//}

	//FORCE_INLINE void setRotationInitial(const Quaternionr &value)
	//{
	//	m_q_initial = value;
	//}

	FORCE_INLINE Matrix3f &getRotationMatrix()
	{
		return m_rot;
	}

	FORCE_INLINE const Matrix3f &getRotationMatrix() const
	{
		return m_rot;
	}

	FORCE_INLINE void setRotationMatrix(const Matrix3f &value)
	{
		m_rot = value;
	}

	FORCE_INLINE Vector3f &getAngularVelocity()
	{
		return m_omega;
	}

	FORCE_INLINE const Vector3f &getAngularVelocity() const
	{
		return m_omega;
	}

	FORCE_INLINE void setAngularVelocity(const Vector3f &value)
	{
		m_omega = value;
	}

	FORCE_INLINE Vector3f &getTorque()
	{
		return m_torque;
	}

	FORCE_INLINE const Vector3f &getTorque() const
	{
		return m_torque;
	}

	FORCE_INLINE void setTorque(const Vector3f &value)
	{
		m_torque = value;
	}

	//FORCE_INLINE double getRestitutionCoeff() const
	//{
	//	return m_restitutionCoeff;
	//}

	//FORCE_INLINE void setRestitutionCoeff(double val)
	//{
	//	m_restitutionCoeff = val;
	//}

	//FORCE_INLINE double getFrictionCoeff() const
	//{
	//	return m_frictionCoeff;
	//}

	//FORCE_INLINE void setFrictionCoeff(double val)
	//{
	//	m_frictionCoeff = val;
	//}

	/*RigidBodyGeometry& getGeometry()
	{
		return m_geometry;
	}*/
};