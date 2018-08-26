#include "DynamicRigidParticleObject.h"

void DynamicRigidParticleObject::Initialize()
{
	InitParticle();
	InitRigidBody();
}

void DynamicRigidParticleObject::InitParticle()
{
	//从哈希网格列表当中得到两个关键的静态值,调用的都是静态值
	int cellNum = SharedData::GetHashCellNum();
	double cellLength = SharedData::GetCellLength();
	//直径
	double diam = radius * 2;
	//设置好当前x方向的粒子个数
	int xLength = floor(length / diam + 0.5f) + 1;
	//设置好当前动态块的长宽高上粒子的数量
	dynamicCubeWLH.SetValue(xLength, xLength, xLength);

	std::cout << "初始化动态刚体" << std::endl;
	std::vector<Vector3f> boundaryParticles;


	Vector3f boundaryMinP(ox, oy, oz);
	Vector3f boundaryMaxP(ox + length, oy + length, oz + length);

	//使用直径来作为每一步向前的变化量
	double xshift = diam;
	double yshift = diam;
	double zshift = diam;

	const int stepsX = xLength;
	const int stepsY = xLength;
	const int stepsZ = xLength;

	Vector3f start = boundaryMinP;

	//六面墙，使用粒子进行墙的创建
	for (int i = 0; i < stepsX; i++)
	{
		for (int j = 0; j < stepsY; j++)
		{
			Vector3f currPos = Vector3f(i*xshift, j*yshift, 0) + boundaryMinP;
			//根据当前粒子的位置求得当前粒子的哈希值
			int hashValue = FunctionKit::PositionMapHash(currPos, cellLength, cellNum, boundaryOffset);
			//创建一个新的粒子
			DynamicRigidParticle* newParticle = new DynamicRigidParticle(currPos, hashValue);
			//将当前粒子装入到粒子列表当中
			dynamicParticles.push_back(newParticle);
		}
	}
	for (int i = 0; i < stepsX; i++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = Vector3f(i*xshift, 0, k*zshift) + boundaryMinP;
			//根据当前粒子的位置求得当前粒子的哈希值
			int hashValue = FunctionKit::PositionMapHash(currPos, cellLength, cellNum, boundaryOffset);
			//创建一个新的粒子
			DynamicRigidParticle* newParticle = new DynamicRigidParticle(currPos, hashValue);
			//将当前粒子装入到粒子列表当中
			dynamicParticles.push_back(newParticle);
		}
	}
	for (int j = 0; j < stepsY; j++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = Vector3f(0, j*yshift, k*zshift) + boundaryMinP;
			//根据当前粒子的位置求得当前粒子的哈希值
			int hashValue = FunctionKit::PositionMapHash(currPos, cellLength, cellNum, boundaryOffset);
			//创建一个新的粒子
			DynamicRigidParticle* newParticle = new DynamicRigidParticle(currPos, hashValue);
			//将当前粒子装入到粒子列表当中
			dynamicParticles.push_back(newParticle);
		}
	}
	for (int i = 0; i < stepsX; i++)
	{
		for (int j = 0; j < stepsY; j++)
		{
			Vector3f currPos = -Vector3f(i*xshift, j*yshift, 0) + boundaryMaxP;
			//根据当前粒子的位置求得当前粒子的哈希值
			int hashValue = FunctionKit::PositionMapHash(currPos, cellLength, cellNum, boundaryOffset);
			//创建一个新的粒子
			DynamicRigidParticle* newParticle = new DynamicRigidParticle(currPos, hashValue);
			//将当前粒子装入到粒子列表当中
			dynamicParticles.push_back(newParticle);
		}
	}
	for (int i = 0; i < stepsX; i++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = -Vector3f(i*xshift, 0, k*zshift) + boundaryMaxP;
			//根据当前粒子的位置求得当前粒子的哈希值
			int hashValue = FunctionKit::PositionMapHash(currPos, cellLength, cellNum, boundaryOffset);
			//创建一个新的粒子
			DynamicRigidParticle* newParticle = new DynamicRigidParticle(currPos, hashValue);
			//将当前粒子装入到粒子列表当中
			dynamicParticles.push_back(newParticle);
		}
	}
	for (int j = 0; j < stepsY; j++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = -Vector3f(0, j*yshift, k*zshift) + boundaryMaxP;
			//根据当前粒子的位置求得当前粒子的哈希值
			int hashValue = FunctionKit::PositionMapHash(currPos, cellLength, cellNum, boundaryOffset);
			//创建一个新的粒子
			DynamicRigidParticle* newParticle = new DynamicRigidParticle(currPos, hashValue);
			//将当前粒子装入到粒子列表当中
			dynamicParticles.push_back(newParticle);
		}
	}



}

void DynamicRigidParticleObject::InitRigidBody()
{
	//首先更新当前刚体的位置，设置为当前的质心
	Vector3f position;
	position.Zero();
	for (int i=0; i<dynamicParticles.size(); i++)
	{
		position += dynamicParticles[i]->m_x;
	}
	position /= dynamicParticles.size();
	m_rigidBody.m_x = position;
	m_rigidBody.m_lastX = position;
	m_rigidBody.m_oldX = position;
	m_rigidBody.m_x0 = position;
	//计算所有的粒子到刚体中心的距离。
	CalculateDiffDistance();

	//更新当前刚体的质量，设置为所有粒子质量的总和
	m_rigidBody.getMass() = dynamicParticles.size()*particleMass;
	//计算当前的惯性张量，这个值为定值
	CalculateInertiaTensor();
}

//初始化用
void DynamicRigidParticleObject::CalculateInertiaTensor()
{
	//初始化刚体粒子之后，根据刚体粒子对惯性张量进行计算，这个值是一个定值
	//由于其生成的物体为一个对称物体，除了矩阵对角线方向以外其他元素都是
	//0，因此只计算是三个量，并且对刚体的存储性质以Vec3形式
	double IXX = 0;
	double IYY = 0;
	double IZZ = 0;
	for (int i = 0; i < dynamicParticles.size(); i++)
	{
		Vector3f relativePosition;
		relativePosition = dynamicParticles[i]->m_x - m_rigidBody.getPosition();
		IXX += (relativePosition.y*relativePosition.y + relativePosition.z*relativePosition.z)*particleMass;
		IYY += (relativePosition.x*relativePosition.x + relativePosition.z*relativePosition.z)*particleMass;
		IZZ += (relativePosition.x*relativePosition.x + relativePosition.y*relativePosition.y)*particleMass;
	}
	//这两个值基本是不变的
	m_rigidBody.m_inertiaTensor.SetValue(IXX, IYY, IZZ);
	m_rigidBody.m_inertiaTensorInverse.SetValue(double(1.0 / IXX), double(1.0 / IYY), double(1.0 / IZZ));
	Matrix3f diagonalInertiaTensor;
	//这里搞一手惯性张量在世界坐标系下的逆矩阵，随时间变化
	m_rigidBody.m_inertiaTensorInverseW = m_rigidBody.m_rot * diagonalInertiaTensor.Diagonal(m_rigidBody.m_inertiaTensorInverse) *m_rigidBody.m_rot.Tranpose();
}

//实时调用更新
void DynamicRigidParticleObject::CalculateRotationMat()
{
	//根据角速度求得对应刚体的方向矩阵
	Vector3f omiga = m_rigidBody.getAngularVelocity();
	Matrix3f temp;
	Matrix3f rotationMatrix = temp.Rotate(omiga.x, omiga.y, omiga.z);
	m_rigidBody.setRotationMatrix(rotationMatrix);
}

void DynamicRigidParticleObject::CalculateDiffDistance()
{
	for (int j = 0; j < (int)dynamicParticles.size(); j++)
	{
		dynamicParticles[j]->diffDistance = dynamicParticles[j]->m_x - m_rigidBody.getPosition();
	}
}

void DynamicRigidParticleObject::UpdateObjectPar()
{
	if (isDynamic())
	{
		Vector3f force, torque;
		force.Zero();
		torque.Zero();
		//遍历当前所有粒子的位置
		for (int j = 0; j < (int)dynamicParticles.size(); j++)
		{
			force += dynamicParticles[j]->m_f;
			torque += (dynamicParticles[j]->m_x - m_rigidBody.getPosition()).Cross(dynamicParticles[j]->m_f);
			//在计算结束后，将当前粒子的力设置为0
			dynamicParticles[j]->m_f.Zero();
		}
		//刷新得到当前刚体的速度
		addForce(force);
		//刷新得到当前刚体的角速度
		addTorque(torque);
		//依据当前角速度更新当前刚体的旋转矩阵
		CalculateRotationMat();
		if(m_rigidBody.getPosition().y <0.2)
		{
			Vector3f newPos(m_rigidBody.getPosition().x, m_rigidBody.getPosition().z, 0.2);
			m_rigidBody.setPosition(newPos);
		}
		else
		{
			//更新当前刚体的位置，依据当前速度对位置进行刷新
			m_rigidBody.setPosition(m_rigidBody.getPosition() + m_rigidBody.getVelocity()*m_h);
		}
		//依据当前的的旋转矩阵，更新当前世界坐标下的惯性张量逆矩阵
		Matrix3f diagonalInertiaTensor;
		m_rigidBody.m_inertiaTensorInverseW = m_rigidBody.m_rot * diagonalInertiaTensor.Diagonal(m_rigidBody.m_inertiaTensorInverse) *m_rigidBody.m_rot.Tranpose();
	}
}

void DynamicRigidParticleObject::UpdateParticlePar()
{
	if (isDynamic() )
	{
		int cellNum = SharedData::GetHashCellNum();
		double cellLength = SharedData::GetCellLength();
        #pragma omp parallel default(shared)
		{
            #pragma omp for schedule(static)  
			for (int j = 0; j < (int)dynamicParticles.size(); j++)
			{
				if (j == 0)
				{
					std::cout << dynamicParticles[j]->m_x0 << std::endl;
				}
				Matrix3f test1 = m_rigidBody.m_rot;
				Vector3f test2 = test1 * dynamicParticles[j]->m_x0;
				//依据当前刚体的位置得到当前每个粒子的位置
				dynamicParticles[j]->m_x = m_rigidBody.m_rot * dynamicParticles[j]->diffDistance + getPosition(); //设置为上一帧的位置
				dynamicParticles[j]->m_x0 = dynamicParticles[j]->m_x;
				//更新当前刚体粒子的哈希值
				int hashValue = FunctionKit::PositionMapHash(dynamicParticles[j]->m_x, cellLength,
					cellNum, boundaryOffset);
				//在计算结束后，将当前粒子的力设置为0
				dynamicParticles[j]->hashValue = hashValue;

				//依据当前刚体的位置和角速度得到当前粒子的速度
				Vector3f angularV = getAngularVelocity();
				dynamicParticles[j]->m_v = angularV.Cross(dynamicParticles[j]->m_x - getPosition()) + getVelocity();
			}
		}
	}
	
}

void DynamicRigidParticleObject::Calculation()
{
	Vector3f G(0.0, 0.0, -9.8f);
	//在当前帧计算之前先将刚体的加速度初始化为重力加速度，角速度、速度都不用初始化
	m_rigidBody.setAcceleration(G);
	//先更新当前刚体的各种属性
	UpdateObjectPar();
	//依据刚体更新后的属性对粒子属性进行刷新
	UpdateParticlePar();

}
