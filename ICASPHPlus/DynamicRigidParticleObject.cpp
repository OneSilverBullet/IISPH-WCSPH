#include "DynamicRigidParticleObject.h"

void DynamicRigidParticleObject::Initialize()
{
	InitParticle();
	InitRigidBody();
}

void DynamicRigidParticleObject::InitParticle()
{
	//�ӹ�ϣ�����б��еõ������ؼ��ľ�ֵ̬,���õĶ��Ǿ�ֵ̬
	int cellNum = SharedData::GetHashCellNum();
	double cellLength = SharedData::GetCellLength();
	//ֱ��
	double diam = radius * 2;
	//���úõ�ǰx��������Ӹ���
	int xLength = floor(length / diam + 0.5f) + 1;
	//���úõ�ǰ��̬��ĳ���������ӵ�����
	dynamicCubeWLH.SetValue(xLength, xLength, xLength);

	std::cout << "��ʼ����̬����" << std::endl;
	std::vector<Vector3f> boundaryParticles;


	Vector3f boundaryMinP(ox, oy, oz);
	Vector3f boundaryMaxP(ox + length, oy + length, oz + length);

	//ʹ��ֱ������Ϊÿһ����ǰ�ı仯��
	double xshift = diam;
	double yshift = diam;
	double zshift = diam;

	const int stepsX = xLength;
	const int stepsY = xLength;
	const int stepsZ = xLength;

	Vector3f start = boundaryMinP;

	//����ǽ��ʹ�����ӽ���ǽ�Ĵ���
	for (int i = 0; i < stepsX; i++)
	{
		for (int j = 0; j < stepsY; j++)
		{
			Vector3f currPos = Vector3f(i*xshift, j*yshift, 0) + boundaryMinP;
			//���ݵ�ǰ���ӵ�λ����õ�ǰ���ӵĹ�ϣֵ
			int hashValue = FunctionKit::PositionMapHash(currPos, cellLength, cellNum, boundaryOffset);
			//����һ���µ�����
			DynamicRigidParticle* newParticle = new DynamicRigidParticle(currPos, hashValue);
			//����ǰ����װ�뵽�����б���
			dynamicParticles.push_back(newParticle);
		}
	}
	for (int i = 0; i < stepsX; i++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = Vector3f(i*xshift, 0, k*zshift) + boundaryMinP;
			//���ݵ�ǰ���ӵ�λ����õ�ǰ���ӵĹ�ϣֵ
			int hashValue = FunctionKit::PositionMapHash(currPos, cellLength, cellNum, boundaryOffset);
			//����һ���µ�����
			DynamicRigidParticle* newParticle = new DynamicRigidParticle(currPos, hashValue);
			//����ǰ����װ�뵽�����б���
			dynamicParticles.push_back(newParticle);
		}
	}
	for (int j = 0; j < stepsY; j++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = Vector3f(0, j*yshift, k*zshift) + boundaryMinP;
			//���ݵ�ǰ���ӵ�λ����õ�ǰ���ӵĹ�ϣֵ
			int hashValue = FunctionKit::PositionMapHash(currPos, cellLength, cellNum, boundaryOffset);
			//����һ���µ�����
			DynamicRigidParticle* newParticle = new DynamicRigidParticle(currPos, hashValue);
			//����ǰ����װ�뵽�����б���
			dynamicParticles.push_back(newParticle);
		}
	}
	for (int i = 0; i < stepsX; i++)
	{
		for (int j = 0; j < stepsY; j++)
		{
			Vector3f currPos = -Vector3f(i*xshift, j*yshift, 0) + boundaryMaxP;
			//���ݵ�ǰ���ӵ�λ����õ�ǰ���ӵĹ�ϣֵ
			int hashValue = FunctionKit::PositionMapHash(currPos, cellLength, cellNum, boundaryOffset);
			//����һ���µ�����
			DynamicRigidParticle* newParticle = new DynamicRigidParticle(currPos, hashValue);
			//����ǰ����װ�뵽�����б���
			dynamicParticles.push_back(newParticle);
		}
	}
	for (int i = 0; i < stepsX; i++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = -Vector3f(i*xshift, 0, k*zshift) + boundaryMaxP;
			//���ݵ�ǰ���ӵ�λ����õ�ǰ���ӵĹ�ϣֵ
			int hashValue = FunctionKit::PositionMapHash(currPos, cellLength, cellNum, boundaryOffset);
			//����һ���µ�����
			DynamicRigidParticle* newParticle = new DynamicRigidParticle(currPos, hashValue);
			//����ǰ����װ�뵽�����б���
			dynamicParticles.push_back(newParticle);
		}
	}
	for (int j = 0; j < stepsY; j++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = -Vector3f(0, j*yshift, k*zshift) + boundaryMaxP;
			//���ݵ�ǰ���ӵ�λ����õ�ǰ���ӵĹ�ϣֵ
			int hashValue = FunctionKit::PositionMapHash(currPos, cellLength, cellNum, boundaryOffset);
			//����һ���µ�����
			DynamicRigidParticle* newParticle = new DynamicRigidParticle(currPos, hashValue);
			//����ǰ����װ�뵽�����б���
			dynamicParticles.push_back(newParticle);
		}
	}



}

void DynamicRigidParticleObject::InitRigidBody()
{
	//���ȸ��µ�ǰ�����λ�ã�����Ϊ��ǰ������
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
	//�������е����ӵ��������ĵľ��롣
	CalculateDiffDistance();

	//���µ�ǰ���������������Ϊ���������������ܺ�
	m_rigidBody.getMass() = dynamicParticles.size()*particleMass;
	//���㵱ǰ�Ĺ������������ֵΪ��ֵ
	CalculateInertiaTensor();
}

//��ʼ����
void DynamicRigidParticleObject::CalculateInertiaTensor()
{
	//��ʼ����������֮�󣬸��ݸ������ӶԹ����������м��㣬���ֵ��һ����ֵ
	//���������ɵ�����Ϊһ���Գ����壬���˾���Խ��߷�����������Ԫ�ض���
	//0�����ֻ�����������������ҶԸ���Ĵ洢������Vec3��ʽ
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
	//������ֵ�����ǲ����
	m_rigidBody.m_inertiaTensor.SetValue(IXX, IYY, IZZ);
	m_rigidBody.m_inertiaTensorInverse.SetValue(double(1.0 / IXX), double(1.0 / IYY), double(1.0 / IZZ));
	Matrix3f diagonalInertiaTensor;
	//�����һ�ֹ�����������������ϵ�µ��������ʱ��仯
	m_rigidBody.m_inertiaTensorInverseW = m_rigidBody.m_rot * diagonalInertiaTensor.Diagonal(m_rigidBody.m_inertiaTensorInverse) *m_rigidBody.m_rot.Tranpose();
}

//ʵʱ���ø���
void DynamicRigidParticleObject::CalculateRotationMat()
{
	//���ݽ��ٶ���ö�Ӧ����ķ������
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
		//������ǰ�������ӵ�λ��
		for (int j = 0; j < (int)dynamicParticles.size(); j++)
		{
			force += dynamicParticles[j]->m_f;
			torque += (dynamicParticles[j]->m_x - m_rigidBody.getPosition()).Cross(dynamicParticles[j]->m_f);
			//�ڼ�������󣬽���ǰ���ӵ�������Ϊ0
			dynamicParticles[j]->m_f.Zero();
		}
		//ˢ�µõ���ǰ������ٶ�
		addForce(force);
		//ˢ�µõ���ǰ����Ľ��ٶ�
		addTorque(torque);
		//���ݵ�ǰ���ٶȸ��µ�ǰ�������ת����
		CalculateRotationMat();
		if(m_rigidBody.getPosition().y <0.2)
		{
			Vector3f newPos(m_rigidBody.getPosition().x, m_rigidBody.getPosition().z, 0.2);
			m_rigidBody.setPosition(newPos);
		}
		else
		{
			//���µ�ǰ�����λ�ã����ݵ�ǰ�ٶȶ�λ�ý���ˢ��
			m_rigidBody.setPosition(m_rigidBody.getPosition() + m_rigidBody.getVelocity()*m_h);
		}
		//���ݵ�ǰ�ĵ���ת���󣬸��µ�ǰ���������µĹ������������
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
				//���ݵ�ǰ�����λ�õõ���ǰÿ�����ӵ�λ��
				dynamicParticles[j]->m_x = m_rigidBody.m_rot * dynamicParticles[j]->diffDistance + getPosition(); //����Ϊ��һ֡��λ��
				dynamicParticles[j]->m_x0 = dynamicParticles[j]->m_x;
				//���µ�ǰ�������ӵĹ�ϣֵ
				int hashValue = FunctionKit::PositionMapHash(dynamicParticles[j]->m_x, cellLength,
					cellNum, boundaryOffset);
				//�ڼ�������󣬽���ǰ���ӵ�������Ϊ0
				dynamicParticles[j]->hashValue = hashValue;

				//���ݵ�ǰ�����λ�úͽ��ٶȵõ���ǰ���ӵ��ٶ�
				Vector3f angularV = getAngularVelocity();
				dynamicParticles[j]->m_v = angularV.Cross(dynamicParticles[j]->m_x - getPosition()) + getVelocity();
			}
		}
	}
	
}

void DynamicRigidParticleObject::Calculation()
{
	Vector3f G(0.0, 0.0, -9.8f);
	//�ڵ�ǰ֡����֮ǰ�Ƚ�����ļ��ٶȳ�ʼ��Ϊ�������ٶȣ����ٶȡ��ٶȶ����ó�ʼ��
	m_rigidBody.setAcceleration(G);
	//�ȸ��µ�ǰ����ĸ�������
	UpdateObjectPar();
	//���ݸ�����º�����Զ��������Խ���ˢ��
	UpdateParticlePar();

}
