#include "SPHKernel.h"

//���ھ�̬�������г�ʼ��
double CubicKernel::m_radius;
double CubicKernel::m_k;
double CubicKernel::m_l;
double CubicKernel::m_W_zero;

double Poly6Kernel::m_radius;
double Poly6Kernel::m_k;
double Poly6Kernel::m_l;
double Poly6Kernel::m_m;
double Poly6Kernel::m_W_zero;

double SpikyKernel::m_radius;
double SpikyKernel::m_k;
double SpikyKernel::m_l;
double SpikyKernel::m_W_zero;



//��Щ�˺���Ҳ���ã���Ϊ�ο�
double SPHKernal::Poly6Kernel(double h, double r)
{
	if (r>h)
	{
		return 0;
	}
	return double(315.0/(64.0*PI*pow(h, 9))*pow((h*h-r*r),3));
}

double SPHKernal::SpikyKernel(double h, double r)
{
	if (r > h)
	{
		return 0;
	}
	return float(15/PI*pow(h, 6)*pow((h-r), 3));
}

//����˺��� ����SPHSURVIVALKIT������������������û��ʹ��������˺�������
Vector3f SPHKernal::SpikyKernelGrad(double h, Vector3f origin, double distance)
{
	Vector3f result(0, 0, 0);
	if (distance > h)
	{
		return result;
	}
	return origin*(-45 / (PI*pow(h, 6))*pow((h - distance), 2));
	//������ɣ�SPH SURVIVAL KIT���У����һ��r��r��
 }

double SPHKernal::PolyNormal(double h, double r)
{
	if (r > h)
	{
		return 0;
	}
	return float(-r*r*r/(2*h*h*h)+r*r/(h*h)+h/(2*r)-1);
}

double SPHKernal::PolyNormalGrad(double h, double r)
{
	if (r > h)
	{
		return 0;
	}
	return float(45 / (PI*pow(h, 6))*(h - r));
}
