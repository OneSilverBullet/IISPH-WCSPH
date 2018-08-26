#include "SPHKernel.h"

//对于静态变量进行初始化
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



//这些核函数也有用，作为参考
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

//这个核函数 按照SPHSURVIVALKIT所做，在真正程序当中没有使用这个，核函数待定
Vector3f SPHKernal::SpikyKernelGrad(double h, Vector3f origin, double distance)
{
	Vector3f result(0, 0, 0);
	if (distance > h)
	{
		return result;
	}
	return origin*(-45 / (PI*pow(h, 6))*pow((h - distance), 2));
	//这里存疑，SPH SURVIVAL KIT当中，最后一个r是r拔
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
