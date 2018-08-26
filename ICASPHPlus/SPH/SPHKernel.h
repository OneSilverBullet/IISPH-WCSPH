#pragma once
#include <math.h>
#include "ShareData.h"
#include "Vector3f.h"
#include "FunctionKit.h"

//////////////////////////////////////////////////////////////////////////
//SPHKernal类
//作用：包含了SPH的基本核函数，但是在具体计算当中没有使用该类的核函数。
//////////////////////////////////////////////////////////////////////////
class SPHKernal
{
public:

	//使用POLY6核函数对粒子密度进行计算优化
	static double Poly6Kernel(double h, double r);
	//使用Spiky Kernel对压力梯度进行计算
	static double SpikyKernel(double h, double r);
    //使用Spiky Kernel的梯度值
	static Vector3f SpikyKernelGrad(double h, Vector3f origin, double distance);
	//针对黏度，我们使用另外一种PolyNormal核函数
	static double PolyNormal(double h, double r);
	//PolyNormal的拉普拉斯算子，计算黏度的核函数的拉普拉斯算子
	static double PolyNormalGrad(double h, double r);
};

//////////////////////////////////////////////////////////////////////////
//CubicKernelOne类
//作用：针对每个粒子进行操作的核函数
//////////////////////////////////////////////////////////////////////////
class CubicKernelOne
{
public:
	static double W(double m_radius, const double r)
	{
		double h3 = m_radius*m_radius*m_radius;
		double 	m_k = 8.0 / (PI*h3);
		double res = 0.0;
		const double q = r / m_radius;
		if (q <= 1.0)
		{
			if (q <= 0.5)
			{
				const double q2 = q*q;
				const double q3 = q2*q;
				res = m_k * (6.0*q3 - 6.0*q2 + 1.0);
			}
			else
			{
				res = m_k * (2.0*pow(1.0 - q, 3));
			}
		}
		return res;
	}

	static double W(double radius, const Vector3f &r)
	{
		return W(radius, r.Norm());
	}

	//这个函数在WCSPH的压力计算当中使用
	static Vector3f gradW(double m_radius, const Vector3f &r)
	{
		Vector3f res;
		double h3 = m_radius*m_radius*m_radius;
		double m_l = 48.0 / (PI*h3);
		const double rl = r.Norm();
		const double q = rl / m_radius;
		if (q <= 1.0)
		{
			if (rl > 1.0e-6)
			{
				const Vector3f gradq = r * ((double) 1.0 / (rl*m_radius));
				if (q <= 0.5)
				{
					res = m_l*q*((double) 3.0*q - (double) 2.0)*gradq;
				}
				else
				{
					const double factor = 1.0 - q;
					res = m_l*(-factor*factor)*gradq;
				}
			}
		}
		else
			res.Zero();

		return res;
	}

	static double W_zero(double radius)
	{
		Vector3f a(0,0,0);
		return W(radius, a);
	}
};


//////////////////////////////////////////////////////////////////////////
//CubicKernel类
//作用：包含了Cubic SPH的基本核函数，在计算当中使用了这个类的相应的核函数。
//////////////////////////////////////////////////////////////////////////
class CubicKernel
{
protected:
	static double m_radius;
	static double m_k;
	static double m_l;
	static double m_W_zero;
public:
	static double getRadius() { return m_radius; }
	static void setRadius(double val)
	{
		m_radius = val;
		static const double pi = static_cast<double>(PI);

		const double h3 = m_radius*m_radius*m_radius;
		m_k = 8.0 / (pi*h3);
		m_l = 48.0 / (pi*h3);
		Vector3f a(0, 0, 0);
		m_W_zero = W(a);
	}

public:
	static double W(const double r)
	{
		double res = 0.0;
		const double q = r / m_radius;
		if (q <= 1.0)
		{
			if (q <= 0.5)
			{
				const double q2 = q*q;
				const double q3 = q2*q;
				res = m_k * (6.0*q3 - 6.0*q2 + 1.0);
			}
			else
			{
				res = m_k * (2.0*pow(1.0 - q, 3));
			}
		}
		return res;
	}

	static double W(const Vector3f &r)
	{
		return W(r.Norm());
	}

	//这个函数在WCSPH的压力计算当中使用
	static Vector3f gradW(const Vector3f &r)
	{
		Vector3f res;
		const double rl = r.Norm();
		const double q = rl / m_radius;
		if (q <= 1.0)
		{
			if (rl > 1.0e-6)
			{
				const Vector3f gradq = r * ((double) 1.0 / (rl*m_radius));
				if (q <= 0.5)
				{
					res = m_l*q*((double) 3.0*q - (double) 2.0)*gradq;
				}
				else
				{
					const double factor = 1.0 - q;
					res = m_l*(-factor*factor)*gradq;
				}
			}
		}
		else
			res.Zero();

		return res;
	}

	static double W_zero()
	{
		return m_W_zero;
	}
};

//////////////////////////////////////////////////////////////////////////
//Poly6KernelOne类
//作用：包含了Poly6Kernel的基本核函数，适用于每个不同半径粒子的使用
//////////////////////////////////////////////////////////////////////////
class Poly6KernelOne
{
public:

	/**
	* W(r,h) = (315/(64 pi h^9))(h^2-|r|^2)^3
	*        = (315/(64 pi h^9))(h^2-r*r)^3
	*/
	static double W(const double h, const double r)
	{
		double res = 0.0;
		const double r2 = r*r;
		const double radius2 = h*h;
		const double m_k = 315.0 / (64.0*PI*pow(h, 9));
		if (r2 <= radius2)
		{
			res = pow(radius2 - r2, 3)*m_k;
		}
		return res;
	}

	static double W(const double h, Vector3f &r)
	{
		double res = 0.0;
		const double r2 = r.SquaredNorm();
		const double radius2 = h*h;
		const double m_k = 315.0 / (64.0*PI*pow(h, 9));
		if (r2 <= radius2)
		{
			res = pow(radius2 - r2, 3)*m_k;
		}
		return res;
	}


	/**
	* grad(W(r,h)) = r(-945/(32 pi h^9))(h^2-|r|^2)^2
	*              = r(-945/(32 pi h^9))(h^2-r*r)^2
	*/
	static Vector3f gradW(const double h, const Vector3f &r)
	{
		Vector3f res;
		const double r2 = r.SquaredNorm();
		const double radius2 = h*h;
		const double m_l = -945.0 / (32.0*PI*pow(h, 9));
		if (r2 <= radius2)
		{
			double tmp = radius2 - r2;
			res = m_l * tmp * tmp * r;
		}
		else
			res.Zero();

		return res;
	}

	/**
	* laplacian(W(r,h)) = (-945/(32 pi h^9))(h^2-|r|^2)(-7|r|^2+3h^2)
	*                   = (-945/(32 pi h^9))(h^2-r*r)(3 h^2-7 r*r)
	*/
	static double laplacianW(const double h, const Vector3f &r)
	{
		double res;
		const double r2 = r.SquaredNorm();
		const double radius2 = h*h;
		const double m_l = -945.0 / (32.0*PI*pow(h, 9));
		const double m_m = m_l;
		if (r2 <= radius2)
		{
			double tmp = radius2 - r2;
			double tmp2 = 3 * radius2 - 7 * r2;
			res = m_m * tmp  * tmp2;
		}
		else
			res = (double) 0.;

		return res;
	}

	//W对h求偏导数
	static double WPartialDerivativeH(const double h, const Vector3f r)
	{
		const double m_k = 315.0 / (64.0*PI*pow(h, 10));
		const double r2 = r.SquaredNorm();
		const double h2 = h*h;
		double res;
		if (r2 <= h2)
		{
			double tmp = -3.0*pow(h, 6) + 15 * r2*h*h*h*h - 21 * h*h*r2*r2 + 9 * pow(r2, 3);
			res = m_k*tmp;
		}
		else
			res = (double) 0.;

		return res;
	}

	//W对h求偏导数
	static double WPartialDerivativeH(const double h, const double r)
	{
		const double m_k = 315.0 / (64.0*PI*pow(h, 10));
		const double r2 = r*r;
		const double h2 = h*h;
		double res;
		if (r2 <= h2)
		{
			double tmp = -3.0*pow(h, 6) + 15 * r*r*h*h*h*h - 21 * h*h*r*r*r*r + 9 * pow(r, 6);
			res = m_k*tmp;
		}
		else
			res = (double) 0.;

		return res;
	}

	//对w_zero求偏导数
	static double WZeroPartialDerivativeH(const double h)
	{
		double  m_W_zero = WPartialDerivativeH(h, 0);
		return m_W_zero;
	}

	static double W_zero(const double h)
	{
		Vector3f temp;
		temp.Zero();
		double m_W_zero = W(h, temp);
		return m_W_zero;
	}
};


//////////////////////////////////////////////////////////////////////////
//Poly6Kernel类
//作用：包含了Poly6Kernel的基本核函数，在计算当没有用了这个类的相应的核函数。
//////////////////////////////////////////////////////////////////////////
class Poly6Kernel
{
protected:
	static double m_radius;
	static double m_k;
	static double m_l;
	static double m_m;
	static double m_W_zero;
public:
	static double getRadius() { return m_radius; }
	static void setRadius(double val)
	{
		m_radius = val;
		static const double pi = static_cast<double>(PI);
		m_k = 315.0 / (64.0*pi*pow(m_radius, 9));
		m_l = -945.0 / (32.0*pi*pow(m_radius, 9));
		m_m = m_l;
		Vector3f temp;
		temp.Zero();
		m_W_zero = W(temp);
	}

public:

	/**
	* W(r,h) = (315/(64 pi h^9))(h^2-|r|^2)^3
	*        = (315/(64 pi h^9))(h^2-r*r)^3
	*/
	static double W(const double r)
	{
		double res = 0.0;
		const double r2 = r*r;
		const double radius2 = m_radius*m_radius;
		if (r2 <= radius2)
		{
			res = pow(radius2 - r2, 3)*m_k;
		}
		return res;
	}

	static double W( Vector3f &r)
	{
		double res = 0.0;
		const double r2 = r.SquaredNorm();
		const double radius2 = m_radius*m_radius;
		if (r2 <= radius2)
		{
			res = pow(radius2 - r2, 3)*m_k;
		}
		return res;
	}


	/**
	* grad(W(r,h)) = r(-945/(32 pi h^9))(h^2-|r|^2)^2
	*              = r(-945/(32 pi h^9))(h^2-r*r)^2
	*/
	static Vector3f gradW(const Vector3f &r)
	{
		Vector3f res;
		const double r2 = r.SquaredNorm();
		const double radius2 = m_radius*m_radius;
		if (r2 <= radius2)
		{
			double tmp = radius2 - r2;
			res = m_l * tmp * tmp * r;
		}
		else
			res.Zero();

		return res;
	}

	/**
	* laplacian(W(r,h)) = (-945/(32 pi h^9))(h^2-|r|^2)(-7|r|^2+3h^2)
	*                   = (-945/(32 pi h^9))(h^2-r*r)(3 h^2-7 r*r)
	*/
	static double laplacianW(const Vector3f &r)
	{
		double res;
		const double r2 = r.SquaredNorm();
		const double radius2 = m_radius*m_radius;
		if (r2 <= radius2)
		{
			double tmp = radius2 - r2;
			double tmp2 = 3 * radius2 - 7 * r2;
			res = m_m * tmp  * tmp2;
		}
		else
			res = (double) 0.;

		return res;
	}

	static double W_zero()
	{
		return m_W_zero;
	}
};


//////////////////////////////////////////////////////////////////////////
//SpikyKernel类
//作用：包含了SpikyKernel的基本核函数，在计算当没有用了这个类的相应的核函数。
//////////////////////////////////////////////////////////////////////////
class SpikyKernel
{
protected:
	static double m_radius;
	static double m_k;
	static double m_l;
	static double m_W_zero;
public:
	static double getRadius() { return m_radius; }
	static void setRadius(double val)
	{
		m_radius = val;
		const double radius6 = pow(m_radius, 6);
		static const double pi = static_cast<double>(P1);
		m_k = 15.0 / (pi*radius6);
		m_l = -45.0 / (pi*radius6);
		Vector3f temp;
		temp.Zero();
		m_W_zero = W(temp);
	}

public:

	/**
	* W(r,h) = 15/(pi*h^6) * (h-r)^3
	*/
	static double W(const double r)
	{
		double res = 0.0;
		const double r2 = r*r;
		const double radius2 = m_radius*m_radius;
		if (r2 <= radius2)
		{
			const double hr3 = pow(m_radius - r, 3);
			res = m_k * hr3;
		}
		return res;
	}

	static double W(const Vector3f &r)
	{
		double res = 0.0;
		const double r2 = r.SquaredNorm();
		const double radius2 = m_radius*m_radius;
		if (r2 <= radius2)
		{
			const double hr3 = pow(m_radius - sqrt(r2), 3);
			res = m_k * hr3;
		}
		return res;
	}


	/**
	* grad(W(r,h)) = -r(45/(pi*h^6) * (h-r)^2)
	*/
	static Vector3f gradW(const Vector3f &r)
	{
		Vector3f res;
		const double r2 = r.SquaredNorm();
		const double radius2 = m_radius*m_radius;
		if (r2 <= radius2)
		{
			const double r_l = sqrt(r2);
			const double hr = m_radius - r_l;
			const double hr2 = hr*hr;
			res = m_l * hr2 * r * (static_cast<double>(1.0) / r_l);
		}
		else
			res.Zero();

		return res;
	}

	static double W_zero()
	{
		return m_W_zero;
	}
};