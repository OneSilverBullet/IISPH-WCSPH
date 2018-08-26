#pragma once
//////////////////////////////////////////////////////////////////////////
//Vector3i类
//作用：数据结构，用于测算流体方块和方格空间的长宽高的单位数量
//////////////////////////////////////////////////////////////////////////

class Vector3i
{
public :
	int x, y, z;

public:
	Vector3i() { x = 0; y = 0; z = 0; }
	Vector3i(int xValue, int yValue, int zValue) {
		x = xValue;
		y = yValue;
		z = zValue;
	}
	Vector3i(Vector3i& v)
	{
		operator=(v); //执行等于参数
	}
	~Vector3i(){}

	Vector3i& operator=( Vector3i& value);
	Vector3i& operator=(const Vector3i& a);
	Vector3i& operator+( Vector3i& value);
	Vector3i& operator-( Vector3i& value);
	Vector3i& operator*( Vector3i& value);
	Vector3i& operator*( int value);
	bool operator==( Vector3i& value);
	bool operator!=( Vector3i& value);

	int GetX();
	int GetY();
	int GetZ();
	void SetValue(int xV, int yV, int zV);
};
