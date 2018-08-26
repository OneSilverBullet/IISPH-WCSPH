#include "Vector3i.h"

Vector3i & Vector3i::operator=(Vector3i & value)
{
	this->x = value.GetX();
	this->y = value.GetY();
	this->z = value.GetZ();
	return *this;
}

Vector3i& Vector3i::operator=(const Vector3i& value)
{
	this->x = value.x;
	this->y = value.y;
	this->z = value.z;
	return *this;
}

Vector3i & Vector3i::operator+(Vector3i & value)
{
	this->x += value.GetX();
	this->y += value.GetY();
	this->z += value.GetZ();
	return *this;
}

Vector3i & Vector3i::operator-(Vector3i & value)
{
	this->x -= value.GetX();
	this->y -= value.GetY();
	this->z -= value.GetZ();
	return *this;
}

Vector3i & Vector3i::operator*(Vector3i & value)
{
	this->x *= value.GetX();
	this->y *= value.GetY();
	this->z *= value.GetZ();
	return *this;
}

Vector3i & Vector3i::operator*(int value)
{
	this->x *= value;
	this->y *= value;
	this->z *= value;
	return *this;
}


bool Vector3i::operator==(Vector3i & value)
{
	return (this->x == value.GetX() && this->y == value.GetY() && this->z == value.GetZ());
}

bool Vector3i::operator!=(Vector3i & value)
{
	return !operator==(value);
}

int Vector3i::GetX()
{
	return x;
}

int Vector3i::GetY()
{
	return y;
}

int Vector3i::GetZ()
{
	return z;
}

void Vector3i::SetValue(int xV, int yV, int zV)
{
	x = xV;
	y = yV;
	z = zV;
}
