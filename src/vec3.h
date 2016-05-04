#ifndef VEC3_H_INCL
#define VEC3_H_INCL

#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace la {

template<typename T>
struct vec3 {

	static const vec3 I;
	static const vec3 J;
	static const vec3 K;
	static const vec3 E;
	static const vec3 O;

	T x, y, z;

	explicit vec3(T a = 0) { set(a); }
	vec3(T X, T Y, T Z) { set(X, Y, Z); }
	vec3(const vec3& v) { set(v); }

	vec3& set(T X, T Y, T Z) {
		x = X;
		y = Y;
		z = Z;
		return *this;
	}

	vec3& set(T a = 0) { return set(a, a, a); }
	vec3& set(const vec3& v) { return set(v.x, v.y, v.z); }

	std::string str() const {
		std::ostringstream ss;
		ss << "[ " << x << " " << y << " " << z << " ]";
		return ss.str();
	}

	friend std::ostream& operator<<(std::ostream& out, const vec3& v) { return out << v.str(); }

	vec3& operator=(const vec3& v) { return set(v); }
	vec3 operator-() const { return *this * (-1); }

	vec3 operator+(const vec3& v) const { return vec3(x + v.x, y + v.y, z + v.z); }
	vec3 operator-(const vec3& v) const { return vec3(x - v.x, y - v.y, z - v.z); }
	vec3 operator*(const vec3& v) const { return vec3(x * v.x, y * v.y, z * v.z); }
	vec3 operator/(const vec3& v) const { return vec3(x / v.x, y / v.y, z / v.z); }

	vec3 operator*(T a) const { return vec3(x * a, y * a, z * a); }
	vec3 operator/(T a) const { return vec3(x / a, y / a, z / a); }

	vec3& operator+=(const vec3& v) { return set(*this + v); }
	vec3& operator-=(const vec3& v) { return set(*this - v); }
	vec3& operator*=(const vec3& v) { return set(*this * v); }
	vec3& operator/=(const vec3& v) { return set(*this / v); }

	vec3& operator*=(T a) { return set(*this * a); }
	vec3& operator/=(T a) { return set(*this / a); }

	bool operator>(const vec3& v) const { return x > v.x && y > v.y && z > v.z; }
	bool operator>=(const vec3& v) const { return x >= v.x && y >= v.y && z >= v.z; }
	bool operator<(const vec3& v) const { return x < v.x && y < v.y && z < v.z; }
	bool operator<=(const vec3& v) const { return x <= v.x && y <= v.y && z <= v.z; }

	T dot(const vec3& v) const { return x * v.x + y * v.y + z * v.z; }
	vec3 cross(const vec3& v) const { return vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }

	T lensq() const { return dot(*this); }
	T len() const { return sqrt(lensq()); }
	T distsq(const vec3& v) const { return (*this - v).lensq(); }
	T dist(const vec3& v) const { return sqrt(distsq(v)); }
	vec3 norm() const { return *this / len(); }
	T angle(const vec3& v) const { return acos(dot(v) / (len() * v.len())); }

	T projlen(const vec3& v) const { return dot(v) / v.len(); }
	vec3 proj(const vec3& v) const { return v * (dot(v) / v.lensq()); }
	vec3 reflect(const vec3& v) const { return proj(v) * 2 - *this; }
	vec3 lerp(const vec3& v, T t) const { return *this + (v - *this) * t; }

	vec3 refract(const vec3& n, T r) {
		T dotn = dot(n);
		T d = 1 - r * r * (1 - dotn * dotn);
		if (d < 0) {
			return vec3();
		} else {
			return *this * r - n * (r * dotn + sqrt(d));
		}
	}



	static vec3 rand() {
		return vec3(
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX));
	}

	static vec3 randunit() {
		vec3 v;
		do { v = vec3::rand(); } while (v.lensq() > 1);
		return v.norm();
	}

};

typedef vec3<float> v3;

template<typename T>
const vec3<T> vec3<T>::I(1, 0, 0);
template<typename T>
const vec3<T> vec3<T>::J(0, 1, 0);
template<typename T>
const vec3<T> vec3<T>::K(0, 0, 1);
template<typename T>
const vec3<T> vec3<T>::E(1, 1, 1);
template<typename T>
const vec3<T> vec3<T>::O(0, 0, 0);

}

#endif