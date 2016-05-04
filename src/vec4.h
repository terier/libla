#ifndef VEC4_H_INCL
#define VEC4_H_INCL

#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace la {

template<typename T>
struct vec4 {

	static const vec4 I;
	static const vec4 J;
	static const vec4 K;
	static const vec4 L;
	static const vec4 E;
	static const vec4 O;

	T x, y, z, w;

	explicit vec4(T a = 0) { set(a); }
	vec4(T X, T Y, T Z, T W = 1) { set(X, Y, Z, W); }
	vec4(const vec4& v) { set(v); }

	vec4& set(T X, T Y, T Z, T W) {
		x = X;
		y = Y;
		z = Z;
		w = W;
		return *this;
	}

	vec4& set(T a = 0) { return set(a, a, a, a); }
	vec4& set(const vec4& v) { return set(v.x, v.y, v.z, v.w); }

	std::string str() const {
		std::ostringstream ss;
		ss << "[ " << x << " " << y << " " << z << " " << w << " ]";
		return ss.str();
	}

	friend std::ostream& operator<<(std::ostream& out, const vec4& v) { return out << v.str(); }

	vec4& operator=(const vec4& v) { return set(v); }
	vec4 operator-() const { return *this * (-1); }

	vec4 operator+(const vec4& v) const { return vec4(x + v.x, y + v.y, z + v.z, w + v.w); }
	vec4 operator-(const vec4& v) const { return vec4(x - v.x, y - v.y, z - v.z, w - v.w); }
	vec4 operator*(const vec4& v) const { return vec4(x * v.x, y * v.y, z * v.z, w * v.w); }
	vec4 operator/(const vec4& v) const { return vec4(x / v.x, y / v.y, z / v.z, w / v.w); }

	vec4 operator*(T a) const { return vec4(x * a, y * a, z * a, w * a); }
	vec4 operator/(T a) const { return vec4(x / a, y / a, z / a, w / a); }

	vec4& operator+=(const vec4& v) { return set(*this + v); }
	vec4& operator-=(const vec4& v) { return set(*this - v); }
	vec4& operator*=(const vec4& v) { return set(*this * v); }
	vec4& operator/=(const vec4& v) { return set(*this / v); }

	vec4& operator*=(T a) { return set(*this * a); }
	vec4& operator/=(T a) { return set(*this / a); }

	bool operator>(const vec4& v) const { return x > v.x && y > v.y && z > v.z && w > v.w; }
	bool operator>=(const vec4& v) const { return x >= v.x && y >= v.y && z >= v.z && w >= v.w; }
	bool operator<(const vec4& v) const { return x < v.x && y < v.y && z < v.z && w < v.w; }
	bool operator<=(const vec4& v) const { return x <= v.x && y <= v.y && z <= v.z && w <= v.w; }

	T dot(const vec4& v) const { return x * v.x + y * v.y + z * v.z + w * v.w; }

	T lensq() const { return dot(*this); }
	T len() const { return sqrt(lensq()); }
	T distsq(const vec4& v) const { return (*this - v).lensq(); }
	T dist(const vec4& v) const { return sqrt(distsq(v)); }
	vec4 norm() const { return *this / len(); }
	vec4 regular() const { return *this / w; }
	T angle(const vec4& v) const { return acos(dot(v) / (len() * v.len())); }

	T projlen(const vec4& v) const { return dot(v) / v.len(); }
	vec4 proj(const vec4& v) const { return v * (dot(v) / v.lensq()); }
	vec4 reflect(const vec4& v) const { return proj(v) * 2 - *this; }
	vec4 lerp(const vec4& v, T t) const { return *this + (v - *this) * t; }

	vec4 refract(const vec4& n, T r) {
		T dotn = dot(n);
		T d = 1 - r * r * (1 - dotn * dotn);
		if (d < 0) {
			return vec4();
		} else {
			return *this * r - n * (r * dotn + sqrt(d));
		}
	}



	static vec4 rand() {
		return vec4(
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX));
	}

	static vec4 randunit() {
		return rand().norm();
	}

};

typedef vec4<float> v4;

template<typename T>
const vec4<T> vec4<T>::I(1, 0, 0, 0);
template<typename T>
const vec4<T> vec4<T>::J(0, 1, 0, 0);
template<typename T>
const vec4<T> vec4<T>::K(0, 0, 1, 0);
template<typename T>
const vec4<T> vec4<T>::L(0, 0, 0, 1);
template<typename T>
const vec4<T> vec4<T>::E(1, 1, 1, 1);
template<typename T>
const vec4<T> vec4<T>::O(0, 0, 0, 0);

}

#endif