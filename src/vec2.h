#ifndef VEC2_H_INCL
#define VEC2_H_INCL

#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace la {

template<typename T>
struct vec2 {

	static const vec2 I;
	static const vec2 J;
	static const vec2 E;
	static const vec2 O;

	T x, y;

	explicit vec2(T a = 0) { set(a); }
	vec2(T X, T Y) { set(X, Y); }
	vec2(const vec2& v) { set(v); }

	vec2& set(T X, T Y) {
		x = X;
		y = Y;
		return *this;
	}

	vec2& set(T a = 0) { return set(a, a); }
	vec2& set(const vec2& v) { return set(v.x, v.y); }

	std::string str() const {
		std::ostringstream ss;
		ss << "[ " << x << " " << y << " ]";
		return ss.str();
	}

	friend std::ostream& operator<<(std::ostream& out, const vec2& v) { return out << v.str(); }

	vec2& operator=(const vec2& v) { return set(v); }
	vec2 operator-() const { return *this * (-1); }

	vec2 operator+(const vec2& v) const { return vec2(x + v.x, y + v.y); }
	vec2 operator-(const vec2& v) const { return vec2(x - v.x, y - v.y); }
	vec2 operator*(const vec2& v) const { return vec2(x * v.x, y * v.y); }
	vec2 operator/(const vec2& v) const { return vec2(x / v.x, y / v.y); }

	vec2 operator*(T a) const { return vec2(x * a, y * a); }
	vec2 operator/(T a) const { return vec2(x / a, y / a); }

	vec2& operator+=(const vec2& v) { return set(*this + v); }
	vec2& operator-=(const vec2& v) { return set(*this - v); }
	vec2& operator*=(const vec2& v) { return set(*this * v); }
	vec2& operator/=(const vec2& v) { return set(*this / v); }

	vec2& operator*=(T a) { return set(*this * a); }
	vec2& operator/=(T a) { return set(*this / a); }

	bool operator>(const vec2& v) const { return x > v.x && y > v.y; }
	bool operator>=(const vec2& v) const { return x >= v.x && y >= v.y; }
	bool operator<(const vec2& v) const { return x < v.x && y < v.y; }
	bool operator<=(const vec2& v) const { return x <= v.x && y <= v.y; }

	T dot(const vec2& v) const { return x * v.x + y * v.y; }
	T wedge(const vec2& v) const { return x * v.y - y * v.x; }

	T lensq() const { return dot(*this); }
	T len() const { return sqrt(lensq()); }
	T distsq(const vec2& v) const { return (*this - v).lensq(); }
	T dist(const vec2& v) const { return sqrt(distsq(v)); }
	vec2 norm() const { return *this / len(); }

	T angle(const vec2& v) const { return acos(dot(v) / (len() * v.len())); }
	T projlen(const vec2& v) const { return dot(v) / v.len(); }
	vec2 proj(const vec2& v) const { return v * (dot(v) / v.lensq()); }
	vec2 reflect(const vec2& v) const { return proj(v) * 2 - *this; }
	vec2 lerp(const vec2& v, T t) const { return *this + (v - *this) * t; }

	vec2 refract(const vec2& n, T r) {
		T dotn = dot(n);
		T d = 1 - r * r * (1 - dotn * dotn);
		if (d < 0) {
			return vec2();
		} else {
			return *this * r - n * (r * dotn + sqrt(d));
		}
	}



	static vec2 rand() {
		return vec2(
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX));
	}

	static vec2 randunit() {
		float phi = ::rand() / (float) RAND_MAX;
		return vec2(cos(phi), sin(phi));
	}

};

typedef vec2<float> v2;

template<typename T>
const vec2<T> vec2<T>::I(1, 0);
template<typename T>
const vec2<T> vec2<T>::J(0, 1);
template<typename T>
const vec2<T> vec2<T>::E(1, 1);
template<typename T>
const vec2<T> vec2<T>::O(0, 0);

}

#endif