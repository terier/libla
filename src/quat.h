#ifndef QUAT_H_INCL
#define QUAT_H_INCL

#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "vec3.h"
#include "vec4.h"
#include "mat4.h"

namespace la {

template<typename T>
struct quat {

	static const quat I;
	static const quat J;
	static const quat K;

	T x, y, z, w;

	quat(T a = 0) { set(a); }
	quat(T X, T Y, T Z, T W) { set(X, Y, Z, W); }
	quat(const quat& q) { set(q); }
	quat(const vec4<T>& v) { set(v); }
	quat(const vec3<T>& v) { set(0, v.x, v.y, v.z); }
	quat(T a, const vec3<T>& v) { set(a, v.x, v.y, v.z); }

	quat& set(T X, T Y, T Z, T W) {
		x = X;
		y = Y;
		z = Z;
		w = W;
		return *this;
	}

	quat& set(T a = 0) { return set(a, 0, 0, 0); }
	quat& set(const quat& q) { return set(q.x, q.y, q.z, q.w); }
	quat& set(const vec4<T>& v) { return set(v.x, v.y, v.z, v.w); }

	std::string str() const {
		std::ostringstream ss;
		ss << x << " + " << y << "i + " << z << "j + " << w << "k";
		return ss.str();
	}

	friend std::ostream& operator<<(std::ostream& out, const quat& q) { return out << q.str(); }

	quat& operator=(const quat& q) { return set(q); }
	quat operator-() const { return *this * (-1); }

	quat operator+(const quat& q) const { return quat(x + q.x, y + q.y, z + q.z, w + q.w); }
	quat operator-(const quat& q) const { return quat(x - q.x, y - q.y, z - q.z, w - q.w); }
	quat operator*(const quat& q) const {
		return quat(
			x * q.x - y * q.y - z * q.z - w * q.w,
			x * q.y + y * q.x + z * q.w - w * q.z,
			x * q.z - y * q.w + z * q.x + w * q.y,
			x * q.w + y * q.z - z * q.y + w * q.x);
	}
	quat operator/(const quat& q) const { return *this * q.inv(); }

	quat operator*(T a) const { return quat(x * a, y * a, z * a, w * a); }
	quat operator/(T a) const { return quat(x / a, y / a, z / a, w / a); }

	quat& operator+=(const quat& q) { return set(*this + q); }
	quat& operator-=(const quat& q) { return set(*this - q); }
	quat& operator*=(const quat& q) { return set(*this * q); }
	quat& operator/=(const quat& q) { return set(*this / q); }

	quat& operator*=(T a) { return set(*this * a); }
	quat& operator/=(T a) { return set(*this / a); }

	quat conj() const { return quat(x, -y, -z, -w); }
	T lensq() const { return x * x + y * y + z * z + w * w; }
	T len() const { return sqrt(lensq()); }
	T distsq(const quat& q) const { return (*this - q).lensq(); }
	T dist(const quat& q) const { return sqrt(distsq(q)); }
	quat norm() const { return *this / len(); }
	quat inv() const { return conj() / lensq(); }
	T real() const { return x; }
	vec3<T> imag() const { return vec3<T>(y, z, w); }
	vec4<T> vec() const { return vec4<T>(x, y, z, w); }

	vec3<T> euler() const {
		return vec3<T>(
			atan2(2 * (x * y + z * w), 1 - 2 * (y * y + z * z)),
			asin(2 * (x * z - w * y)),
			atan2(2 * (x * w + y * z), 1 - 2 * (z * z + w * w)));
	}

	quat exp() const {
		vec3<T> v = imag();
		T vn = v.len();
		return quat(cos(vn), v * (sin(vn) / vn)) * exp(x);
	}

	mat4<T> mat() const {
		return mat4<T>(
			x, y, z, w,
			-y, x, -w, z,
			-z, w, x, -y,
			-w, -z, y, x);
	}

	mat4<T> rmat() const {
		return mat4<T>(
			1 - 2 * (z * z - w * w), 2 * (y * z - w * x), 2 * (y * w + z * x), 0,
			2 * (y * z + w * x), 1 - 2 * (y * y - w * w), 2 * (z * w - y * x), 0,
			2 * (y * w - z * x), 2 * (z * w + y * x), 1 - 2 * (y * y - z * z), 0,
			0, 0, 0, 1);

	}

	vec3<T> rotate(const vec3<T>& v) { return (*this * quat(v) * inv()).imag(); }



	static quat rand() {
		return quat(
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX));
	}

	static quat randnorm() { return rand().norm(); }

	static quat randrot() {
		quat r;
		do { r = rand(); } while (r.lensq() > 1);
		return r.norm();
	}

	static quat rotX(T a) { return quat(cos(a / 2), sin(a / 2), 0, 0); }
	static quat rotY(T a) { return quat(cos(a / 2), 0, sin(a / 2), 0); }
	static quat rotZ(T a) { return quat(cos(a / 2), 0, 0, sin(a / 2)); }
	static quat rot(const vec3<T> v) { return rotZ(v.z) * rotY(v.y) * rotZ(v.z); }
	static quat rot(T a, const vec3<T> v) { return quat(cos(a / 2), v * sin(a / 2)); }

};

typedef quat<float> quatf;

template<typename T>
const quat<T> quat<T>::I(0, 1, 0, 0);
template<typename T>
const quat<T> quat<T>::J(0, 0, 1, 0);
template<typename T>
const quat<T> quat<T>::K(0, 0, 0, 1);

}

#endif