#ifndef MAT3_H_INCL
#define MAT3_H_INCL

#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "vec3.h"

namespace la {

template<typename T>
struct mat3 {

	static const mat3 I;
	static const mat3 E;

	T M[9];

	explicit mat3(T a = 0) { set(a); }
	mat3(T* m) { set(m); }
	mat3(T m11, T m12, T m13, T m21, T m22, T m23, T m31, T m32, T m33) {
		set(m11, m12, m13, m21, m22, m23, m31, m32, m33);
	}
	mat3(const mat3& m) { set(m); }

	mat3& set(T a) { return set(mat3::E * a); }
	mat3& set(const T* m) {
		for (int i = 0; i < 9; ++i) {
			M[i] = m[i];
		}
		return *this;
	}
	
	mat3& set(T m11, T m12, T m13, T m21, T m22, T m23, T m31, T m32, T m33) {
		M[0] = m11;
		M[1] = m12;
		M[2] = m13;
		M[3] = m21;
		M[4] = m22;
		M[5] = m23;
		M[6] = m31;
		M[7] = m32;
		M[8] = m33;
		return *this;
	}

	mat3& set(const mat3& m) { return set(m.M); }

	std::string str() const {
		std::ostringstream ss;
		ss << "[ " << M[0] << ", " << M[1] << ", " << M[2]
			<< "; " << M[3] << ", " << M[4] << ", " << M[5]
			<< "; " << M[6] << ", " << M[7] << ", " << M[8] << " ]";
		return ss.str();
	}

	friend std::ostream& operator<<(std::ostream& out, const mat3& m) { return out << m.str(); }

	mat3& operator=(const mat3& m) { return set(m); }
	mat3 operator-() const { return *this * (-1); }

	mat3 operator+(const mat3& m) const {
		return mat3(
			M[0] + m.M[0], M[1] + m.M[1], M[2] + m.M[2],
			M[3] + m.M[3], M[4] + m.M[4], M[5] + m.M[5],
			M[6] + m.M[6], M[7] + m.M[7], M[8] + m.M[8]);
	}

	mat3 operator-(const mat3& m) const {
		return mat3(
			M[0] - m.M[0], M[1] - m.M[1], M[2] - m.M[2],
			M[3] - m.M[3], M[4] - m.M[4], M[5] - m.M[5],
			M[6] - m.M[6], M[7] - m.M[7], M[8] - m.M[8]);
	}

	mat3 operator*(T a) const {
		return mat3(
			M[0] * a, M[1] * a, M[2] * a,
			M[3] * a, M[4] * a, M[5] * a,
			M[6] * a, M[7] * a, M[8] * a);
	}

	mat3 operator/(T a) const {
		return *this * (1 / a);
	}

	mat3 operator*(const mat3& m) const {
		return mat3(
			M[0] * m.M[0] + M[1] * m.M[3] + M[2] * m.M[6],
			M[0] * m.M[1] + M[1] * m.M[4] + M[2] * m.M[7],
			M[0] * m.M[2] + M[1] * m.M[5] + M[2] * m.M[8],
			M[3] * m.M[0] + M[4] * m.M[3] + M[5] * m.M[6],
			M[3] * m.M[1] + M[4] * m.M[4] + M[5] * m.M[7],
			M[3] * m.M[2] + M[4] * m.M[5] + M[5] * m.M[8],
			M[6] * m.M[0] + M[7] * m.M[3] + M[8] * m.M[6],
			M[6] * m.M[1] + M[7] * m.M[4] + M[8] * m.M[7],
			M[6] * m.M[2] + M[7] * m.M[5] + M[8] * m.M[8]);
	}

	vec3<T> operator*(const vec3<T>& v) const {
		return vec3<T>(
			M[0] * v.x + M[1] * v.y + M[2] * v.z,
			M[3] * v.x + M[4] * v.y + M[5] * v.z,
			M[6] * v.x + M[7] * v.y + M[8] * v.z);
	}

	vec3<T> operator/(const vec3<T>& v) const { return inv() * v; }

	T det() const {
		return M[0] * M[4] * M[8]
			+ M[1] * M[5] * M[6]
			+ M[2] * M[3] * M[7]
			- M[0] * M[5] * M[7]
			- M[1] * M[3] * M[8]
			- M[2] * M[4] * M[6];
	}

	mat3 inv() const {
		return mat3(
			M[4] * M[8] - M[5] * M[7],
			M[2] * M[7] - M[1] * M[8],
			M[1] * M[5] - M[2] * M[4],
			M[5] * M[6] - M[3] * M[8],
			M[0] * M[8] - M[2] * M[6],
			M[2] * M[3] - M[0] * M[5],
			M[3] * M[7] - M[4] * M[6],
			M[1] * M[6] - M[0] * M[7],
			M[0] * M[4] - M[1] * M[3]) / det();
	}

	mat3 transpose() const { return mat3(M[0], M[3], M[6], M[1], M[4], M[7], M[2], M[5], M[8]); }
	T tr() const { return M[0] + M[4] + M[8]; }
	vec3<T> solve(const vec3<T>& b) const { return inv() * b; }

	vec3<T> col1() const { return vec3<T>(M[0], M[3], M[6]); }
	vec3<T> col2() const { return vec3<T>(M[1], M[4], M[7]); }
	vec3<T> col3() const { return vec3<T>(M[2], M[5], M[8]); }
	vec3<T> row1() const { return vec3<T>(M[0], M[1], M[2]); }
	vec3<T> row2() const { return vec3<T>(M[3], M[4], M[5]); }
	vec3<T> row3() const { return vec3<T>(M[6], M[7], M[8]); }



	static mat3 rand() {
		return mat3(
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX));
	}

	static mat3 rotX(T phi) {
		T s = sin(phi);
		T c = cos(phi);
		return mat3(
			1, 0, 0,
			0, c, -s,
			0, s, c);
	}

	static mat3 rotY(T phi) {
		T s = sin(phi);
		T c = cos(phi);
		return mat3(
			c, 0, s,
			0, 1, 0,
			-s, 0, c);
	}

	static mat3 rotZ(T phi) {
		T s = sin(phi);
		T c = cos(phi);
		return mat3(
			c, -s, 0,
			s, c, 0,
			0, 0, 1);
	}

	static mat3 rot(const vec3<T>& rot) {
		return rotX(rot.x) * rotY(rot.y) * rotZ(rot.z);
	}

	static mat3 scale(const vec3<T>& v) {
		return mat3(
			v.x, 0, 0,
			0, v.y, 0,
			0, 0, v.z);
	}

	static mat3 fromRows(const vec3<T>& v1, const vec3<T>& v2, const vec3<T>& v3) {
		return mat3(
			v1.x, v1.y, v1.z,
			v2.x, v2.y, v2.z,
			v3.x, v3.y, v3.z);
	}

	static mat3 fromColumns(const vec3<T>& v1, const vec3<T>& v2, const vec3<T>& v3) {
		return mat3(
			v1.x, v2.x, v3.x,
			v1.y, v2.y, v3.y,
			v1.z, v2.z, v3.z);
	}

};

typedef mat3<float> m3;

template<typename T>
const mat3<T> mat3<T>::I(1, 0, 0, 0, 1, 0, 0, 0, 1);
template<typename T>
const mat3<T> mat3<T>::E(1, 1, 1, 1, 1, 1, 1, 1, 1);

}

#endif