#ifndef MAT4_H_INCL
#define MAT4_H_INCL

#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "vec4.h"
#include "vec3.h"

namespace la {

template<typename T>
struct mat4 {

	static const mat4 I;
	static const mat4 E;

	T M[16];

	explicit mat4(T a = 0) { set(a); }
	mat4(T* m) { set(m); }
	mat4(T m11, T m12, T m13, T m14, T m21, T m22, T m23, T m24, T m31, T m32, T m33, T m34, T m41, T m42, T m43, T m44) {
		set(m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44);
	}
	mat4(const mat4& m) { set(m); }

	mat4& set(T a) { return set(mat4::E * a); }
	
	mat4& set(const T* m) {
		for (int i = 0; i < 16; ++i) {
			M[i] = m[i];
		}
		return *this;
	}

	mat4& set(T m11, T m12, T m13, T m14, T m21, T m22, T m23, T m24, T m31, T m32, T m33, T m34, T m41, T m42, T m43, T m44) {
		M[ 0] = m11;
		M[ 1] = m12;
		M[ 2] = m13;
		M[ 3] = m14;
		M[ 4] = m21;
		M[ 5] = m22;
		M[ 6] = m23;
		M[ 7] = m24;
		M[ 8] = m31;
		M[ 9] = m32;
		M[10] = m33;
		M[11] = m34;
		M[12] = m41;
		M[13] = m42;
		M[14] = m43;
		M[15] = m44;
		return *this;
	}

	mat4& set(const mat4& m) { return set(m.M); }

	std::string str() const {
		std::ostringstream ss;
		ss << "[ " << M[0] << ", " << M[1] << ", " << M[2] << ", " << M[3]
			<< "; " << M[4] << ", " << M[5] << ", " << M[6] << ", " << M[7]
			<< "; " << M[8] << ", " << M[9] << ", " << M[10] << ", " << M[11]
			<< "; " << M[12] << ", " << M[13] << ", " << M[14] << ", " << M[15] << " ]";
		return ss.str();
	}

	friend std::ostream& operator<<(std::ostream& out, const mat4& m) { return out << m.str(); }

	mat4& operator=(const mat4& m) { return set(m); }
	mat4 operator-() const { return *this * (-1); }

	mat4 operator+(const mat4& m) const {
		return mat4(
			M[ 0] + m.M[ 0], M[ 1] + m.M[ 1], M[ 2] + m.M[ 2], M[ 3] + m.M[ 3],
			M[ 4] + m.M[ 4], M[ 5] + m.M[ 5], M[ 6] + m.M[ 6], M[ 7] + m.M[ 7],
			M[ 8] + m.M[ 8], M[ 9] + m.M[ 9], M[10] + m.M[10], M[11] + m.M[11],
			M[12] + m.M[12], M[13] + m.M[13], M[14] + m.M[14], M[15] + m.M[15]);
	}

	mat4 operator-(const mat4& m) const {
		return mat4(
			M[ 0] - m.M[ 0], M[ 1] - m.M[ 1], M[ 2] - m.M[ 2], M[ 3] - m.M[ 3],
			M[ 4] - m.M[ 4], M[ 5] - m.M[ 5], M[ 6] - m.M[ 6], M[ 7] - m.M[ 7],
			M[ 8] - m.M[ 8], M[ 9] - m.M[ 9], M[10] - m.M[10], M[11] - m.M[11],
			M[12] - m.M[12], M[13] - m.M[13], M[14] - m.M[14], M[15] - m.M[15]);
	}

	mat4 operator*(T a) const {
		return mat4(
			M[ 0] * a, M[ 1] * a, M[ 2] * a, M[ 3] * a,
			M[ 4] * a, M[ 5] * a, M[ 6] * a, M[ 7] * a,
			M[ 8] * a, M[ 9] * a, M[10] * a, M[11] * a,
			M[12] * a, M[13] * a, M[14] * a, M[15] * a);
	}

	mat4 operator/(T a) const {
		return *this * (1 / a);
	}

	mat4 operator*(const mat4& m) const {
		return mat4(
			M[ 0] * m.M[ 0] + M[ 1] * m.M[ 4] + M[ 2] * m.M[ 8] + M[ 3] * m.M[12],
			M[ 0] * m.M[ 1] + M[ 1] * m.M[ 5] + M[ 2] * m.M[ 9] + M[ 3] * m.M[13],
			M[ 0] * m.M[ 2] + M[ 1] * m.M[ 6] + M[ 2] * m.M[10] + M[ 3] * m.M[14],
			M[ 0] * m.M[ 3] + M[ 1] * m.M[ 7] + M[ 2] * m.M[11] + M[ 3] * m.M[15],
			M[ 4] * m.M[ 0] + M[ 5] * m.M[ 4] + M[ 6] * m.M[ 8] + M[ 7] * m.M[12],
			M[ 4] * m.M[ 1] + M[ 5] * m.M[ 5] + M[ 6] * m.M[ 9] + M[ 7] * m.M[13],
			M[ 4] * m.M[ 2] + M[ 5] * m.M[ 6] + M[ 6] * m.M[10] + M[ 7] * m.M[14],
			M[ 4] * m.M[ 3] + M[ 5] * m.M[ 7] + M[ 6] * m.M[11] + M[ 7] * m.M[15],
			M[ 8] * m.M[ 0] + M[ 9] * m.M[ 4] + M[10] * m.M[ 8] + M[11] * m.M[12],
			M[ 8] * m.M[ 1] + M[ 9] * m.M[ 5] + M[10] * m.M[ 9] + M[11] * m.M[13],
			M[ 8] * m.M[ 2] + M[ 9] * m.M[ 6] + M[10] * m.M[10] + M[11] * m.M[14],
			M[ 8] * m.M[ 3] + M[ 9] * m.M[ 7] + M[10] * m.M[11] + M[11] * m.M[15],
			M[12] * m.M[ 0] + M[13] * m.M[ 4] + M[14] * m.M[ 8] + M[15] * m.M[12],
			M[12] * m.M[ 1] + M[13] * m.M[ 5] + M[14] * m.M[ 9] + M[15] * m.M[13],
			M[12] * m.M[ 2] + M[13] * m.M[ 6] + M[14] * m.M[10] + M[15] * m.M[14],
			M[12] * m.M[ 3] + M[13] * m.M[ 7] + M[14] * m.M[11] + M[15] * m.M[15]);
	}

	vec4<T> operator*(const vec4<T>& v) const {
		return vec4<T>(
			M[ 0] * v.x + M[ 1] * v.y + M[ 2] * v.z + M[ 3] * v.w,
			M[ 4] * v.x + M[ 5] * v.y + M[ 6] * v.z + M[ 7] * v.w,
			M[ 8] * v.x + M[ 9] * v.y + M[10] * v.z + M[11] * v.w,
			M[12] * v.x + M[13] * v.y + M[14] * v.z + M[15] * v.w);
	}

	vec4<T> operator/(const vec4<T>& v) const { return inv() * v; }

	T det() const {
		return
			M[ 3] * M[ 6] * M[ 9] * M[12] - M[ 2] * M[ 7] * M[ 9] * M[12] - M[ 3] * M[ 5] * M[10] * M[12] + M[ 1] * M[ 7] * M[10] * M[12]+
			M[ 2] * M[ 5] * M[11] * M[12] - M[ 1] * M[ 6] * M[11] * M[12] - M[ 3] * M[ 6] * M[ 8] * M[13] + M[ 2] * M[ 7] * M[ 8] * M[13]+
			M[ 3] * M[ 4] * M[10] * M[13] - M[ 0] * M[ 7] * M[10] * M[13] - M[ 2] * M[ 4] * M[11] * M[13] + M[ 0] * M[ 6] * M[11] * M[13]+
			M[ 3] * M[ 5] * M[ 8] * M[14] - M[ 1] * M[ 7] * M[ 8] * M[14] - M[ 3] * M[ 4] * M[ 9] * M[14] + M[ 0] * M[ 7] * M[ 9] * M[14]+
			M[ 1] * M[ 4] * M[11] * M[14] - M[ 0] * M[ 5] * M[11] * M[14] - M[ 2] * M[ 5] * M[ 8] * M[15] + M[ 1] * M[ 6] * M[ 8] * M[15]+
			M[ 2] * M[ 4] * M[ 9] * M[15] - M[ 0] * M[ 6] * M[ 9] * M[15] - M[ 1] * M[ 4] * M[10] * M[15] + M[ 0] * M[ 5] * M[10] * M[15];
	}

	mat4 inv() const {
		return mat4(
			M[ 6] * M[11] * M[13] - M[ 7] * M[10] * M[13] + M[ 7] * M[ 9] * M[14] - M[ 5] * M[11] * M[14] - M[ 6] * M[ 9] * M[15] + M[ 5] * M[10] * M[15],
			M[ 3] * M[10] * M[13] - M[ 2] * M[11] * M[13] - M[ 3] * M[ 9] * M[14] + M[ 1] * M[11] * M[14] + M[ 2] * M[ 9] * M[15] - M[ 1] * M[10] * M[15],
			M[ 2] * M[ 7] * M[13] - M[ 3] * M[ 6] * M[13] + M[ 3] * M[ 5] * M[14] - M[ 1] * M[ 7] * M[14] - M[ 2] * M[ 5] * M[15] + M[ 1] * M[ 6] * M[15],
			M[ 3] * M[ 6] * M[ 9] - M[ 2] * M[ 7] * M[ 9] - M[ 3] * M[ 5] * M[10] + M[ 1] * M[ 7] * M[10] + M[ 2] * M[ 5] * M[11] - M[ 1] * M[ 6] * M[11],
			M[ 7] * M[10] * M[12] - M[ 6] * M[11] * M[12] - M[ 7] * M[ 8] * M[14] + M[ 4] * M[11] * M[14] + M[ 6] * M[ 8] * M[15] - M[ 4] * M[10] * M[15],
			M[ 2] * M[11] * M[12] - M[ 3] * M[10] * M[12] + M[ 3] * M[ 8] * M[14] - M[ 0] * M[11] * M[14] - M[ 2] * M[ 8] * M[15] + M[ 0] * M[10] * M[15],
			M[ 3] * M[ 6] * M[12] - M[ 2] * M[ 7] * M[12] - M[ 3] * M[ 4] * M[14] + M[ 0] * M[ 7] * M[14] + M[ 2] * M[ 4] * M[15] - M[ 0] * M[ 6] * M[15],
			M[ 2] * M[ 7] * M[ 8] - M[ 3] * M[ 6] * M[ 8] + M[ 3] * M[ 4] * M[10] - M[ 0] * M[ 7] * M[10] - M[ 2] * M[ 4] * M[11] + M[ 0] * M[ 6] * M[11],
			M[ 5] * M[11] * M[12] - M[ 7] * M[ 9] * M[12] + M[ 7] * M[ 8] * M[13] - M[ 4] * M[11] * M[13] - M[ 5] * M[ 8] * M[15] + M[ 4] * M[ 9] * M[15],
			M[ 3] * M[ 9] * M[12] - M[ 1] * M[11] * M[12] - M[ 3] * M[ 8] * M[13] + M[ 0] * M[11] * M[13] + M[ 1] * M[ 8] * M[15] - M[ 0] * M[ 9] * M[15],
			M[ 1] * M[ 7] * M[12] - M[ 3] * M[ 5] * M[12] + M[ 3] * M[ 4] * M[13] - M[ 0] * M[ 7] * M[13] - M[ 1] * M[ 4] * M[15] + M[ 0] * M[ 5] * M[15],
			M[ 3] * M[ 5] * M[ 8] - M[ 1] * M[ 7] * M[ 8] - M[ 3] * M[ 4] * M[ 9] + M[ 0] * M[ 7] * M[ 9] + M[ 1] * M[ 4] * M[11] - M[ 0] * M[ 5] * M[11],
			M[ 6] * M[ 9] * M[12] - M[ 5] * M[10] * M[12] - M[ 6] * M[ 8] * M[13] + M[ 4] * M[10] * M[13] + M[ 5] * M[ 8] * M[14] - M[ 4] * M[ 9] * M[14],
			M[ 1] * M[10] * M[12] - M[ 2] * M[ 9] * M[12] + M[ 2] * M[ 8] * M[13] - M[ 0] * M[10] * M[13] - M[ 1] * M[ 8] * M[14] + M[ 0] * M[ 9] * M[14],
			M[ 2] * M[ 5] * M[12] - M[ 1] * M[ 6] * M[12] - M[ 2] * M[ 4] * M[13] + M[ 0] * M[ 6] * M[13] + M[ 1] * M[ 4] * M[14] - M[ 0] * M[ 5] * M[14],
			M[ 1] * M[ 6] * M[ 8] - M[ 2] * M[ 5] * M[ 8] + M[ 2] * M[ 4] * M[ 9] - M[ 0] * M[ 6] * M[ 9] - M[ 1] * M[ 4] * M[10] + M[ 0] * M[ 5] * M[10]) / det();
	}

	mat4 transpose() const {
		return mat4(
			M[ 0], M[ 4], M[ 8], M[12],
			M[ 1], M[ 5], M[ 9], M[13],
			M[ 2], M[ 6], M[10], M[14],
			M[ 3], M[ 7], M[11], M[15]);
	}

	T tr() const { return M[0] + M[5] + M[10] + M[15]; }
	vec4<T> solve(const vec4<T>& b) const { return inv() * b; }

	vec4<T> col1() const { return vec4<T>(M[ 0], M[ 4], M[ 8], M[12]); }
	vec4<T> col2() const { return vec4<T>(M[ 1], M[ 5], M[ 9], M[13]); }
	vec4<T> col3() const { return vec4<T>(M[ 2], M[ 6], M[10], M[14]); }
	vec4<T> col4() const { return vec4<T>(M[ 3], M[ 7], M[11], M[15]); }
	vec4<T> row1() const { return vec4<T>(M[ 0], M[ 1], M[ 2], M[ 3]); }
	vec4<T> row2() const { return vec4<T>(M[ 4], M[ 5], M[ 6], M[ 7]); }
	vec4<T> row3() const { return vec4<T>(M[ 8], M[ 9], M[10], M[11]); }
	vec4<T> row4() const { return vec4<T>(M[12], M[13], M[14], M[15]); }



	static mat4 rand() {
		return mat4(
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
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

	static mat4 rotX(T phi) {
		T s = sin(phi);
		T c = cos(phi);
		return mat4(
			1, 0, 0, 0,
			0, c, -s, 0,
			0, s, c, 0,
			0, 0, 0, 1);
	}

	static mat4 rotY(T phi) {
		T s = sin(phi);
		T c = cos(phi);
		return mat4(
			c, 0, s, 0,
			0, 1, 0, 0,
			-s, 0, c, 0,
			0, 0, 0, 1);
	}

	static mat4 rotZ(T phi) {
		T s = sin(phi);
		T c = cos(phi);
		return mat4(
			c, -s, 0, 0,
			s, c, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1);
	}

	static mat4 rotW(T phi) {
		T s = sin(phi);
		T c = cos(phi);
		return mat4(
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, c, s,
			0, 0, -s, c);
	}

	static mat4 rot(const vec4<T>& rot) {
		return rotX(rot.x) * rotY(rot.y) * rotZ(rot.z) * rotW(rot.w);
	}

	static mat4 rot(const vec3<T>& rot) {
		return rot(rot.x, rot.y, rot.z, 0);
	}

	static mat4 scale(const vec4<T>& v) {
		return mat4(
			v.x, 0, 0, 0,
			0, v.y, 0, 0,
			0, 0, v.z, 0,
			0, 0, 0, v.z);
	}

	static mat4 scale(const vec3<T>& v) {
		return scale(vec4<T>(v.x, v.y, v.z, 0));
	}

	static mat4 translation(const vec3<T>& t) {
		return mat4(
			1, 0, 0, t.x,
			0, 1, 0, t.y,
			0, 0, 1, t.z,
			0, 0, 0, 1);
	}

	static mat4 projection(const T d) {
		return mat4(
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, d, 0);
	}

	static mat4 perspective(float fovy, float aspect, float near, float far) {
		float a = tan(fovy * 0.5f);
		float b = a / aspect;
		float c = (near + far) / (near - far);
		float d = (2 * near * far) / (near - far);
		return mat4(
			b, 0, 0, 0,
			0, a, 0, 0,
			0, 0, c, d,
			0, 0, -1, 0);
	}

	static mat4 fromRows(const vec4<T>& v1, const vec4<T>& v2, const vec4<T>& v3, const vec4<T>& v4) {
		return mat4(
			v1.x, v1.y, v1.z, v1.w,
			v2.x, v2.y, v2.z, v2.w,
			v3.x, v3.y, v3.z, v3.w,
			v4.x, v4.y, v4.z, v4.w);
	}

	static mat4 fromColumns(const vec4<T>& v1, const vec4<T>& v2, const vec4<T>& v3, const vec4<T>& v4) {
		return mat4(
			v1.x, v2.x, v3.x, v4.x,
			v1.y, v2.y, v3.y, v4.y,
			v1.z, v2.z, v3.z, v4.z,
			v1.w, v2.w, v3.w, v4.w);
	}

};

typedef mat4<float> m4;

template<typename T>
const mat4<T> mat4<T>::I(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
template<typename T>
const mat4<T> mat4<T>::E(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);

}

#endif