#ifndef MAT2_H_INCL
#define MAT2_H_INCL

#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "vec2.h"

namespace la {

template<typename T>
struct mat2 {

	static const mat2 I;
	static const mat2 E;

	T M[4];

	explicit mat2(T a = 0) { set(a); }
	mat2(T m11, T m12, T m21, T m22) { set(m11, m12, m21, m22); }
	mat2(const mat2& m) { set(m); }

	mat2& set(T a) { return set(a, a, a, a); }
	mat2& set(const mat2& m) { return set(m.M[0], m.M[1], m.M[2], m.M[3]); }
	
	mat2& set(T m11, T m12, T m21, T m22) {
		M[0] = m11;
		M[1] = m12;
		M[2] = m21;
		M[3] = m22;
		return *this;
	}

	std::string str() const {
		std::ostringstream ss;
		ss << "[ " << M[0] << ", " << M[1] << "; "
			<< M[2] << ", " << M[3] << " ]";
		return ss.str();
	}

	friend std::ostream& operator<<(std::ostream& out, const mat2& m) { return out << m.str(); }

	mat2& operator=(const mat2& m) { return set(m); }
	mat2 operator-() const { return *this * (-1); }

	mat2 operator+(const mat2& m) const { return mat2(M[0] + m.M[0], M[1] + m.M[1], M[2] + m.M[2], M[3] + m.M[3]); }
	mat2 operator-(const mat2& m) const { return mat2(M[0] - m.M[0], M[1] - m.M[1], M[2] - m.M[2], M[3] - m.M[3]); }

	mat2 operator*(T a) const { return mat2(M[0] * a, M[1] * a, M[2] * a, M[3] * a); }
	mat2 operator/(T a) const { return *this * (1 / a); }

	mat2 operator*(const mat2& m) const {
		return mat2(
			M[0] * m.M[0] + M[1] * m.M[2],
			M[0] * m.M[1] + M[1] * m.M[3],
			M[2] * m.M[0] + M[3] * m.M[2],
			M[2] * m.M[1] + M[3] * m.M[3]);
	}

	vec2<T> operator*(const vec2<T>& v) const {
		return vec2<T>(
			M[0] * v.x + M[1] * v.y,
			M[2] * v.x + M[3] * v.y);
	}

	vec2<T> operator/(const vec2<T>& v) const { return inv() * v; }

	T det() const { return M[0] * M[3] - M[2] * M[1]; }
	mat2 inv() const { return mat2(M[3], -M[1], -M[2], M[0]) / det(); }
	mat2 transpose() const { return mat2(M[0], M[2], M[1], M[3]); }
	T tr() const { return M[0] + M[3]; }
	vec2<T> solve(const vec2<T>& b) const { return inv() * b; }

	vec2<T> col1() const { return vec2<T>(M[0], M[2]); }
	vec2<T> col2() const { return vec2<T>(M[1], M[3]); }
	vec2<T> row1() const { return vec2<T>(M[0], M[1]); }
	vec2<T> row2() const { return vec2<T>(M[2], M[3]); }



	static mat2 rand() {
		return mat2(
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX));
	}

	static mat2 rot(T phi) {
		T s = sin(phi);
		T c = cos(phi);
		return mat2(
			c, -s,
			s, c);
	}

	static mat2 scale(const vec2<T>& v) {
		return mat2(
			v.x, 0,
			0, v.y);
	}

	static mat2 fromRows(const vec2<T>& v1, const vec2<T>& v2) {
		return mat2(
			v1.x, v1.y,
			v2.x, v2.y);
	}

	static mat2 fromColumns(const vec2<T>& v1, const vec2<T>& v2) {
		return mat2(
			v1.x, v2.x,
			v1.y, v2.y);
	}

};

typedef mat2<float> m2;

template<typename T>
const mat2<T> mat2<T>::I(1, 0, 0, 1);
template<typename T>
const mat2<T> mat2<T>::E(1, 1, 1, 1);

}

#endif