#ifndef CPLX_H_INCL
#define CPLX_H_INCL

#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "mat2.h"

namespace la {

template<typename T>
struct cplx {

	static const cplx I;

	T x, y;

	cplx(T a = 0) { set(a); }
	cplx(T X, T Y) { set(X, Y); }
	cplx(const cplx& c) { set(c); }
	cplx(const vec2<T>& v) { set(v); }

	cplx& set(T X, T Y) {
		x = X;
		y = Y;
		return *this;
	}

	cplx& set(T a = 0) { return set(a, 0); }
	cplx& set(const cplx& c) { return set(c.x, c.y); }
	cplx& set(const vec2<T>& v) { return set(v.x, v.z); }

	std::string str() const {
		std::ostringstream ss;
		ss << x << " + " << y << "i";
		return ss.str();
	}

	friend std::ostream& operator<<(std::ostream& out, const cplx& c) { return out << c.str(); }

	cplx& operator=(const cplx& c) { return set(c); }
	cplx operator-() const { return *this * (-1); }

	cplx operator+(const cplx& c) const { return cplx(x + c.x, y + c.y); }
	cplx operator-(const cplx& c) const { return cplx(x - c.x, y - c.y); }
	cplx operator*(const cplx& c) const { return cplx(x * c.x - y * c.y, x * c.y + y * c.x); }
	cplx operator/(const cplx& c) const { return *this * c.inv(); }

	cplx operator*(T a) const { return cplx(x * a, y * a); }
	cplx operator/(T a) const { return cplx(x / a, y / a); }

	cplx& operator+=(const cplx& c) { return set(*this + c); }
	cplx& operator-=(const cplx& c) { return set(*this - c); }
	cplx& operator*=(const cplx& c) { return set(*this * c); }
	cplx& operator/=(const cplx& c) { return set(*this / c); }

	cplx& operator*=(T a) { return set(*this * a); }
	cplx& operator/=(T a) { return set(*this / a); }

	cplx conj() const { return cplx(x, -y); }
	T lensq() const { return x * x + y * y; }
	T len() const { return sqrt(lensq()); }
	T distsq(const cplx& c) const { return (*this - c).lensq(); }
	T dist(const cplx& c) const { return sqrt(distsq(c)); }
	cplx norm() const { return *this / len(); }
	cplx inv() const { return conj() / lensq(); }
	T real() const { return x; }
	T imag() const { return y; }
	T arg() const { return atan2(y, x); }
	cplx polar() const { return cplx(len(), arg()); }
	cplx exp() const { return cplx(cos(y), sin(y)) * exp(x); }
	mat2<T> mat() const { return mat2<T>(x, -y, y, x); }



	static cplx rand() {
		return cplx(
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX));
	}

	static cplx randnorm() { return rand().norm(); }

	static cplx fromPolar(const vec2<T>& v) { return cplx(v.x * cos(v.y), v.y * sin(v.y)); }

};

typedef cplx<float> cplxf;

template<typename T>
const cplx<T> cplx<T>::I(0, 1);

}

#endif