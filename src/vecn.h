#ifndef VECN_H_INCL
#define VECN_H_INCL

#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace la {

template<typename T>
struct vecn {

	T* M;
	int length;

	vecn() {
		length = 0;
		M = NULL;
	}
	vecn(int len) {
		length = len;
		M = new T[length];
	}
	vecn(int len, T a = 0) {
		length = len;
		M = new T[length];
		set(a);
	}
	vecn(const vecn& v) {
		length = v.length;
		M = new T[length];
		set(v);
	}

	~vecn() {
		delete [] M;
	}

	vecn& set(T a = 0) {
		for (int i = 0; i < length; i++) {
			M[i] = a;
		}
		return *this;
	}
	vecn& set(const vecn& v) {
		if (v.length == length) {
			for (int i = 0; i < length; i++) {
				M[i] = v.M[i];
			}
		} else {
			delete [] M;
			length = v.length;
			M = new T[length];
			set(v);
		}
		return *this;
	}

	std::string str() const {
		std::ostringstream ss;
		ss << "[ ";
		for (int i = 0; i < length; i++) {
			ss << M[i] << " ";
		}
		ss << "]";
		return ss.str();
	}

	friend std::ostream& operator<<(std::ostream& out, const vecn& v) { return out << v.str(); }

	vecn& operator=(const vecn& v) { return set(v); }
	vecn operator-() const { return *this * (-1); }

	vecn operator+(const vecn& v) const {
		if (length == v.length) {
			vecn n(length);
			for (int i = 0; i < length; i++) {
				n.M[i] = M[i] + v.M[i];
			}
			return vecn;
		} else {
			return vecn();
		}
	}
	vecn operator-(const vecn& v) const {
		if (length == v.length) {
			vecn n(length);
			for (int i = 0; i < length; i++) {
				n.M[i] = M[i] - v.M[i];
			}
			return vecn;
		} else {
			return vecn();
		}
	}
	vecn operator*(const vecn& v) const {
		if (length == v.length) {
			vecn n(length);
			for (int i = 0; i < length; i++) {
				n.M[i] = M[i] * v.M[i];
			}
			return vecn;
		} else {
			return vecn();
		}
	}
	vecn operator/(const vecn& v) const {
		if (length == v.length) {
			vecn n(length);
			for (int i = 0; i < length; i++) {
				n.M[i] = M[i] / v.M[i];
			}
			return vecn;
		} else {
			return vecn();
		}
	}

	vecn operator*(T a) const {
		if (length == v.length) {
			vecn n(length);
			for (int i = 0; i < length; i++) {
				n.M[i] = M[i] * a;
			}
			return vecn;
		} else {
			return vecn();
		}
	}
	vecn operator/(T a) const { return *this * (1 / a); }

	vecn& operator+=(const vecn& v) {
		int l = length < v.length ? length : v.length;
		for (int i = 0; i < l; i++) {
			M[i] += v.M[i];
		}
		return *this;
	}
	vecn& operator-=(const vecn& v) {
		int l = length < v.length ? length : v.length;
		for (int i = 0; i < l; i++) {
			M[i] -= v.M[i];
		}
		return *this;
	}
	vecn& operator*=(const vecn& v) {
		int l = length < v.length ? length : v.length;
		for (int i = 0; i < l; i++) {
			M[i] *= v.M[i];
		}
		return *this;
	}
	vecn& operator/=(const vecn& v) {
		int l = length < v.length ? length : v.length;
		for (int i = 0; i < l; i++) {
			M[i] /= v.M[i];
		}
		return *this;
	}

	vecn& operator*=(T a) {
		for (int i = 0; i < length; i++) {
			M[i] *= a;
		}
		return *this;
	}
	vecn& operator/=(T a) {
		for (int i = 0; i < length; i++) {
			M[i] /= a;
		}
		return *this;
	}

	bool operator>(const vecn& v) const { return x > v.x && y > v.y && z > v.z; }
	bool operator>=(const vecn& v) const { return x >= v.x && y >= v.y && z >= v.z; }
	bool operator<(const vecn& v) const { return x < v.x && y < v.y && z < v.z; }
	bool operator<=(const vecn& v) const { return x <= v.x && y <= v.y && z <= v.z; }

	T dot(const vecn& v) const { return x * v.x + y * v.y + z * v.z; }
	vecn cross(const vecn& v) const { return vecn(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }

	T lensq() const { return dot(*this); }
	T len() const { return sqrt(lensq()); }
	T distsq(const vecn& v) const { return (*this - v).lensq(); }
	T dist(const vecn& v) const { return sqrt(distsq(v)); }
	vecn norm() const { return *this / len(); }
	T angle(const vecn& v) const { return acos(dot(v) / (len() * v.len())); }

	T projlen(const vecn& v) const { return dot(v) / v.len(); }
	vecn proj(const vecn& v) const { return v * (dot(v) / v.lensq()); }
	vecn reflect(const vecn& v) const { return proj(v) * 2 - *this; }
	vecn lerp(const vecn& v, T t) const { return *this + (v - *this) * t; }

	vecn refract(const vecn& n, T r) {
		T dotn = dot(n);
		T d = 1 - r * r * (1 - dotn * dotn);
		if (d < 0) {
			return vecn();
		} else {
			return *this * r - n * (r * dotn + sqrt(d));
		}
	}



	static vecn rand() {
		return vecn(
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX),
			(T)(::rand() / (float) RAND_MAX));
	}

	static vecn randunit() {
		float phi = ::rand() / (float) RAND_MAX;
		float theta = ::rand() / (float) RAND_MAX;
		return vecn(
			cos(phi) * cos(theta),
			sin(phi) * cos(theta),
			sin(theta));
	}

};

}

#endif