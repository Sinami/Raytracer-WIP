#ifndef VECTOR3_H
#define VECTOR3_H

#include <cmath>
#include <iostream>

using std::sqrt;

class vector3 {
    public:
    vector3() : e{0,0,0} {}
    vector3(double e0, double e1, double e2) : e{e0, e1, e2} {}

    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    vector3 operator-() const { return vector3(-e[0], -e[1], -e[2]); }
    double operator[](int i) const { return e[i]; }
    double& operator[](int i) { return e[i]; }

    vector3& operator+=(const vector3 &v) {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    vector3& operator*=(const double t) {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    vector3& operator/=(const double t) {
        return *this *= 1/t;
    }

    double length() const {
        return sqrt(length_squared());
    }

    double length_squared() const {
        return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
    }

    public:
    double e[3];
};

using point3 = vector3;
using color = vector3;



// vector3 Utiity Functions

inline std::ostream& operator<<(std::ostream &out, const vector3 &v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vector3 operator+(const vector3 &u, const vector3 &v){
    return vector3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vector3 operator-(const vector3 &u, const vector3 &v){
    return vector3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vector3 operator*(const vector3 &u, const vector3 &v){
    return vector3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vector3 operator*(double t, const vector3 &v) {
    return vector3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline vector3 operator*(const vector3 &v, double t) {
    return t * v;
}

inline vector3 operator/(vector3 v, double t) {
    return (1/t) * v;
}

inline double dot(const vector3 &u, const vector3 &v) {
    return u.e[0] * v.e[0]
        +  u.e[1] * v.e[1]
        +  u.e[2] * v.e[2];
}

inline vector3 cross(const vector3 &u, const vector3 &v) {
    return vector3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                   u.e[2] * v.e[0] - u.e[0] * v.e[2],
                   u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vector3 unit_vector(vector3 v) {
    return v / v.length();
}

#endif