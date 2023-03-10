#ifndef __TYPES_H__
#define __TYPES_H__

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

template <int32_t DIM, typename T> 
class vec {
public:
    vec() {
        for (int32_t i = 0; i < DIM; ++i)
            data_[i] = T();
    }
    T& operator[](const int32_t i) {assert(i<DIM&&i>=0); return data_[i];}
    const T& operator[](const int32_t i) const {assert(i<DIM&&i>=0); return data_[i];}
private:
    T data_[DIM];
};

// operator - : v-v
template <int32_t DIM, typename T> 
vec<DIM,T> operator-(vec<DIM, T> lhs, const vec<DIM, T>& rhs) {
    for (int32_t i=0; i < DIM; ++i) {lhs[i]-=rhs[i];};
    return lhs;
}
// operator - : -v
template <int32_t DIM, typename T> 
vec<DIM,T> operator-(vec<DIM, T> lhs) {
    for (int32_t i=0; i < DIM; ++i) {lhs[i]=-lhs[i];};
    return lhs;
}
// operator + : v+v
template <int32_t DIM, typename T> 
vec<DIM,T> operator+(vec<DIM, T> lhs, const vec<DIM, T>& rhs) {
    for (int32_t i=0; i < DIM; ++i) {lhs[i]+=rhs[i];};
    return lhs;
}
// operator += : v+=v
template <int32_t DIM, typename T> 
vec<DIM,T>& operator+=(vec<DIM, T>& lhs, const vec<DIM, T>& rhs) {
    for (int32_t i=0; i < DIM; ++i) {lhs[i]+=rhs[i];};
    return lhs;
}
// operator += : v+=k
template <int32_t DIM, typename T, typename U> 
vec<DIM,T>& operator+=(vec<DIM, T>& lhs, const U &rhs) {
    for (int32_t i=0; i < DIM; ++i) {lhs[i]+=rhs;};
    return lhs;
}
// operator * : v*v
template <int32_t DIM, typename T> 
T operator*(const vec<DIM, T> lhs, const vec<DIM, T>& rhs) {
    T ret = T();
    for (int32_t i=0; i < DIM; ++i) {ret += lhs[i]*rhs[i];};
    return ret;
}
// operator * : v*k
template<int32_t DIM,typename T,typename U> 
vec<DIM,T> operator*(const vec<DIM,T> &lhs, const U& rhs) {
    vec<DIM,T> ret;
    for (int32_t i=DIM; i--; ret[i]=lhs[i]*rhs);
    return ret;
}
// operator * : k*v
template<int32_t DIM,typename T,typename U> 
vec<DIM,T> operator*(const U& lhs, const vec<DIM,T> &rhs) {
    vec<DIM,T> ret;
    for (int32_t i=DIM; i--; ret[i]=rhs[i]*lhs);
    return ret;
}
// operator / : v/k
template<int32_t DIM,typename T,typename U> 
vec<DIM,T> operator/(const vec<DIM,T> &lhs, const U& rhs) {
    assert(rhs!=0);
    vec<DIM,T> ret;
    for (int32_t i=DIM; i--; ret[i]=lhs[i]/rhs);
    return ret;
}

template <int32_t DIM, typename T> 
std::ostream& operator<<(std::ostream& out, const vec<DIM,T>& v) {
    for(unsigned int i=0; i<DIM; i++) {
        out << v[i] << " " ;
    }
    return out ;
}

typedef vec<2, float> Vec2f;
typedef vec<3, float> Vec3f;
typedef vec<3, float> Vec3i;
typedef vec<4, float> Vec4f;

template <typename T> class vec<2,T> {
public:
    vec() : x(T()), y(T()) {}
    vec(T x_, T y_) :x(x_), y(y_) {}
    // template <class U> vec<2,T>(const vec<2,U> &v);
    T& operator[](const int32_t i) { assert(i<2&&i>=0); return i<=0 ? x : y; }
    const T& operator[](const int32_t i) const { assert(i<2&&i>=0); return i<=0 ? x : y; }
    T x, y;
};

template <typename T> class vec<3,T> {
public:
    vec() : x(T()), y(T()), z(T()) {}
    vec(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}
    T& operator[](const int32_t i) { assert(i>=0&&i<3); return i<=0 ? x : (1==i ? y : z); }
    const T& operator[](const int32_t i) const { assert(i>=0&&i<3); return i<=0 ? x : (1==i ? y : z); }
    float norm() { return std::sqrt(x*x+y*y+z*z); }
    vec<3,T> & normalize(T l=1) { *this = (*this)*(l/norm()); return *this; }
    T x,y,z;
};

template <typename T> 
vec<3,T> cross(vec<3,T> v1, vec<3,T> v2) {
    return vec<3,T>(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}


template <typename T> class vec<4,T> {
public:
    vec() : x(T()), y(T()), z(T()), w(T()) {}
    vec(T x_, T y_, T z_, T w_) : x(x_), y(y_), z(z_), w(w_) {}
    T& operator[](const int32_t i) { assert(i>=0&&i<4); return i<=1 ? (i==0?x:y) : (i==2?z:w); }
    const T& operator[](const int32_t i) const { assert(i>=0&&i<4); return i<=1 ? (i==0?x:y) : (i==2?z:w); }
    float norm() { return std::sqrt(x*x+y*y+z*z+w*w); }
    vec<3,T> & normalize(T l=1) { *this = (*this)*(l/norm()); return *this; }
    T x,y,z,w;
};

#endif //__TYPES_H__