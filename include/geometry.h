#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "types.h"
#include "material.h"

struct Sphere {
    Vec3f center;
    float radius;
    Material mat;

    Sphere(const Vec3f &c, const float &r, const Material &m):
        center(c), radius(r), mat(m) {}
    bool ray_intersact(const Vec3f &orig, const Vec3f &dir, float &t0) const {
        Vec3f L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

// I is toward the surface
Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}

// I is toward the surface
Vec3f refract(const Vec3f &I, const Vec3f &N, const float &refractive_index) { 
    // Snell's law
    float cosi = - std::clamp(I*N, -1.f, 1.f);
    float etai = 1, etat = refractive_index;
    Vec3f n = N;
    if (cosi < 0) {
    // if the ray is inside the object, 
    // swap the indices and invert the normal to get the correct result
    cosi = -cosi;
        std::swap(etai, etat); n = -N;
    }
    float eta = etai / etat;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    //return k < 0 ? Vec3f(0,0,0) : I*eta + n*(eta * cosi - sqrtf(k));
    // consider total reflection
    return k < 0 ? I - N*2.f*(I*n) : I*eta + n*(eta * cosi - sqrtf(k));
}
#endif //__GEOMETRY_H__