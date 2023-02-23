#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "types.h"
#include "material.h"

class Sphere {
public:
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

#endif //__GEOMETRY_H__