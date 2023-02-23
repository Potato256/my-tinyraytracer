#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include "types.h"

class Material {
public:
    //Material(Vec3f &color): diffuse_color(color) {}
    Material(const Vec2f &a, const Vec3f &color, const float &e): 
        albedo(a), diffuse_color(color), specular_exponent(e) {}
    Material(): albedo(1, 0), diffuse_color(), specular_exponent() {}
    Vec2f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
};

#endif // __MATERIAL_H__