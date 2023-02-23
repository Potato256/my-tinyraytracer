#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include "types.h"

class Material {
public:
    //Material(Vec3f &color): diffuse_color(color) {}
    Material(const Vec4f &a, const Vec3f &color, const float &e, const float &eta): 
        albedo(a), diffuse_color(color), specular_exponent(e), refract_index(eta) {}
    Material(): 
        albedo(1, 0, 0, 0), diffuse_color(), specular_exponent(), refract_index() {}
    // order: diffuse, specular, reflect, refract
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
    float refract_index;
};

#endif // __MATERIAL_H__