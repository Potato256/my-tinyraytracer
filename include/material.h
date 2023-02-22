#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include "types.h"

class Material {
public:
    Material(Vec3f &color): diffuse_color(color) {}
    Material(const Vec3f &color): diffuse_color(color) {}
    Material(): diffuse_color() {}
    Vec3f diffuse_color;
};

#endif // __MATERIAL_H__