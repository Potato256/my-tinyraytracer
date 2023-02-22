#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "geometry.h"
#include "types.h"
#include "material.h"
#include "light.h"

bool scene_intersect(
    const Vec3f &orig, 
    const Vec3f &dir, 
    const std::vector<Sphere> &spheres,
    Vec3f &hit,
    Vec3f &normal,
    Material& mat
) {
    float sphere_dist = std::numeric_limits<float>::max();
    for (int i = 0; i < spheres.size(); ++i) {
        float L_dist = std::numeric_limits<float>::max();
        if (spheres[i].ray_intersact(orig, dir, L_dist) &&
             L_dist < sphere_dist) {
            //std::cout<<"hit\n";
            sphere_dist = L_dist; 
            hit = orig + dir * sphere_dist;
            normal = (hit - spheres[i].center).normalize();
            mat = spheres[i].mat;
        }
    }
    return sphere_dist < 1000;
}

Vec3f cast_ray(
    const Vec3f &orig,
    const Vec3f &dir,
    const std::vector<Sphere> spheres
) {
    Vec3f hit, normal;
    Material mat;
    if (!scene_intersect(orig, dir, spheres, hit, normal, mat)) {
        // background color
        return Vec3f(0.2f, 0.7f, 0.8f);
    }
    // "shader" program
    return mat.diffuse_color;
}

void render(const std::vector<Sphere> &spheres) {
    const int width     = 1024;
    const int height    = 768;
    Vec3f cam_pos = Vec3f(0.0f, 0.0f, 0.0f);
    Vec3f cam_view = Vec3f(0.0f, 0.0f, -1.0f);
    float cam_to_screen = 1.0f;
    float hfov = M_PI/2.;

    std::vector<Vec3f> framebuffer(width*height);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            float x = (2*(j+0.5)/width - 1) * tan(hfov/2) * 
                width/float(height);
            float y = (2*(i+0.5)/height - 1) * tan(hfov/2);
            Vec3f dir = Vec3f(x, y, 0) + cam_view;
            dir.normalize();
            framebuffer[i*width+j] = 
                cast_ray(cam_pos, dir, spheres);
        }
    }

    // For more file formats, it's recommended to use a third-party libary, 
    // like stb
    std::ofstream ofs;
    ofs.open("./out.ppm", std::ofstream::out | std::ofstream::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (int i = 0; i < height*width; ++i) {
        for(int j = 0; j < 3; ++j) {
            ofs << (uint8_t) (255 * std::clamp(framebuffer[i][j], 0.0f, 1.0f));
        }
    }
    ofs.close();
}

int main() {
    std::vector<Sphere> spheres;

    Material      ivory(Vec3f(0.4, 0.4, 0.3));
    Material red_rubber(Vec3f(0.3, 0.1, 0.1));

    spheres.push_back(Sphere(Vec3f(-3,    0,   -16), 2,      ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2, red_rubber));
    spheres.push_back(Sphere(Vec3f( 1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(Vec3f( 7,    5,   -18), 4,      ivory));
    
    render(spheres);
    return 0;
}