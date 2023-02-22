#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "geometry.h"
#include "types.h"

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const Sphere &sphere) {
    float sphere_dist = std::numeric_limits<float>::max();
    if (!sphere.ray_intersact(orig, dir, sphere_dist)) {
        return Vec3f(0.2f, 0.7f, 0.8f);
    }
    return Vec3f(0.4f, 0.4f, 0.3f);
}

void render(const Sphere sphere) {
    const int width     = 1024;
    const int height    = 768;
    Vec3f cam_pos = Vec3f(0.0f, 0.0f, 0.0f);
    Vec3f cam_view = Vec3f(0.0f, 0.0f, -1.0f);
    float cam_to_screen = 1.0f;
    float wfov = 90.f;

    std::vector<Vec3f> framebuffer(width*height);
    for (int32_t i = 0; i < height; ++i) {
        for (int32_t j = 0; j < width; ++j) {
            float x = (2*(j+0.5)/width - 1) * tan(wfov/2);
            float y = (2*(i+0.5)/height - 1) * tan(wfov/2) * 
                height/float(width);
            Vec3f dir = Vec3f(x, y, 0) + cam_view;
            dir.normalize();
            framebuffer[i*width+j] = 
                cast_ray(cam_pos, dir, sphere);
        }
    }

    // For more file formats, it's recommended to use a third-party libary, 
    // like stb
    std::ofstream ofs;
    ofs.open("./out.ppm", std::ofstream::out | std::ofstream::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (int32_t i = 0; i < height*width; ++i) {
        for(int32_t j = 0; j < 3; ++j) {
            ofs << (uint8_t) (255 * std::clamp(framebuffer[i][j], 0.0f, 1.0f));
        }
    }
    ofs.close();
}

int main() {
    std::vector<Sphere> spheres;
    Sphere sphere(Vec3f(-3, 0, -16), 2);
    render(sphere);
    return 0;
}