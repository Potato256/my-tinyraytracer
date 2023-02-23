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
    const std::vector<Light> lights,
    const std::vector<Sphere> spheres
) {
    Vec3f hit, normal;
    Material mat;
    if (!scene_intersect(orig, dir, spheres, hit, normal, mat)) {
        // background color
        return Vec3f(0.2f, 0.7f, 0.8f);
    }
    // "shader" program
    float diffuse_light_intensity = 0.0f;
    float specular_light_intensity = 0.0f;
    for (int i = 0; i < lights.size(); ++i) {
        Vec3f light_dir = lights[i].position - hit;
        float light_d = std::max(0.1f, light_dir.norm());
        float attenuate = 1; //1 / (light_d * light_d);
        //std::cout<<light_d<<std::endl;
        light_dir.normalize();
        diffuse_light_intensity += 
            lights[i].intensity * std::max(0.0f, light_dir * normal) * attenuate;
        Vec3f half_vector = (light_dir + (-dir)) / 2.0f;
        half_vector.normalize();
      
        specular_light_intensity +=
            lights[i].intensity * powf(std::max(0.0f, half_vector * normal), 
                mat.specular_exponent) * attenuate;
    }
    return mat.diffuse_color * diffuse_light_intensity * mat.albedo[0] +
           Vec3f(1.0f, 1.0f, 1.0f) * specular_light_intensity * mat.albedo[1];
}

void render(
    const std::vector<Light> &lights,
    const std::vector<Sphere> &spheres
) {
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
            float y = -(2*(i+0.5)/height - 1) * tan(hfov/2);
            Vec3f dir = Vec3f(x, y, 0) + cam_view;
            dir.normalize();
            framebuffer[i*width+j] = 
                cast_ray(cam_pos, dir, lights, spheres);
        }
    }

    // For more file formats, it's recommended to use a third-party libary, 
    // like stb
    std::ofstream ofs;
    ofs.open("./out.ppm", std::ofstream::out | std::ofstream::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (int i = 0; i < height*width; ++i) {
        for(int j = 0; j < 3; ++j) {
            Vec3f &c = framebuffer[i];
            float max = std::max(c[0], std::max(c[1], c[2]));
            if (max>1) c = c*(1./max);
            ofs << (uint8_t) (255 * std::clamp(framebuffer[i][j], 0.0f, 1.0f));
        }
    }
    ofs.close();
}

int main() {
    std::vector<Sphere> spheres;
    Material      ivory(Vec2f(0.6,  0.3), Vec3f(0.4, 0.4, 0.3),   50.f);
    Material red_rubber(Vec2f(0.9,  0.1), Vec3f(0.3, 0.1, 0.1),   10.f);

    spheres.push_back(Sphere(Vec3f(-3,    0,   -16), 2,      ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2, red_rubber));
    spheres.push_back(Sphere(Vec3f( 1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(Vec3f( 7,    5,   -18), 4,      ivory));
    
    std::vector<Light> lights;
    lights.push_back(Light(Vec3f(-10, 10,  10), 1.50));
    lights.push_back(Light(Vec3f( 15, 20, -10), 1.80));
    lights.push_back(Light(Vec3f( 15, 10,  15), 1.70));

    render(lights, spheres);
    return 0;
}