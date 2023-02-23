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
    const std::vector<Sphere> spheres,
    int32_t depth=0
) {
    Vec3f hit, normal;
    Material mat;
    if (depth>=5 || !scene_intersect(orig, dir, spheres, hit, normal, mat)) {
        // background color
        return Vec3f(0.2f, 0.7f, 0.8f);
    }
    // "shader" program
    float diffuse_intensity = 0.0f;
    float specular_intensity = 0.0f;
    for (int i = 0; i < lights.size(); ++i) {
        Vec3f light_dir = lights[i].position - hit;
        // The back side will not be shaded
        if (light_dir * normal < 0)
            continue;
        //std::cout<<"1\n";
        float light_d = std::max(0.01f, light_dir.norm());
        float attenuate = 1; //1 / (light_d * light_d);
        //std::cout<<light_d<<std::endl;
        light_dir.normalize();
        // shadow ray
        Vec3f shadow_dir = light_dir;
        Vec3f shadow_hit, shadow_normal;
        Material shadow_mat;
        if(scene_intersect(hit+0.001*normal, shadow_dir, spheres, 
            shadow_hit, shadow_normal, shadow_mat)) {
            if((shadow_hit-hit).norm() < light_d)
                continue;
        }
        //std::cout<<"2\n";
        // diffuse shading
        diffuse_intensity += 
            lights[i].intensity * light_dir * normal * attenuate;
        Vec3f half_vector = (light_dir + (-dir)) / 2.0f;
        half_vector.normalize();
        // specular shading
        specular_intensity +=
            lights[i].intensity * powf(half_vector * normal, 
                mat.specular_exponent) * attenuate;
    }
    // reflect shading
    Vec3f reflect_orig = hit + 0.001 * normal;
    Vec3f reflect_dir = reflect(dir, normal).normalize();
    Vec3f reflect_intensity = cast_ray(reflect_orig, reflect_dir, 
        lights, spheres, depth+1);
    // sum up
    return mat.diffuse_color * diffuse_intensity * mat.albedo[0] +
           Vec3f(1.0f, 1.0f, 1.0f) * specular_intensity * mat.albedo[1] +
           reflect_intensity * mat.albedo[2];
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
    Material      ivory(Vec3f(0.6,  0.3, 0.1), Vec3f(0.4, 0.4, 0.3),   50.);
    Material red_rubber(Vec3f(0.9,  0.1, 0.0), Vec3f(0.3, 0.1, 0.1),   10.);
    Material     mirror(Vec3f(0.0, 10.0, 0.8), Vec3f(1.0, 1.0, 1.0), 1425.);

    spheres.push_back(Sphere(Vec3f(-3,    0,   -16), 2,      ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2,      mirror));
    spheres.push_back(Sphere(Vec3f( 1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(Vec3f( 7,    5,   -18), 4,     mirror));
    spheres.push_back(Sphere(Vec3f( 0,    8,   -20), 1,      ivory));
    
    std::vector<Light> lights;
    lights.push_back(Light(Vec3f(-10, 10,  10), 1.50));
    lights.push_back(Light(Vec3f( 15, 20, -10), 1.80));
    lights.push_back(Light(Vec3f( 15, 10,  15), 1.70));

    render(lights, spheres);
    return 0;
}