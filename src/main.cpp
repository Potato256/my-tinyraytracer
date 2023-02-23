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

#define AA
#define AA4
//#define AA16

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

    float checkerboard_dist = std::numeric_limits<float>::max();
    if(fabs(dir.y)>1e-3) {
        float d = -(orig.y+6)/dir.y; // y=-6
        Vec3f pt = orig + dir*d;
        if (d>0 && fabs(pt.x)<15 && pt.z<-5 && pt.z>-35 && d<sphere_dist)
        {
            checkerboard_dist = d;
            hit = pt;
            normal = Vec3f(0,1,0);
            mat.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ?
                Vec3f(1,1,1) : Vec3f(1, .7, .3);
            mat.diffuse_color = mat.diffuse_color * .3;
        }
    }
    return std::min(checkerboard_dist, sphere_dist) < 1000;
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
    if (depth>=4 || !scene_intersect(orig, dir, spheres, hit, normal, mat)) {
        // background color
        return Vec3f(0.2f, 0.7f, 0.8f);
    }
    // "shader" program
    float diffuse_intensity = 0.0f;
    float specular_intensity = 0.0f;
    for (int i = 0; i < lights.size(); ++i) {
        Vec3f light_dir = lights[i].position - hit;

        // The back side will not be shaded
        // if (light_dir * normal < 0)
        //    continue;
        
        //std::cout<<"1\n";
        float light_d = std::max(0.01f, light_dir.norm());
        float attenuate = 1; //1 / (light_d * light_d);
        //std::cout<<light_d<<std::endl;
        light_dir.normalize();
        // shadow ray
        Vec3f shadow_orig = normal*dir>0 ? hit-0.001*normal : hit+0.001*normal;
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
            lights[i].intensity * std::max(0.0f, light_dir * normal) * attenuate;
        Vec3f half_vector = (light_dir + (-dir)) / 2.0f;
        half_vector.normalize();
        // specular shading
        specular_intensity +=
            lights[i].intensity * powf(std::max(0.0f, half_vector * normal), 
                mat.specular_exponent) * attenuate;
    }

    // reflect shading
    Vec3f reflect_intensity = Vec3f();
    if (mat.albedo[2] > 1e-3) {
        Vec3f reflect_dir = reflect(dir, normal).normalize();
        Vec3f reflect_orig = reflect_dir*normal>0 ? hit+normal*1e-3 : hit-normal*1e-3;
        reflect_intensity = cast_ray(reflect_orig, reflect_dir, 
            lights, spheres, depth+1);
    }
    
    // refract shading
    Vec3f refract_intensity = Vec3f();
    if (mat.albedo[3] > 1e-3) {
        Vec3f refract_dir = refract(dir, normal, mat.refract_index).normalize();
        Vec3f refract_orig = refract_dir*normal>0 ? hit+normal*1e-3 : hit-normal*1e-3;
        refract_intensity = cast_ray(refract_orig, refract_dir, 
            lights, spheres, depth+1);
    }

    // sum up
    return mat.diffuse_color * diffuse_intensity * mat.albedo[0] +
           Vec3f(1.0f, 1.0f, 1.0f) * specular_intensity * mat.albedo[1] +
           reflect_intensity * mat.albedo[2] +
           refract_intensity * mat.albedo[3];
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

    std::vector<Vec3f> framebuffer(width*height, Vec3f());
#ifndef AA 
#pragma omp parallel for
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            float x = (2*(j+0.5)/width - 1) * tan(hfov/2) * 
                width/float(height);
            float y = -(2*(i+0.5)/height - 1) * tan(hfov/2);
            Vec3f dir = Vec3f(x, y, 0) + cam_view;
            dir.normalize();
            framebuffer[i*width+j] += 
                cast_ray(cam_pos, dir, lights, spheres);
        }
    }
#endif // AA
#ifdef AA4
#pragma omp parallel for
    for (int i = 0; i < 2*height; ++i) {
        for (int j = 0; j < 2*width; ++j) {
            float x = (2*(j+0.5)/(width*2) - 1) * tan(hfov/2) * 
                width/float(height);
            float y = -(2*(i+0.5)/(height*2) - 1) * tan(hfov/2);
            Vec3f dir = Vec3f(x, y, 0) + cam_view;
            dir.normalize();
            framebuffer[(i/2)*width+(j/2)] += 
                cast_ray(cam_pos, dir, lights, spheres) / 4;
        }
    }
#endif // AA4
#ifdef AA16
#pragma omp parallel for
    for (int i = 0; i < 4*height; ++i) {
        for (int j = 0; j < 4*width; ++j) {
            float x = (2*(j+0.5)/(width*4) - 1) * tan(hfov/2) * 
                width/float(height);
            float y = -(2*(i+0.5)/(height*4) - 1) * tan(hfov/2);
            Vec3f dir = Vec3f(x, y, 0) + cam_view;
            dir.normalize();
            framebuffer[(i/4)*width+(j/4)] += 
                cast_ray(cam_pos, dir, lights, spheres) / 16;
        }
    }
#endif // AA16

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
    Material      ivory(Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3),   50. ,1.0);
    Material      glass(Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8),  125. ,1.5);
    Material red_rubber(Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1),   10. ,1.0);
    Material     mirror(Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425. ,1.0);

    spheres.push_back(Sphere(Vec3f(-3,    0,   -16), 2,      ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2,      glass));
    spheres.push_back(Sphere(Vec3f(3.0, -1.5, -12), 0.5,      glass));
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