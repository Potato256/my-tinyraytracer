#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "geometry.h"

void render() {
    const int width     = 1024;
    const int height    = 768;
    std::vector<Vec3f> framebuffer(width*height);
    for (int32_t i = 0; i < height; ++i) {
        for (int32_t j = 0; j < width; ++j) {
            framebuffer[i*width+j] = 
                Vec3f(i/float(height), j/float(width), 0);
        }
    }

    std::ofstream ofs;
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (int32_t i = 0; i < height*width; ++i) {
        for(int32_t j = 0; j < 3; ++j) {
            ofs << (uint8_t) (255 * std::clamp(framebuffer[i][j], 0.0f, 1.0f));
        }
    }
    ofs.close();
}

int main() {
    render();
    return 0;
}