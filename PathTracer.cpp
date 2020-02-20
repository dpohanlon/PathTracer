#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <random>
#include <mutex>

#include <omp.h>

#include "Vector3.hpp"
#include "Sphere.hpp"
#include "Ray.h"
#include "Vec.h"
#include "Point.h"
#include "Utilities.hpp"
#include "Constants.hpp"
#include "Image.h"

using namespace std;

const int kMaxRayDepth = 10;

typedef Vector3<float> Vec3f;

std::default_random_engine rng;
std::uniform_real_distribution<float> uniform(0, 1);

// g++ -Ofast -std=c++11 RayTracer.cpp Ray.cpp Vector3.cpp Sphere.cpp Image.cpp -o RayTracer -lomp; ./RayTracer; convert rt_out.ppm rt_out.png; open rt_out.png

float fresnelMix(const float a, const float b, const float m)
{
    return b * m + a * (1 - m);
}

// Construct Scene -> Place Objects
// Render -> For each pixel - construct primary ray - Trace
// Trace -> Test for intersections (all?)
// Shade -> Apply shader (color)

void createSamplingCoords(const Vec3f & n, Vec3f & nx, Vec3f & nz)
{
    // Use surface normal n to create a right handed coordinate system spanned by
    // Nx and Nz, where Nx is perpendicular to n

    if (std::fabs(n.x) > std::fabs(n.y)) {
        nx = Vec3f(n.z, 0, -n.x).normalize();
    } else {
        nx = Vec3f(0, -n.z, n.y).normalize();
    }

    nz = Utilities::cross(n, nx);
}

Vec3f sampleHemisphere(float r1, float r2)
{
    float sinTheta = std::sqrt(1. - r1 * r1);
    float phi = 2 * M_PI * r2;

    float x = sinTheta * std::cos(phi);
    float z = sinTheta * std::sin(phi);

    return Vec3f(x, r1, z);
}

Vec3f sampleHemisphere(void)
{
    return sampleHemisphere(uniform(rng), uniform(rng));
}

Vec3f cast(Vec3f & origin, Vec3f & direction, int depth)
{
    if (depth > kMaxRayDepth) return 0;

    Vec3f normal; // CHANGE ME

    Vec3f nx, ny;

    createSamplingCoords(normal, nx, ny);

    Vec3f hSamp = sampleHemisphere();

    float sx = hSamp.x * nx.x + hSamp.y * normal.x + hSamp.z * ny.x;
    float sy = hSamp.x * nx.y + hSamp.y * normal.y + hSamp.z * ny.y;
    float sz = hSamp.x * nx.z + hSamp.y * normal.z + hSamp.z * ny.z;

    Vec3f hSampWorld(sx, sy, sz);

    return 0;
}

// void render(Image & img, const vector<Sphere<double> > & obj, const vector<Sphere<double> > & lights)
// {
//
//     double invWidth = 1. / img.width;
//     double invHeight = 1. / img.height;
//
//     double fov = 30;
//
//     double aspectRatio = img.width / img.height;
//
//     double angle = tan(M_PI * 0.5 * fov / 180.);
//
//     int samples = 256;
//
//     random_device rd;
//     mt19937 mt(rd());
//     normal_distribution<double> nx(0, 1.0);
//     normal_distribution<double> ny(0, 1.0);
//
//     std::mutex imgMutex;
//
//     #pragma omp parallel for
//     for (int y = 0; y < img.height; ++y) {
//         for (int x = 0; x < img.width; ++x) {
//
//             Vector3<int> total(0);
//
//             for (int s = 0; s < samples; s++) {
//                 double rayX = (2. * ((x + nx(mt) + 0.5) * invWidth) - 1.) * angle * aspectRatio;
//                 double rayY = (1. - 2. * ((y + ny(mt) + 0.5) * invHeight)) * angle;
//
//                 Vec rayDir = Vec(rayX, rayY, 1); // -> camera faces z direction
//
//                 Ray ray = Ray(Point(0.), rayDir);
//                 Vector3<double> pixelColor = trace(ray, obj, lights, 0);
//
//                 total += Vector3<int>(min(255, int(pixelColor.x * 255)),
//                                       min(255, int(pixelColor.y * 255)),
//                                       min(255, int(pixelColor.z * 255)));
//             }
//
//             total.x = total.x / samples;
//             total.y = total.y / samples;
//             total.z = total.z / samples;
//
//             imgMutex.lock();
//             img.pixels[y][x] = total;
//             imgMutex.unlock();
//         }
//     }
// }

int main(int argc, char const *argv[])
{

    omp_set_num_threads(4);

    Image img(1024, 1024);
    // Image img(640, 640);

    std::vector<Sphere<float> > objects;
    std::vector<Sphere<float> > lights;

    objects.push_back(Sphere<float>(Vec3f(0, 1, 30), 3, Vec3f(0, 1, 0), 0.0, 0.75));
    objects.push_back(Sphere<float>(Vec3f(5, 3, 30), 4, Vec3f(0, 0, 1), 0.0, 0.75));
    objects.push_back(Sphere<float>(Vec3f(0, -10003, 30), 10000, Vec3f(1), 0.0, 0.8));

    objects.push_back(Sphere<float>(Vec3f(5, 12, -2), 10, Vec3f(1), 0.0, 0.0, Vec3f(0.15)));

    // render(img, objects, lights);

    img.write("rt_out.ppm");

    return 0;
}
