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

const int kMaxRayDepth = 3;

typedef Vector3<float> Vec3f;

struct Intersection
{
    Sphere<float> * hitObject = nullptr;
    float tNear = Constants::kInfinity;
};

struct Scene
{
    std::vector<Sphere<float> *> objects;
    std::vector<Sphere<float> *> lights;
};

std::default_random_engine rng;
std::uniform_real_distribution<float> uniform(0, 1);

// g++ -Ofast -std=c++11 RayTracer.cpp Ray.cpp Vector3.cpp Sphere.cpp Image.cpp -o RayTracer -lomp; ./RayTracer; convert rt_out.ppm rt_out.png; open rt_out.png

float fresnelMix(const float a, const float b, const float m)
{
    return b * m + a * (1 - m);
}

void createSamplingCoords(const Vec3f & n, Vec3f & nx, Vec3f & nz)
{
    // Use surface hitNormal n to create a right handed coordinate system spanned by
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
    float sinTheta = std::sqrt(1. - r1 * r1); // r1 = cosTheta
    float phi = 2 * M_PI * r2;

    float x = sinTheta * std::cos(phi);
    float z = sinTheta * std::sin(phi);

    return Vec3f(x, r1, z);
}

Vec3f sampleHemisphere(void)
{
    return sampleHemisphere(uniform(rng), uniform(rng));
}

Intersection trace(Vec3f & origin, Vec3f & direction, Scene & scene)
{
    Intersection isect;

    for (auto obj : scene.objects) {

        float tNear = Constants::kInfinity;
        float tFar = Constants::kInfinity;

        // if (obj->intersect(origin, direction, &tNear, &tFar))
        // {
        //     std::cout << tNear << " " << isect.tNear << std::endl;
        // }

        bool isIntersected = obj->intersect(origin, direction, tNear, tFar);

        if (tNear > tFar) std::swap(tNear, tFar);

        if (tNear < 0) {
            tNear = tFar;
            if (tNear < 0) {
                isIntersected = false;
            }
        }

        if (isIntersected && tNear < isect.tNear) {
            isect.hitObject = obj;
            isect.tNear = tNear;
        }
    }

    return isect;
}

Vec3f cast(Vec3f & origin, Vec3f & direction, Scene & scene, int depth)
{
    if (depth > kMaxRayDepth) return 0;

    Vec3f hitColour;

    Intersection isect = trace(origin, direction, scene);

    if (isect.hitObject != nullptr) {

        // std::cout << isect.hitObject << std::endl;

        Vec3f hitPoint = origin + direction * isect.tNear;
        Vec3f hitNormal = isect.hitObject->hitNormal(hitPoint);

        Vec3f directLight;

        for (auto light : scene.lights) {

            Intersection isectShadow;
            Vec3f lightDirection;
            Vec3f lightIntensity;
            float distance;

            light->illuminate(hitPoint, lightDirection, lightIntensity, isectShadow.tNear);

            Vec3f shadOrigin = hitPoint + hitNormal + Constants::kEpsilon;
            Vec3f shadDirection = -lightDirection;

            Intersection block = trace(shadOrigin, shadDirection, scene);
            bool visible = block.hitObject == nullptr || block.tNear > isectShadow.tNear;

            directLight += visible * lightIntensity * std::max(float(0.), Utilities::dot(hitNormal, -lightDirection));
        }

        directLight /= float(M_PI);

        // Vec3f nx, ny;
        //
        // createSamplingCoords(hitNormal, nx, ny);
        //
        // int nSamples = 1000;

        Vec3f indirectLight;

        // for (int i = 0; i < nSamples; i++) {
        //
        //     Vec3f hSamp = sampleHemisphere();
        //
        //     float sx = hSamp.x * nx.x + hSamp.y * hitNormal.x + hSamp.z * ny.x;
        //     float sy = hSamp.x * nx.y + hSamp.y * hitNormal.y + hSamp.z * ny.y;
        //     float sz = hSamp.x * nx.z + hSamp.y * hitNormal.z + hSamp.z * ny.z;
        //
        //     Vec3f hSampWorld = Vec3f(sx, sy, sz);
        //     Vec3f newOrigin = hitPoint + hSampWorld * Constants::kEpsilon;
        //
        //     float cosTheta = hSamp.y;
        //
        //     indirectLight += cosTheta * cast(newOrigin, hSampWorld, scene, depth + 1);
        //
        // }
        //
        // indirectLight /= nSamples;

        hitColour = (directLight + 2. * indirectLight) * isect.hitObject->surfaceColor;

    }

    return hitColour;
}

void render(Image & img, Scene & scene)
{

    double invWidth = 1. / img.width;
    double invHeight = 1. / img.height;

    double fov = 90;

    double aspectRatio = img.width / img.height;

    double angle = tan(M_PI * 0.5 * fov / 180.);

    std::mutex imgMutex;

    // random_device rd;
    // mt19937 mt(rd());
    // normal_distribution<double> nx(0, 1.0);
    // normal_distribution<double> ny(0, 1.0);

    #pragma omp parallel for
    for (int y = 0; y < img.height; ++y) {
        std::cout << y << std::endl;
        for (int x = 0; x < img.width; ++x) {

            Vector3<int> total(0);

            double rayX = (2. * ((x + 0.5) * invWidth) - 1.) * angle * aspectRatio;
            double rayY = (1. - 2. * ((y + 0.5) * invHeight)) * angle;

            Vec3f rayDirection = Vec3f(rayX, rayY, 1); // -> camera faces z direction
            Vec3f rayOrigin = 0;

            Vec3f pixelColor = cast(rayOrigin, rayDirection, scene, 0);

            total += Vector3<int>(min(255, int(pixelColor.x * 255)),
                                  min(255, int(pixelColor.y * 255)),
                                  min(255, int(pixelColor.z * 255)));

            imgMutex.lock();
            img.pixels[y][x] = total;
            imgMutex.unlock();

            // Vector3<int> total(0);
            //
            // int samples = 100;
            //
            // for (int s = 0; s < samples; s++) {
            //     double rayX = (2. * ((x + nx(mt) + 0.5) * invWidth) - 1.) * angle * aspectRatio;
            //     double rayY = (1. - 2. * ((y + ny(mt) + 0.5) * invHeight)) * angle;
            //
            //     Vec3f rayDirection = Vec3f(rayX, rayY, 1); // -> camera faces z direction
            //     Vec3f rayOrigin = 0;
            //
            //     Vec3f pixelColor = cast(rayOrigin, rayDirection, scene, 0);
            //
            //     total += Vector3<int>(min(255, int(pixelColor.x * 255)),
            //                           min(255, int(pixelColor.y * 255)),
            //                           min(255, int(pixelColor.z * 255)));
            // }
            //
            // total.x = total.x / samples;
            // total.y = total.y / samples;
            // total.z = total.z / samples;
            //
            // imgMutex.lock();
            // img.pixels[y][x] = total;
            // imgMutex.unlock();

        }
    }
}

int main(int argc, char const *argv[])
{
    omp_set_num_threads(4);

    Image img(1024, 1024);
    // Image img(128, 128);

    Scene scene;

    // scene.objects.push_back(new Sphere<float>(Vec3f(0, 1, 10), 0.1, Vec3f(0, 1, 0), 0.0, 0.75));
    scene.objects.push_back(new Sphere<float>(Vec3f(0, 2.0, 4.), 1.0, Vec3f(0, 0, 1.0), 0.0, 1.0));
    scene.objects.push_back(new Sphere<float>(Vec3f(0, -10003, 30), 10000, Vec3f(1), 0.0, 1.0));

    scene.lights.push_back(new Sphere<float>(Vec3f(0, 0, -10), 10, Vec3f(1), 0.0, 0.0, Vec3f(1.0)));
    // scene.objects.push_back(new Sphere<float>(Vec3f(5, 5, -150.), 5, Vec3f(0, 0, 1), 0.0, 0.75, Vec3f(1.0)));

    render(img, scene);

    img.write("rt_out.ppm");

    return 0;
}
