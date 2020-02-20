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

#include "noise.h"
#include "noiseutils.h"

using namespace std;

const int kMaxRayDepth = 10;

// g++ -Ofast -std=c++11 RayTracer.cpp Ray.cpp Vector3.cpp Sphere.cpp Image.cpp -o RayTracer -lomp; ./RayTracer; convert rt_out.ppm rt_out.png; open rt_out.png 

// module::Perlin noiseModule;
// utils::NoiseMap heightMap;
// utils::NoiseMapBuilderPlane heightMapBuilder;

float fresnelMix(const float a, const float b, const float m)
{
    return b * m + a * (1 - m);
}

// Construct Scene -> Place Objects
// Render -> For each pixel - construct primary ray - Trace
// Trace -> Test for intersections (all?)
// Shade -> Apply shader (color)

Vector3<double> trace(Ray & r, const vector<Sphere<double> > & obj, const vector<Sphere<double> > & lights, int depth)
{
    // std::cout << "Depth " << depth << std::endl;
    // std::cout << std::endl;

    Vector3<double> minDist = Constants::kInfinity;

    // This will eventually be a pointer to the base class
    const Sphere<double> * nearObj = nullptr;

    for (auto & o : obj) {

        if (o.intersect(r)) {
            Vector3<double> thisDist = r(r.t);

            if (thisDist.z < minDist.z){
                minDist = thisDist;
                nearObj = &o;
            }
        }
    }

    // If there are no intersections, return the background colour
    if (!nearObj) return Vector3<double>(0.);

    // Point pHit = r(r.t); // Doesn't work w/intersection -> why?
    Point pHit = minDist;

    Vec nHit = pHit - nearObj->center; // normal to hit point, only for sphere!
    nHit.normalize();

    bool insideObj = false;
    if (Utilities::dot(r.direction, nHit) > 0) {
        nHit = -nHit;
        insideObj = true;
    }

    Vector3<double> surfaceColor = Vector3<double>(0.);

    // std::cout << nearObj->transparency  << ", " << nearObj->reflection << std::endl;

    if ((nearObj->transparency > Constants::kEpsilon ||
         nearObj->reflection > Constants::kEpsilon) &&
         depth < kMaxRayDepth) {

        float facingRatio = -Utilities::dot(r.direction, nHit);

        // float fresnel =  fresnelMix(pow(1 - facingRatio, 3), 1, 0.1); // Broken?
        float fresnel = 1.0;

        Vector3<double> reflDir = r.direction - nHit * 2 * Utilities::dot(r.direction, nHit);
        reflDir.normalize();

        Ray reflRay = Ray(pHit + nHit + Constants::kEpsilon, reflDir);

        Vector3<double> reflection = trace(reflRay, obj, lights, depth + 1);

        Vector3<double> refraction = Vector3<double>(0);

        if (nearObj->transparency > Constants::kEpsilon) {
            // std::cout << "transparent" << std::endl;
            float ior = 1.1;
            float eta = insideObj ? ior : 1 / ior;
            float cosi = Utilities::dot(r.direction, nHit);
            float k =  1. - eta * eta * (1. - cosi * cosi);

            Vector3<double> refracDir = r.direction * eta + nHit * (eta * cosi - sqrt(k));
            refracDir.normalize();

            Ray refracRay = Ray(pHit - nHit + Constants::kEpsilon, refracDir);

            Vector3<double> refraction = trace(refracRay, obj, lights, depth + 1);
        }

        surfaceColor = nearObj->surfaceColor * (nearObj->reflection * reflection * fresnel + refraction * (1 - fresnel) * nearObj->transparency);
        // surfaceColor = refraction;

        // std::cout << "surf " << surfaceColor << std::endl;
        // std::cout << reflection << std::endl;

    }
    // else {

        for (auto & oi : obj) {
            if (!oi.emissive) continue;
            // std::cout << "EMISSIVE" << std::endl;
            Ray lightRay = Ray(Point(oi.center), Vec(oi.center - pHit));

            bool blocked = false;

            // std::cout << blocked << std::endl;

            for (auto & oj : obj) {
                if (oj.is(oj)) continue;
                if (oj.intersect(lightRay)) {
                    // std::cout << "Blocked by obj at " << oj.center << std::endl;
                    // Attenuate according to transparency
                    // In this case just zero the light
                    blocked = true;
                }
            }

            if (!blocked) {
                // std::cout << "Not blocked" << std::endl;
                // random_device rd;
                // mt19937 mt(rd());
                // normal_distribution<double> n(0, 0.05);

                // nHit += Vec(n(mt), n(mt), 0);
                // nHit.normalize();

                // nHit += 0.2 * Vec(0, 0, heightMap.GetValue(x, y));

                // nHit.normalize();

                double lambert = std::max(0., Utilities::dot(lightRay.direction, nHit));
                // std::cout << "Lambert: " << lambert << std::endl;
                Vec phongDir = 2 * std::max(0., Utilities::dot(lightRay.direction, nHit)) * nHit - lightRay.direction;
                phongDir.normalize();

                double i_s = 0.8;
                double k_s = 0.2;
                double alpha = 10.;

                double phong = pow(std::max(0., Utilities::dot(phongDir, Vec(0, 0, -1))), alpha) * i_s * k_s;

                // Missing refl term?
                // Check which terms are summed
                surfaceColor += nearObj->surfaceColor * nearObj->reflection * (lambert + phong) + oi.emissionColor;

                // std::cout << "Surface color: " << surfaceColor << std::endl;
            }
        }

    // }

    return surfaceColor + nearObj->emissionColor;
}

void render(Image & img, const vector<Sphere<double> > & obj, const vector<Sphere<double> > & lights)
{
    // create ray in camera space (origin + direction)
    // transform into world space
    // construct ray + near and far clipping plane
    // trace

    int nPixels = img.getNPixels();

    double invWidth = 1. / img.width;
    double invHeight = 1. / img.height;

    double fov = 30;

    double aspectRatio = img.width / img.height;

    double angle = tan(M_PI * 0.5 * fov / 180.);

    int samples = 256;
    // int samples = 64;

    random_device rd;
    mt19937 mt(rd());
    normal_distribution<double> nx(0, 1.0);
    normal_distribution<double> ny(0, 1.0);

    std::mutex imgMutex;

    #pragma omp parallel for
    for (int y = 0; y < img.height; ++y) {
        // std::cout << y << std::endl;
        for (int x = 0; x < img.width; ++x) {

            Vector3<int> total(0);

            for (int s = 0; s < samples; s++) {
                double rayX = (2. * ((x + nx(mt) + 0.5) * invWidth) - 1.) * angle * aspectRatio;
                double rayY = (1. - 2. * ((y + ny(mt) + 0.5) * invHeight)) * angle;

                Vec rayDir = Vec(rayX, rayY, 1); // -> camera faces z direction

                Ray ray = Ray(Point(0.), rayDir);
                Vector3<double> pixelColor = trace(ray, obj, lights, 0);

                total += Vector3<int>(min(255, int(pixelColor.x * 255)),
                                      min(255, int(pixelColor.y * 255)),
                                      min(255, int(pixelColor.z * 255)));
            }

            total.x = total.x / samples;
            total.y = total.y / samples;
            total.z = total.z / samples;

            imgMutex.lock();
            img.pixels[y][x] = total;
            imgMutex.unlock();
        }
    }
}

int main(int argc, char const *argv[])
{
    // vector<Sphere<double> > objects = constructScene();

    omp_set_num_threads(4);

    // heightMapBuilder.SetSourceModule(noiseModule);
    // heightMapBuilder.SetDestNoiseMap(heightMap);
    // heightMapBuilder.SetDestSize(3500, 3500);
    // heightMapBuilder.SetBounds(0.0, 35.0, 0.0, 40.0);
    // heightMapBuilder.Build();

    Image img(1024, 1024);
    // Image img(640, 640);

    std::vector<Sphere<double> > objects;
    std::vector<Sphere<double> > lights;

    objects.push_back(Sphere<double>(Vector3<double>(0, 1, 30), 3, Vector3<double>(0, 1, 0), 0.0, 0.75));
    objects.push_back(Sphere<double>(Vector3<double>(5, 3, 30), 4, Vector3<double>(0, 0, 1), 0.0, 0.75));
    objects.push_back(Sphere<double>(Vector3<double>(0, -10003, 30), 10000, Vector3<double>(1), 0.0, 0.8));
    // objects.push_back(Sphere<double>(Vector3<double>(10025, 0, 15), 10000, Vector3<double>(1, 0, 1), 0.0, 0.001));
    objects.push_back(Sphere<double>(Vector3<double>(5, 12, -2), 10, Vector3<double>(1), 0.0, 0.0, Vector3<double>(0.15)));

    // lights.push_back(Sphere<double>(Vector3<double>(500, 500, -500), 100, Vector3<double>(0, 0, 0), 0.0, 1.0, Vector3<double>(100.0, 100.0, 100.0)));
    // lights.push_back(Sphere<double>(Vector3<double>(0, 0, 30), 3, Vector3<double>(0, 1, 0), 0.0, 0.1, Vector3<double>(100.0, 100.0, 100.0)));

    // Vec rayDir = Vec(0, 0, 1); // -> camera faces z direction

    // Ray ray = Ray(Point(0.), rayDir);

    // std::cout << trace(ray, objects, lights, 0) << std::endl;

    // std::cout << lights[0].emissive << std::endl;

    render(img, objects, lights);
    // std::cout << "writing" << std::endl;
    img.write("rt_out.ppm");

    // Sphere<double> s = Sphere<double>(Vector3<double>(0, 0, 150), 10, Vector3<double>(1, 0, 1), 0.25, 0.0);

    // Ray r = Ray(Vector3<double>(0, 0, 0), Vector3<double>(0, 1, 0));

    // cout << s.intersect(r) << endl;
    // std::cout << "exiting" << std::endl;
    // Segfault here when freeing from obj container with > 1 thread?
    return 0;
}
