#ifndef VEC_H
#define VEC_H

#include "Vector3.hpp"

class Vec : public Vector3<float>
{
public:
    // Init base class
    Vec(void) : Vector3<float>() {}
    Vec(float c) : Vector3(c) {}
    Vec(float x, float y, float z) : Vector3(x, y, z) {}
    Vec(Vector3<float> v) : Vector3(v) {}

    float weight = 0;
};

#endif
