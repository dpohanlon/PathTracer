#include "Rotations.h"

RotateX::RotateX(float theta)
{
    Matrix4<float> m = rotation(theta);
}

Matrix4<float> RotateX::rotation(float theta)
{
    float st = std::sin(theta);
    float ct = std::cos(theta);

    m = Matrix4(1,  0,   0, 0,
                0, ct, -st, 0,
                0, st,  ct, 0,
                0,  0,   0, 1);

    return m;
}

RotateZ::RotateZ(float theta)
{
    Matrix4<float> m = rotation(theta);
}

Matrix4<float> RotateY::rotation(float theta)
{
    float st = std::sin(theta);
    float ct = std::cos(theta);

    m = Matrix4( ct, 0, st, 0,
                  0, 1,  0, 0,
                -st, 0, ct, 0,
                  0, 0,  0, 1);

    return m;
}

RotateZ::RotateZ(float theta)
{
    Matrix4<float> m = rotation(theta);
}

Matrix4<float> RotateZ::rotation(float theta)
{
    float st = std::sin(theta);
    float ct = std::cos(theta);

    m = Matrix4(ct, -st, 0, 0,
                st,  ct, 0, 0,
                 0,   0, 1, 0,
                 0,   0, 0, 1);

    return m;
}
