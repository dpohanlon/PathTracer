#ifndef ROTATIONS_H
#define ROTATIONS_H

#include <cmath>

#include "Matrix4.hpp"

class RotateX : public Matrix4<double>
{
public:
    RotateX(double theta);

private:
    Matrix4<double> rotation(double theta);

};

class RotateY : public Matrix4<double>
{
public:
    RotateY(double theta);

private:
    Matrix4<double> rotation(double theta);

};

class RotateZ : public Matrix4<double>
{
public:
    RotateZ(double theta);

private:
    Matrix4<double> rotation(double theta);

};

#endif