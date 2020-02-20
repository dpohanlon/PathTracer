#ifndef IMAGE_H
#define IMAGE_H

#include "Vector3.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

class Image
{
public:
    Image(void);
    Image(int w, int h);

    // void visit(ImageVisitor visitor);
    void write(std::string filename);
    int getNPixels(void);
    int width;
    int height;

    typedef Vector3<int> pixelRGB;

    std::vector<std::vector<pixelRGB> > pixels; 

private:

    std::string formatPPM(void);
};

#endif
