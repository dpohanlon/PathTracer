#include "Image.h"

Image::Image(void) :
    width(640), height(480)
{
    std::vector<pixelRGB> pw = std::vector<pixelRGB>(width);

    for (int i = 0; i < height; ++i) {
        pixels.push_back(pw);
    }
}

Image::Image(int w, int h) :
    width(w), height(h)
{
    std::vector<pixelRGB> pw = std::vector<pixelRGB>(width);

    for (int i = 0; i < height; ++i) {
        pixels.push_back(pw);
    }    
}

void Image::write(std::string filename)
{
    std::ofstream imageFile(filename.c_str());

    if (imageFile.is_open()) {
        imageFile << this->formatPPM();
    }

    imageFile.close();
}

int Image::getNPixels(void)
{
    return width * height;
}

std::string Image::formatPPM(void)
{
    std::stringstream ss;

    ss << "P3" << std::endl
       << width << " " << height << std::endl
       << "255" << std::endl;

    for (auto row : pixels) {
        for (auto p : row) {
            ss << p.x << " "
               << p.y << " "
               << p.z << " ";
        }
        ss << std::endl;
    }

    return ss.str();
}
