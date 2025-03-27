#include <usvg/image.hpp>

#include <fstream>

namespace usvg
{
    Image::Image()
        : mResolution(Eigen::Vector2i(0, 0))
        , mData(nullptr)
    {
    }

    Image::Image(const Eigen::Vector2i& resolution)
        : mResolution(resolution)
    {
        mData = std::make_shared<Array3d>();
        mData->setSize(resolution.x() * resolution.y());
    }

    Image::Image(int width, int height)
        : mResolution(Eigen::Vector2i(width, height))
    {
        mData = std::make_shared<Array3d>();
        mData->setSize(width * height);
    }

    Image::~Image() {}

    void Image::setValue(int px, int py, const Element& value) { setValue((Eigen::Index)py * mResolution.x() + px, value); }

    void Image::setValue(const Eigen::Vector2i& pixel, const Element& value) { setValue(pixel.y() * mResolution.x() + pixel.x(), value); }

    void Image::setValue(Eigen::Index linearIndex, const Element& value) { mData->setValue(linearIndex, value); }

    const typename Array3d::Element& Image::getValue(int px, int py) const { return getValue((Eigen::Index)py * mResolution.x() + px); }

    const typename Array3d::Element& Image::getValue(const Eigen::Vector2i& pixel) const { return getValue(pixel.y() * mResolution.x() + pixel.x()); }

    const typename Array3d::Element& Image::getValue(Eigen::Index linearIndex) const { return mData->getValue(linearIndex); }

    const Eigen::Vector2i& Image::getResolution() const { return mResolution; }

    void Image::setResolution(const Eigen::Vector2i& resolution)
    {
        mResolution = resolution;
        if (mData == nullptr)
            mData = std::make_shared<Array3d>();
        mData->setSize(resolution.x() * resolution.y());
    }

    void Image::setResolution(int resx, int resy)
    {
        setResolution(Eigen::Vector2i(resx, resy));
    }

    Eigen::Index Image::getNumPixels() const { return mResolution.x() * mResolution.y(); }

    void Image::setZero()
    {
        if (mData != nullptr)
            mData->setZero();
    }

    std::shared_ptr<Array3d> Image::getArray() { return mData; }

    std::shared_ptr<const Array3d> Image::getArray() const { return mData; }

    void Image::writeBMP(const std::string& path) const
    {
        unsigned char file[14] = {
            'B', 'M',        // magic
            0, 0, 0, 0,      // size in bytes
            0, 0,            // app data
            0, 0,            // app data
            40 + 14, 0, 0, 0 // start of data offset
        };
        unsigned char info[40] = {
            40, 0, 0, 0,      // info hd size
            0, 0, 0, 0,       // width
            0, 0, 0, 0,       // heigth
            1, 0,             // number color planes
            24, 0,            // bits per pixel
            0, 0, 0, 0,       // compression is none
            0, 0, 0, 0,       // image bits size
            0x13, 0x0B, 0, 0, // horz resoluition in pixel / m
            0x13, 0x0B, 0, 0, // vert resolutions (0x03C3 = 96 dpi, 0x0B13 = 72 dpi)
            0, 0, 0, 0,       // #colors in pallete
            0, 0, 0, 0,       // #important colors
        };

        int w = this->getResolution().x();
        int h = this->getResolution().y();

        int padSize  = (4 - (w * 3) % 4) % 4;
        int sizeData = w * h * 3 + h * padSize;
        int sizeAll  = sizeData + sizeof(file) + sizeof(info);

        file[2] = (unsigned char)(sizeAll);
        file[3] = (unsigned char)(sizeAll >> 8);
        file[4] = (unsigned char)(sizeAll >> 16);
        file[5] = (unsigned char)(sizeAll >> 24);

        info[4] = (unsigned char)(w);
        info[5] = (unsigned char)(w >> 8);
        info[6] = (unsigned char)(w >> 16);
        info[7] = (unsigned char)(w >> 24);

        info[8]  = (unsigned char)(h);
        info[9]  = (unsigned char)(h >> 8);
        info[10] = (unsigned char)(h >> 16);
        info[11] = (unsigned char)(h >> 24);

        info[20] = (unsigned char)(sizeData);
        info[21] = (unsigned char)(sizeData >> 8);
        info[22] = (unsigned char)(sizeData >> 16);
        info[23] = (unsigned char)(sizeData >> 24);

        std::ofstream stream(path, std::ios::binary | std::ios::out);
        stream.write((char*)file, sizeof(file));
        stream.write((char*)info, sizeof(info));

        unsigned char pad[3] = { 0, 0, 0 };

        // devirtualization, since we would otherwise need to make many virtual function calls.
        for (int y = 0; y < h; y++)
        {
            for (int x = 0; x < w; x++)
            {
                Eigen::Vector3d color  = getValue((h - 1 - y) * w + x);
                unsigned char pixel[3] = { 0 };
                pixel[0]               = (unsigned char)(std::min(std::max(0., color.z()), 1.) * 255);
                pixel[1]               = (unsigned char)(std::min(std::max(0., color.y()), 1.) * 255);
                pixel[2]               = (unsigned char)(std::min(std::max(0., color.x()), 1.) * 255);
                stream.write((char*)pixel, 3);
            }
            stream.write((char*)pad, padSize);
        }
    }
}
