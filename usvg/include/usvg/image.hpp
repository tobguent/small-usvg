#pragma once

#include "array.hpp"

namespace usvg
{
    /**
     * @brief Two-dimensional image storing three double values per pixel.
     */
    class Image
    {
    public:
        /**
         * @brief Data type for a single pixel.
         */
        using Element = Array3d::Element;

        /**
         * @brief Empty default constructor.
         */
        Image();

        /**
         * @brief Constructor of image with a certain size.
         * @param resolution Image resolution.
         */
        Image(const Eigen::Vector2i& resolution);

        /**
         * @brief Constructor of image with a certain size.
         * @param width Width of image.
         * @param height Height of image.
         */
        Image(int width, int height);

        /**
         * @brief Destructor.
         */
        ~Image();

        /**
         * @brief Sets the value at a certain pixel.
         * @param px x-coordinate of pixel.
         * @param py y-coordinate of pixel.
         * @param value Value to set.
         */
        void setValue(int px, int py, const Element& value);

        /**
         * @brief Sets the value at a certain pixel.
         * @param pixel Pixel index.
         * @param value Value to set.
         */
        void setValue(const Eigen::Vector2i& pixel, const Element& value);

        /**
         * @brief Sets the value at a certain pixel.
         * @param linearIndex Linear array index.
         * @param value Value to set.
         */
        void setValue(Eigen::Index linearIndex, const Element& value);

        /**
         * @brief Gets the value at a certain pixel.
         * @param px x-coordinate of pixel.
         * @param py y-coordinate of pixel.
         * @return Value at pixel.
         */
        [[nodiscard]] const typename Array3d::Element& getValue(int px, int py) const;

        /**
         * @brief Gets the value at a certain pixel.
         * @param pixel Pixel index.
         * @return Value at pixel.
         */
        [[nodiscard]] const typename Array3d::Element& getValue(const Eigen::Vector2i& pixel) const;

        /**
         * @brief Gets the value at a certain pixel.
         * @param linearIndex Linear array index.
         * @return Value at pixel.
         */
        [[nodiscard]] const typename Array3d::Element& getValue(Eigen::Index linearIndex) const;

        /**
         * @brief Gets the resolution of the image.
         * @return Image resolution.
         */
        [[nodiscard]] const Eigen::Vector2i& getResolution() const;

        /**
         * @brief Sets the resolution of the image. The internal linear image data array gets resized.
         * @param resolution New image resolution.
         */
        void setResolution(const Eigen::Vector2i& resolution);

        /**
         * @brief Sets the resolution of the image.
         * @param resx New width of the image.
         * @param resy New height of the image.
         */
        void setResolution(int resx, int resy);

        /**
         * @brief Gets the total number of pixels.
         * @return Number of pixels.
         */
        [[nodiscard]] inline Eigen::Index getNumPixels() const;

        /**
         * @brief Sets all values to zero.
         */
        void setZero();

        /**
         * @brief Gets the array that stores the data in the image.
         * @return Internal data array.
         */
        [[nodiscard]] std::shared_ptr<Array3d> getArray();

        /**
         * @brief Gets the array that stores the data in the image.
         * @return Internal data array.
         */
        [[nodiscard]] std::shared_ptr<const Array3d> getArray() const;

        /**
         * @brief Writes the image to a BMP file.
         * @param path Path to where the BMP file should be written.
         */
        void writeBMP(const std::string& path) const;

    private:
        /**
         * @brief Image resolution.
         */
        Eigen::Vector2i mResolution;

        /**
         * @brief Data stored in the image.
         */
        std::shared_ptr<Array3d> mData;
    };
}
