#pragma once

#include <stdint.h>

#include "types.hpp"

namespace usvg
{
    /**
     * @brief Type-safe data array.
     * @tparam TElement Type of element in the array.
     */
    template <typename TElement>
    class Array
    {
    public:
        /**
         * @brief Underlying scalar type.
         */
        using Scalar = typename TElement::Scalar;

        /**
         * @brief Number of components in each of the element entries.
         */
        static constexpr int Dimensions = sizeof(TElement) / sizeof(Scalar);

        /**
         * @brief Underlying Eigen Matrix type.
         */
        using Matrix = Eigen::Matrix<Scalar, Dimensions, -1>;

        /**
         * @brief Underlying element vector stored at each array index.
         */
        using Element = TElement;

        /**
         * @brief Gets the underlying vector.
         * @return std::vector that stores the data in the array.
         */
        [[nodiscard]] inline const std::vector<TElement>& getData() const&
        {
            return mData;
        }

        /**
         * @brief Gets the underlying vector.
         * @return std::vector that stores the data in the array.
         */
        [[nodiscard]] inline std::vector<TElement>& getData() &
        {
            return mData;
        }

        /**
         * @brief Gets the underlying vector.
         * @return std::vector that stores the data in the array.
         */
        [[nodiscard]] inline std::vector<TElement>&& getData() &&
        {
            return std::move(mData);
        }

        /**
         * @brief Gets the number of elements.
         * @return Number of elements in the array.
         */
        [[nodiscard]] inline Eigen::Index getSize() const
        {
            return mData.size();
        }

        /**
         * @brief Get the size in bytes of a single element.
         * @return Size in bytes of each element.
         */
        [[nodiscard]] inline Eigen::Index getElementSizeInBytes() const
        {
            return sizeof(Element);
        }

        /**
         * @brief Get the size in bytes of the entire array.
         * @return Total size in bytes.
         */
        [[nodiscard]] inline Eigen::Index getSizeInBytes() const
        {
            return getSize() * getElementSizeInBytes();
        }

        /**
         * @brief Sets the number of elements.
         * @param size New number of elements.
         */
        inline void setSize(Eigen::Index size)
        {
            mData.resize(size);
        }

        /**
         * @brief Gets the element at a specific index in its internal format.
         * @param index Array index.
         * @return Element at a certain array index.
         */
        [[nodiscard]] inline const Element& getValue(Eigen::Index index) const
        {
            return mData[index];
        }

        /**
         * @brief Sets the element at a specific index with all its components in its internal format.
         * @param index Array index.
         * @param value New element to set at the array index.
         */
        inline void setValue(Eigen::Index index, const Element& value)
        {
            mData[index] = value;
        }

        /**
         * @brief Sets all values to zero.
         */
        inline void setZero()
        {
            std::memset(mData.data(), 0, getSizeInBytes());
        }

        /**
         * @brief Deletes all elements.
         */
        inline void clear()
        {
            mData.clear();
        }

        /**
         * @brief Appends a value at the end.
         * @param value Value to append.
         */
        void append(const Element& value)
        {
            setSize(getSize() + 1);
            setValue(getSize() - 1, value);
        }

        /**
         * @brief Appends an entire array at the end.
         * @param other Array to append.
         */
        void append(const Array& other)
        {
            std::size_t oldsize = getSize();
            setSize(oldsize + other.getSize());
            for (int i = 0; i < other.getSize(); ++i)
                setValue(oldsize + i, other.getValue(i));
        }

        /**
         * @brief Reverses the order of the elements.
         */
        void reverse()
        {
            const Eigen::Index n = getSize();
            for (Eigen::Index i = 0; i < n / 2; ++i)
            {
                Element temp = getValue(i);
                setValue(i, getValue(n - 1 - i));
                setValue(n - 1 - i, temp);
            }
        }

        /**
         * @brief Removes the last n elements of this vector.
         * @param n Number of elements to remove at the end.
         */
        inline void removeLast(std::size_t n = 1)
        {
            setSize(getSize() - n);
        }

        /**
         * @brief Removes the first n elements of this vector.
         * @param n Number of elements to remove at the front.
         */
        void removeFirst(std::size_t n = 1)
        {
            // Copy each element n places ahead
            for (std::size_t i = n; i < getSize(); ++i)
                setValue(i - n, getValue(i));
            removeLast(n);
        }

        /**
         * @brief Get the first element.
         * @return First element.
         */
        inline const Element& first() const
        {
            return getValue(0);
        }

        /**
         * @brief Get the last element.
         * @return Last element.
         */
        inline const Element& last() const
        {
            return getValue(getSize() - 1);
        }

    protected:
        /**
         * @brief Vector containing the elements of the array.
         */
        std::vector<Element> mData;
    };

    /**
     * @brief Dense one-dimensional array of double values with one component.
     */
    using Array1d = Array<Eigen::Vector1d>;

    /**
     * @brief Dense one-dimensional array of double values with two components.
     */
    using Array2d = Array<Eigen::Vector2d>;

    /**
     * @brief Dense one-dimensional array of double values with three components.
     */
    using Array3d = Array<Eigen::Vector3d>;

    /**
     * @brief Dense one-dimensional array of double values with four components.
     */
    using Array4d = Array<Eigen::Vector4d>;
}
