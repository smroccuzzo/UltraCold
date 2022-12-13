/*--------------------------------------------------------------------------------
*
*    This file is part of the UltraCold project.
*
*    UltraCold is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    any later version.
*    UltraCold is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*    You should have received a copy of the GNU General Public License
*    along with UltraCold.  If not, see <https://www.gnu.org/licenses/>.
*
*--------------------------------------------------------------------------------*/

#include "Vector.hpp"

namespace UltraCold{

    /**
     *
     * @brief Constructor for a one-dimensional array
     * @param n1 *int* The extent of the array
     *
     */

    template<typename T>
    Vector<T>::Vector(int n1)
    {

        number_of_dimensions = 1;
        number_of_elements   = n1;
        extent_0 = n1;
        extent_1 = 0;
        extent_2 = 0;

        elements = (T*) mkl_calloc(number_of_elements,sizeof(T),64);

    }

    /**
     *
     * @brief Constructor for a two-dimensional array
     * @param n1  *int* The extent of the array along the first dimension
     * @param n2  *int* The extent of the array along the second dimension
     *
     */

    template<typename T>
    Vector<T>::Vector(int n1, int n2)
    {

        number_of_dimensions = 2;
        number_of_elements   = n1*n2;
        extent_0 = n1;
        extent_1 = n2;
        extent_2 = 0;

        elements = (T*) mkl_calloc(number_of_elements,sizeof(T),64);

    }

    /**
     *
     * @brief Constructor for a two-dimensional array
     * @param n1 *int* The extent of the array along the first dimension
     * @param n2 *int* The extent of the array along the second dimension
     * @param n3 *int* The extent of the array along the third dimension
     *
     */

    template<typename T>
    Vector<T>::Vector(int n1, int n2, int n3)
    {
        number_of_dimensions = 3;
        number_of_elements   = n1*n2*n3;
        extent_0 = n1;
        extent_1 = n2;
        extent_2 = n3;

        elements = (T*) mkl_calloc(number_of_elements,sizeof(T),64);

    }

    /**
     *
     * @brief Copy constructor
     *
    */

    template <typename T>
    Vector<T>::Vector(const Vector<T>& w)
    {

        number_of_dimensions = w.number_of_dimensions;
        number_of_elements   = w.number_of_elements;
        extent_0 = w.extent_0;
        extent_1 = w.extent_1;
        extent_2 = w.extent_2;

        elements = (T*) mkl_calloc(number_of_elements,sizeof(T),64);

        for (size_t i = 0; i < number_of_elements; ++i)
            elements[i] = w.elements[i];

    }

    /**
     *
     * @brief Copy assignment
     *
    */

    template <typename T>
    Vector<T>& Vector<T>::operator=(const Vector<T> & w)
    {

        number_of_dimensions = w.number_of_dimensions;
        number_of_elements   = w.number_of_elements;
        extent_0 = w.extent_0;
        extent_1 = w.extent_1;
        extent_2 = w.extent_2;

        elements = (T*) mkl_calloc(number_of_elements,sizeof(T),64);

        for (size_t i = 0; i < number_of_elements; ++i)
            elements[i] = w.elements[i];

        return *this;

    }

    /**
     * @brief Move constructor
     *
    */

    template <typename T>
    Vector<T>::Vector(Vector<T>&& w) noexcept
    {
        elements = w.elements;
        number_of_dimensions = w.number_of_dimensions;
        number_of_elements   = w.number_of_elements;
        extent_0 = w.extent_0;
        extent_1 = w.extent_1;
        extent_2 = w.extent_2;

        w.elements = nullptr;
        w.number_of_elements = 0;
        w.number_of_dimensions = 0;
        w.extent_0 = 0;
        w.extent_1 = 0;
        w.extent_2 = 0;

    }

    /**
     * @brief Move assignment
     *
    */

    template <typename T>
    Vector<T>& Vector<T>::operator=(Vector<T>&& w) noexcept
    {
        elements = w.elements;
        number_of_dimensions = w.number_of_dimensions;
        number_of_elements   = w.number_of_elements;
        extent_0 = w.extent_0;
        extent_1 = w.extent_1;
        extent_2 = w.extent_2;

        w.elements = nullptr;
        w.number_of_elements = 0;
        w.number_of_dimensions = 0;
        w.extent_0 = 0;
        w.extent_1 = 0;
        w.extent_2 = 0;

        return *this;
    }

    /**
     *
     * @brief Destructor, releases the memory allocated for the vector
     *
    */

    template<typename T>
    Vector<T>::~Vector()
    {
        mkl_free(elements);
    }

    /**
     * @brief Reinitialize a one-dimensional Vector
     *
     */

    template <typename T>
    void Vector<T>::reinit(int n1)
    {
        number_of_dimensions = 1;
        number_of_elements   = n1;
        extent_0 = n1;
        extent_1 = 0;
        extent_2 = 0;

        mkl_free(elements);
        elements = (T*) mkl_calloc(number_of_elements,sizeof(T),64);

    }

    /**
     * @brief Reinitialize a two-dimensional Vector
     *
     */

    template <typename T>
    void Vector<T>::reinit(int n1,int n2)
    {
        number_of_dimensions = 2;
        number_of_elements   = n1*n2;
        extent_0 = n1;
        extent_1 = n2;
        extent_2 = 0;

        mkl_free(elements);
        elements = (T*) mkl_calloc(number_of_elements,sizeof(T),64);

    }

    /**
     * @brief Reinitialize a three-dimensional Vector
     *
     */

    template <typename T>
    void Vector<T>::reinit(int n1, int n2, int n3)
    {
        number_of_dimensions = 3;
        number_of_elements   = n1*n2*n3;
        extent_0 = n1;
        extent_1 = n2;
        extent_2 = n3;

        mkl_free(elements);
        elements = (T*) mkl_calloc(number_of_elements,sizeof(T),64);

    }

    /**
     *
     * @brief Get the total number of elements
     *
    */

    template<typename T>
    int Vector<T>::size()
    {
        return number_of_elements;
    }

    /**
     *
     * @brief Get the total number of dimensions
     *
    */

    template<typename T>
    int Vector<T>::order()
    {
        return number_of_dimensions;
    }

    /**
     *
     * @brief Get the extent of the Vector along a certain direction. If a wrong index is given, returns -1.
     * @param i *int* direction along which we want to get the extent
     *
    */

    template<typename T>
    int Vector<T>::extent(int i)
    {
        if (i==0)
        {
            return extent_0;
        }
        else if (i==1)
        {
            return extent_1;
        }
        else if (i==2)
        {
            return extent_2;
        }
        else
        {
            return -1;
        }
    }

    /**
     * @brief Get the pointer to the first array element. Use it with care!
     *
     */

    template<typename T>
    T* Vector<T>::data()
    {
        return elements;
    }

    /**
     *
     * @brief Member access operators in Fortran style for a one-dimensional Vector
     * @param i *int* the index of the requested element
     *
    */

    template<typename T>
    T& Vector<T>::operator()(int i)
    {
        assert(i<extent_0);
        assert(i>=0);
        return elements[i];
    }

    /**
     *
     * @brief Member access operators in Fortran style for a two-dimensional Vector
     * @param i *int* first index of the requested element
     * @param j *int* second index of the requested element
     *
     * \note Only the *syntax* is Fortran-style, while the elements are accessed in C-style row-major order.
     *
    */

    template<typename T>
    T& Vector<T>::operator()(int i, int j)
    {
        assert(i>=0);
        assert(j>=0);
        assert(j<extent_1);
        assert(i<extent_0);
        return elements[extent_1*i+j];
    }


    /**
     *
     * @brief Member access operators in Fortran style for a three-dimensional Vector
     * @param i *int* first index of the requested element
     * @param j *int* second index of the requested element
     * @param k *int* third index of the requested element
     *
     * \note Only the *syntax* is Fortran-style, while the elements are accessed in C-style row-major order.
     *
    */

    template<typename T>
    T& Vector<T>::operator()(int i, int j, int k)
    {
        assert(i>=0);
        assert(j>=0);
        assert(k>=0);
        assert(i<extent_0);
        assert(j<extent_1);
        assert(k<extent_2);
        return elements[extent_1*extent_2*i+extent_2*j+k];
    }

    /**
     *
     * @brief  Plain and simple member access operators in C-style
     * @param i *int* the index of the requested element
     *
    */

    template<typename T>
    T& Vector<T>::operator[](int i)
    {
        assert(i<number_of_elements);
        return elements[i];
    }

    // Explicit template instantiations

    template class Vector<double>;
    template class Vector<std::complex<double>>;

} // namespace UltraCold