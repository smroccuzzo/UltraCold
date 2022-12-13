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

#ifndef ULTRACOLD_VECTOR
#define ULTRACOLD_VECTOR

#include <complex>
#include <assert.h>

#include "mkl.h"

namespace UltraCold
{

    /**
    * @brief A class that represents arrays of numerical elements.
    * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
    *
    * The purpose of this class is to provide a useful way to store arrays of numerical elements in contiguous memory
    * spaces and to use optimized mathematical operations acting on them. \n
    * These include:
    *
    *  - Optimized memory allocation through Intel's MKL memry functions,
    *  - Fourier transforms (both forward and backward),
    *  - Fortran-style access to the array elements
    *
    * \note Explicit instantiations are provided only for types types <code>double</code>
    * and <code>std::complex<double></code>, which is equivalent to the C complex type defined
    * in <code>complex.h</code>.
    *
    */

    template<typename T>
    class Vector
    {

    public:

        // Constructors and movements
        Vector() = default;                      // Default constructor
        Vector(int);                             // One-dimensional Vector
        Vector(int,int);                         // Two-dimensional Vector
        Vector(int,int,int);                     // Three-dimensional Vector
        Vector(const Vector<T>&);                // Copy constructor
        Vector& operator=(const Vector<T>&);     // Copy assignment
        Vector(Vector<T>&&) noexcept;            // Move constructor
        Vector& operator=(Vector<T>&&) noexcept; // Move assignment
        ~Vector();                               // Destructor

        // Re-initializers

        void reinit(int);            // Reinitialize a one-dimensional Vector
        void reinit(int,int);        // Reinitialize a two-dimensional Vector
        void reinit(int,int,int);    // Reinitialize a three-dimensional Vector

        // Basic member functions

        int size();                  // Get the total number of elements
        int order();                 // Get the number of dimensions
        int extent(int);             // Get the extent along a certain direction
        T* data();                   // Get the pointer to the first element of the array
        T& operator() (int);         // Fortran-style access, one-dimensional Vector
        T& operator() (int,int);     // Fortran-style access, two-dimensional Vector
        T& operator() (int,int,int); // Fortran-style access, three-dimensional Vector
        T& operator[] (int);         // C-style access


    private:

        T* elements = nullptr;      // Pointer to the first element of the Vector
        int number_of_dimensions;   // Order of the Vector, i.e. the number of its dimensions
        int number_of_elements;     // Total number of elements
        int extent_0;               // Number of elements along direction 0
        int extent_1;               // Number of elements along direction 1
        int extent_2;               // Number of elements along direction 2

    };

} // namespace UltraCold

#endif // ULTRACOLD_VECTOR