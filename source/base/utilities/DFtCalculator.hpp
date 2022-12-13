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

#ifndef ULTRACOLD_MKL_FOURIER_TRANSFORM
#define ULTRACOLD_MKL_FOURIER_TRANSFORM

#include <complex>
#include <assert.h>

#include "Vector.hpp"
#include "mkl.h"

namespace UltraCold
{

    namespace MKLWrappers
    {

        /**
         * @brief Class to calculate Fourier transforms using Intel's MKL DFT functions
         * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
         *
         * The purpose of this class is to provide a simple way for calculating Fourier transforms using Intel's MKL
         * Discrete Fourier Transform functions. \n
         * The class is used as follows
         * \code{.cpp}
         * 
         * DFtCalculator dft(forward_domain_vector,backward_domain_vector);
         * dft.compute_forward(); // For a forward transform
         * 
         * // or 
         * 
         * dft.compute_backward(); // for a backward transform
         * 
         * \endcode 
         * 
         * Here, <code>forward_domain_vector</code>  is the real-space vector, while
         * <code>backward_domain_vector</code> is the Fourier-space vector.
         * \note The two Vectors must be different (out-of-place transform). Moreover, for performance and portability
         * reasons, a <code>DFtCalculator</code> object must be associated with **each** forward-backward pair of
         * vectors you want to use.
         * 
         */

        class DFtCalculator
        {
        
        public:

            // Constructors
            DFtCalculator(Vector<std::complex<double>>& forward_domain_vector,
                          Vector<std::complex<double>>& backward_domain_vector);
            DFtCalculator(Vector<double>& forward_domain_vector,
                          Vector<std::complex<double>>& backward_domain_vector);
            DFtCalculator() = default;

            // Destructor
            ~DFtCalculator();

            // Reinit functions
            void reinit(Vector<std::complex<double>>& forward_domain_vector,
                        Vector<std::complex<double>>& backward_domain_vector);
            void reinit(Vector<double>& forward_domain_vector,
                        Vector<std::complex<double>>& backward_domain_vector);

            // Calculator functions
            void compute_forward();
            void compute_backward();

        private:

            // The descriptor
            DFTI_DESCRIPTOR_HANDLE dft_descriptor;

            // Extents
            long forward_domain_n1,forward_domain_n2,forward_domain_n3;
            long backward_domain_n1,backward_domain_n2,backward_domain_n3;

            // status variable
            long status;

            // pointer to first memory location of input/output Vectors
            void* forward_domain_pointer  = nullptr;
            void* backward_domain_pointer = nullptr;

            // Transform type
            bool transform_is_real_complex = false;
            bool transform_is_2d = false;
            bool transform_is_3d = false;

            // Strides for 2d and 3d transforms. Needed only for real-complex transforms
            long strides_forward_domain_2d[3];
            long strides_forward_domain_3d[4];
            long strides_backward_domain_2d[3];
            long strides_backward_domain_3d[4];

    };
        
    } // namespace MKLWrappers
} // namespace UltraCold


#endif // ULTRACOLD_MKL_FOURIER_TRANSFORM