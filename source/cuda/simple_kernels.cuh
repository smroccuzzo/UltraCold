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

#ifndef ULTRACOLD_CUDA_SIMPLE_KERNELS
#define ULTRACOLD_CUDA_SIMPLE_KERNELS

#include "cuComplex.h"
#include "cufft.h"

namespace UltraCold
{
    namespace SimpleKernels
    {

        /**
         * @brief This namespace contains the simplest CUDA kernels used in GPU-accelerated UltraCold classes.
         *
         * These kernels are just useful CUDA implementations of simple functions used in the solver classes of
         * UltraCold exploiting GPU acceleration. They are not meant to be used outside of such classes
         *
         * */

        ///////////////////
        // Device kernels
        ///////////////////

        __device__ cuDoubleComplex complex_exponential(cuDoubleComplex input);

        ///////////////////
        // Global kernels
        ///////////////////

        __global__ void square_vector( cuDoubleComplex* result,
                                       cuDoubleComplex* input,
                                       int size);
        __global__ void square_vector( double* result,
                                       cuDoubleComplex* input,
                                       int size);
        __global__ void square_vector(double* result,
                                      double* input,
                                      int size);
        __global__ void vector_average(double* result,
                                       double* input1,
                                       cuDoubleComplex* input2,
                                       int size);
        __global__ void vector_multiplication(cuDoubleComplex* result,
                                              cuDoubleComplex* input,
                                              int size);
        __global__ void vector_multiplication(cuDoubleComplex* result,
                                              double* input,
                                              int size);
        __global__ void vector_multiplication(double* result,
                                              double* input,
                                              int size);
        __global__ void vector_multiplication(double* result,
                                              cuDoubleComplex* input1,
                                              cuDoubleComplex* input2,
                                              int size);
        __global__ void vector_multiplication(double* result,
                                              cuDoubleComplex* v1,
                                              double* v2,
                                              int size);
        __global__ void rescale(cuDoubleComplex* result,
                                double input,
                                int size);
        __global__ void rescale(cuDoubleComplex* result,
                                cuDoubleComplex* input1,
                                double input2,
                                int size);
        __global__ void extract_layer_x(double* d_layer,
                                        double* d_array,
                                        int layer_index,
                                        int nx,
                                        int ny,
                                        int nz);
        __global__ void extract_layer_y(double* d_layer,
                                        double* d_array,
                                        int layer_index,
                                        int nx,
                                        int ny,
                                        int nz);
        __global__ void extract_layer_z(double* d_layer,
                                        double* d_array,
                                        int layer_index,
                                        int nx,
                                        int ny,
                                        int nz);
        __global__ void low_pass_filter(cuDoubleComplex* ft_filtered,
                                        cuDoubleComplex* ft_to_filter,
                                        double* kmod2,
                                        double momentum_cutoff,
                                        int size);
        __global__ void density_threshold(cuDoubleComplex* psi,
                                          double* density,
                                          double n_threshold,
                                          int size);
    }
}

#endif // ULTRACOLD_CUDA_SIMPLE_KERNELS