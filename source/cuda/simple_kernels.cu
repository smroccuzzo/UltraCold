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

#include "simple_kernels.cuh"
#define PI 3.1415926535897932384626433
#define TWOPI (2*PI)

namespace UltraCold
{
    namespace SimpleKernels
    {

        ///////////////////////////////////////////////////////
        // Pure device kernels
        ///////////////////////////////////////////////////////

        /**
         * @brief A useful complex exponential function
         **/

        __device__ cuDoubleComplex complex_exponential(cuDoubleComplex input)
        {
            cuDoubleComplex res;
            double t = expf (input.x);
            sincos (input.y, &res.y, &res.x);
            res.x *= t;
            res.y *= t;
            return res;
        }

        ////////////////////////////////////////////////////
        // Global kernels
        ////////////////////////////////////////////////////

        /**
         *
         * @brief Calculate the square of a complex vector, storing the result in another complex vector
         *
         * */

        __global__ void square_vector(cuDoubleComplex* result,
                                      cuDoubleComplex* input,
                                      int size)
        {
            int index = blockIdx.x * blockDim.x + threadIdx.x;
            int stride = blockDim.x * gridDim.x;
            for (int i = index; i < size; i += stride)
            {
                result[i].x = input[i].x*input[i].x +
                              input[i].y*input[i].y;
                result[i].y = 0.0;
            }
        }
        /**
         *
         * @brief Calculate the square of a complex vector
         *
         * */

        __global__ void square_vector(double* result,
                                      cuDoubleComplex* input,
                                      int size)
        {
            int index = blockIdx.x * blockDim.x + threadIdx.x;
            int stride = blockDim.x * gridDim.x;
            for (int i = index; i < size; i += stride)
            {
                result[i] = input[i].x*input[i].x +
                            input[i].y*input[i].y;
            }
        }

        /**
         *
         * @brief Calculate the square of a real vector
         *
         * */

        __global__ void square_vector(double* result,
                                      double* input,
                                      int size)
        {
            int index = blockIdx.x * blockDim.x + threadIdx.x;
            int stride = blockDim.x * gridDim.x;
            for (int i = index; i < size; i += stride)
            {
                result[i] = input[i]*input[i];
            }
        }

        /**
         *
         * @brief Multiply a real vector times the square of a complex vector
         *
         * */

        __global__ void vector_average(double* result,
                                       double* input1,
                                       cuDoubleComplex* input2,
                                       int size)
        {
            int index = blockIdx.x * blockDim.x + threadIdx.x;
            int stride = blockDim.x * gridDim.x;
            for (int i = index; i < size; i += stride)
            {
                result[i] = input1[i]*(input2[i].x*input2[i].x+input2[i].y*input2[i].y);
            }
        }

        /**
         *
         * @brief Multiply two complex vectors. Overwrite the first one
         *
         * */

        __global__ void vector_multiplication(cuDoubleComplex* result,
                                              cuDoubleComplex* input,
                                              int size)
        {
            cuDoubleComplex temp;
            int index = blockIdx.x * blockDim.x + threadIdx.x;
            int stride = blockDim.x * gridDim.x;
            for (int i = index; i < size; i += stride)
            {
                temp = result[i];
                result[i].x = temp.x*input[i].x - temp.y*input[i].y;
                result[i].y = temp.y*input[i].x + temp.x*input[i].y;
            }
        }

        /**
         *
         * @brief Multiply a complex and a real vector. Overwrite the complex one
         *
         * */

        __global__ void vector_multiplication(cuDoubleComplex* result,
                                              double* input,
                                              int size)
        {
            int index = blockIdx.x * blockDim.x + threadIdx.x;
            int stride = blockDim.x * gridDim.x;
            for (int i = index; i < size; i += stride)
            {
                result[i].x = result[i].x*input[i];
                result[i].y = result[i].y*input[i];
            }
        }

        /**
          *
          * @brief Multiply two real vectors. Overwrite the first one
          *
          * */

        __global__ void vector_multiplication(double* result,
                                              double* input,
                                              int size)
        {
            int index = blockIdx.x * blockDim.x + threadIdx.x;
            int stride = blockDim.x * gridDim.x;
            for (int i = index; i < size; i += stride)
            {
                result[i] = result[i]*input[i];
            }
        }

        /**
         *
         * @brief Multiply two complex vectors in the case in which the result is a real one
         *
         * */

        __global__ void vector_multiplication(double* result,
                                              cuDoubleComplex* input1,
                                              cuDoubleComplex* input2,
                                              int size)
        {
            int index = blockIdx.x * blockDim.x + threadIdx.x;
            int stride = blockDim.x * gridDim.x;
            for (int i = index; i < size; i += stride)
            {
                result[i] = input1[i].x * input2[i].x;
            }
        }

        /**
         *
         * @brief Rescale a vector for a given input scalar
         *
         * */

        __global__ void rescale(cuDoubleComplex* result,
                                double input,
                                int size)
        {

            int index = blockIdx.x * blockDim.x + threadIdx.x;
            int stride = blockDim.x * gridDim.x;
            for (int i = index; i < size; i += stride)
            {
                result[i].x = result[i].x * input;
                result[i].y = result[i].y * input;
            }

        }

        /**
         *
         * @brief Rescale a vector for a given input scalar, storing the result in another vector
         *
         * */

        __global__ void rescale(cuDoubleComplex* result,
                                cuDoubleComplex* input1,
                                double input2,
                                int size)
        {

            int index = blockIdx.x * blockDim.x + threadIdx.x;
            int stride = blockDim.x * gridDim.x;
            for (int i = index; i < size; i += stride)
            {
                result[i].x = input1[i].x * input2;
                result[i].y = input1[i].y * input2;
            }

        }

        /**
         * @brief Extract layers orthogonal to the x-axis from three dimensional arrays
         *
         * \note Differently from solver and simple kernels, it is absolutely indispensable to call this kernel
         * **on a two-dimensional grid of blocks**
         *
         * */

        __global__ void extract_layer_x(double* d_layer,
                                        double* d_array,
                                        int layer_index,
                                        int nx,
                                        int ny,
                                        int nz)
        {
            int j = blockIdx.x * blockDim.x + threadIdx.x;
            int k = blockIdx.y * blockDim.y + threadIdx.y;
            if(j < ny && k < nz)
                d_layer[nz*j+k] = d_array[nz*ny*layer_index+nz*j+k];
        }

        /**
         * @brief Extract layers orthogonal to the y-axis from three dimensional arrays
         *
         * \note Differently from solver and simple kernels, it is absolutely indispensable to call this kernel
         * **on a two-dimensional grid of blocks**
         *
         * */

        __global__ void extract_layer_y(double* d_layer,
                                        double* d_array,
                                        int layer_index,
                                        int nx,
                                        int ny,
                                        int nz)
        {
            int i = blockIdx.x * blockDim.x + threadIdx.x;
            int k = blockIdx.y * blockDim.y + threadIdx.y;
            if(i < nx && k < nz)
                d_layer[nz*i+k] = d_array[nz*ny*i+nz*layer_index+k];
        }

        /**
         * @brief Extract layers orthogonal to the z-axis from three dimensional arrays
         *
         * \note Differently from solver and simple kernels, it is absolutely indispensable to call this kernel
         * **on a two-dimensional grid of blocks**
         *
         * */

        __global__ void extract_layer_z(double* d_layer,
                                        double* d_array,
                                        int layer_index,
                                        int nx,
                                        int ny,
                                        int nz)
        {
            int i = blockIdx.x * blockDim.x + threadIdx.x;
            int j = blockIdx.y * blockDim.y + threadIdx.y;
            if(i < nx && j < ny)
                d_layer[ny*i+j] = d_array[nz*ny*i+nz*j+layer_index];
        }

        /**
         * @brief Low pass filter, kill all momenta above a certain threshold
         * @param
         * */

        __global__ void low_pass_filter(cuDoubleComplex* ft_filtered,
                                        cuDoubleComplex* ft_to_filter,
                                        double* kmod2,
                                        double momentum_cutoff,
                                        int size)
        {
            int index = blockIdx.x * blockDim.x + threadIdx.x;
            if(index<size && sqrt(kmod2[index]) < momentum_cutoff)
                ft_filtered[index] = ft_to_filter[index];
        }

        /**
         * @brief Apply a density threshold on a known wave function
         *
         * */

        __global__ void density_threshold(cuDoubleComplex* psi,
                                          double* density,
                                          double n_threshold,
                                          int size)
        {
            int index = blockIdx.x * blockDim.x + threadIdx.x;
            if(index<size)
            {
                double current_density = psi[index].x*psi[index].x+psi[index].y*psi[index].y;
                density[index] = 0.0;
                if(current_density < n_threshold)
                    density[index] = 1.0;
            }
        }
    }

}
