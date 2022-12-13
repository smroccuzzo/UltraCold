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

#include "cudaGPSolver.cuh"
#include "mesh_fourier_space.hpp"
#include "DataWriter.hpp"
#include "cub/cub.cuh"
#include "cufft.h"
#include "simple_kernels.cuh"
#include "solver_kernels.cuh"

#define PI 3.1415926535897932384626433
#define TWOPI (2*PI)

namespace UltraCold
{
    namespace cudaSolvers
    {

        /**
         * @brief Constructor for 1d problems
         * */

        GPSolver::GPSolver(Vector<double> &x,
                           Vector<std::complex<double>> &psi_0,
                           Vector<double> &Vext,
                           double scattering_length)
        {

            // Check the order and extent of the Vectors provided
            assert(x.order()==1);
            assert(psi_0.order() == 1);
            assert(Vext.order() == 1);
            nx=x.extent(0);
            assert(psi_0.extent(0) == nx);
            assert(Vext.extent(0) == nx);
            problem_is_1d=true;
            npoints=nx;

            // Initialize the thread grid, i.e. choose the number of cuda threads per block and the number of blocks.
            blockSize = 512;
            gridSize = (npoints + blockSize - 1) / blockSize;

            // Allocate memory for all device arrays
            cudaMalloc(&external_potential_d,nx*sizeof(double));
            cudaMalloc(&kmod2_d,             nx*sizeof(double));
            cudaMalloc(&density_d,           nx*sizeof(double));
            cudaMalloc(&wave_function_d,     nx*sizeof(cuDoubleComplex));
            cudaMalloc(&hpsi_d,              nx*sizeof(cuDoubleComplex));
            cudaMalloc(&ft_wave_function_d,  nx*sizeof(cuDoubleComplex));

            // Allocate space for device and managed scalars
            cudaMalloc(&scattering_length_d,sizeof(double));
            cudaMallocManaged(&norm_d,              sizeof(double));
            cudaMallocManaged(&initial_norm_d,      sizeof(double));
            cudaMallocManaged(&chemical_potential_d,sizeof(double));

            // Get the first necessary copies of input data from host to device
            cudaMemcpy(external_potential_d,Vext.data(),       npoints*sizeof(double),         cudaMemcpyHostToDevice);
            cudaMemcpy(wave_function_d,     psi_0.data(),      npoints*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
            cudaMemcpy(scattering_length_d, &scattering_length,1      *sizeof(double),         cudaMemcpyHostToDevice);

            // Initialize the mesh in Fourier space, and copy it to the device
            Vector<double> kx(nx);
            Vector<double> kmod2(nx);
            create_mesh_in_Fourier_space(x,kx);
            for (size_t i = 0; i < nx; ++i) kmod2(i) = std::pow(kx(i),2);
            cudaMemcpy(kmod2_d,  kmod2.data(),npoints*sizeof(double),cudaMemcpyHostToDevice);

            // Initialize space steps
            dx = x(1)-x(0);
            dv = dx;

            // Initialize the device reduce kernel
            cub::DeviceReduce::Sum(temporary_storage_d,size_temporary_storage,density_d,norm_d,npoints);
            cudaDeviceSynchronize();

            // Allocate temporary storage memory, required for reduction kernels
            cudaMalloc(&temporary_storage_d,size_temporary_storage);
            cudaDeviceSynchronize();

            // Calculate initial norm
            calculate_density(density_d,wave_function_d,npoints);
            cudaDeviceSynchronize();
            cub::DeviceReduce::Sum(temporary_storage_d,size_temporary_storage,density_d,norm_d,npoints);
            cudaDeviceSynchronize();
            norm_d[0]=norm_d[0]*dv;
            initial_norm_d[0]=norm_d[0];
            std::cout << "Initial norm: " << initial_norm_d[0] << std::endl;

            // Initialize the wave function to return as a result
            result_wave_function.reinit(nx);

            // Initialize the host and device vectors containing the mesh axis. This can be useful in particular for
            // data output
            x_axis.reinit(nx);
            x_axis=x;
            kx_axis.reinit(nx);
            kx_axis=kx;
            cudaMalloc(&x_axis_d,nx*sizeof(double));
            cudaMemcpy(x_axis_d,x_axis.data(),nx*sizeof(double),cudaMemcpyHostToDevice);
            cudaMalloc(&kx_axis_d,nx*sizeof(double));
            cudaMemcpy(kx_axis_d,kx_axis.data(),nx*sizeof(double),cudaMemcpyHostToDevice);
            r2mod.reinit(nx);
            for(int i = 0; i < nx; ++i)
                r2mod[i] = std::pow(x(i),2);
            cudaMalloc(&r2mod_d,nx*sizeof(double));
            cudaMemcpy(r2mod_d,r2mod.data(),nx*sizeof(double),cudaMemcpyHostToDevice);

        }

        /**
         * @brief Constructor for 2d problems
         */

        GPSolver::GPSolver(Vector<double> &x,
                           Vector<double> &y,
                           Vector<std::complex<double>> &psi_0,
                           Vector<double> &Vext,
                           double scattering_length)
        {

            // Check the order and extent of the Vectors provided
            assert(x.order()==1);
            assert(y.order()==1);
            assert(psi_0.order() == 2);
            assert(Vext.order() == 2);
            nx=x.extent(0);
            ny=y.extent(0);
            assert(psi_0.extent(0) == nx);
            assert(psi_0.extent(1) == ny);
            assert(Vext.extent(0) == nx);
            assert(Vext.extent(1) == ny);
            problem_is_2d=true;
            npoints=nx*ny;

            // Initialize the thread grid, i.e. choose the number of cuda threads per block and the number of blocks.
            blockSize = 512;
            gridSize = (npoints + blockSize - 1) / blockSize;

            // Allocate memory for all device arrays
            cudaMalloc(&external_potential_d,npoints*sizeof(double));
            cudaMalloc(&kmod2_d,             npoints*sizeof(double));
            cudaMalloc(&density_d,           npoints*sizeof(double));
            cudaMalloc(&wave_function_d,     npoints*sizeof(cuDoubleComplex));
            cudaMalloc(&hpsi_d,              npoints*sizeof(cuDoubleComplex));
            cudaMalloc(&ft_wave_function_d,  npoints*sizeof(cuDoubleComplex));

            // Allocate space for device and managed scalars
            cudaMalloc(&scattering_length_d,sizeof(double));
            cudaMallocManaged(&norm_d,              sizeof(double));
            cudaMallocManaged(&initial_norm_d,      sizeof(double));
            cudaMallocManaged(&chemical_potential_d,sizeof(double));

            // Get the first necessary copies of input data from host to device
            cudaMemcpy(external_potential_d,Vext.data(),       npoints*sizeof(double),         cudaMemcpyHostToDevice);
            cudaMemcpy(wave_function_d,     psi_0.data(),      npoints*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
            cudaMemcpy(scattering_length_d, &scattering_length,1      *sizeof(double),         cudaMemcpyHostToDevice);

            // Initialize the mesh in Fourier space, and copy it to the device
            Vector<double> kx(nx);
            Vector<double> ky(ny);
            Vector<double> kmod2(nx,ny);
            create_mesh_in_Fourier_space(x,y,kx,ky);
            for (size_t i = 0; i < nx; ++i)
                for (size_t j = 0; j < ny; ++j)
                    kmod2(i,j) = std::pow(kx(i),2) +
                                 std::pow(ky(j),2);
            cudaMemcpy(kmod2_d,kmod2.data(),npoints*sizeof(double),cudaMemcpyHostToDevice);

            // Initialize space steps
            dx = x(1)-x(0);
            dy = y(1)-y(0);
            dv = dx*dy;

            // Initialize the device reduce kernel
            cub::DeviceReduce::Sum(temporary_storage_d,size_temporary_storage,density_d,norm_d,npoints);
            cudaDeviceSynchronize();

            // Allocate temporary storage memory, required for reduction kernels
            cudaMalloc(&temporary_storage_d,size_temporary_storage);
            cudaDeviceSynchronize();

            // Calculate initial norm
            calculate_density(density_d,wave_function_d,npoints);
            cudaDeviceSynchronize();
            cub::DeviceReduce::Sum(temporary_storage_d,size_temporary_storage,density_d,norm_d,npoints);
            cudaDeviceSynchronize();
            norm_d[0]=norm_d[0]*dv;
            initial_norm_d[0]=norm_d[0];
            std::cout << "Initial norm: " << initial_norm_d[0] << std::endl;

            // Initialize the wave function to return as a result
            result_wave_function.reinit(nx,ny);

            // Initialize the host vectors containing the mesh axis. This can be useful in particular for data output
            x_axis.reinit(nx);
            y_axis.reinit(ny);
            x_axis=x;
            y_axis=y;
            kx_axis.reinit(nx);
            ky_axis.reinit(ny);
            kx_axis=kx;
            ky_axis=ky;
            cudaMalloc(&x_axis_d,nx*sizeof(double));
            cudaMalloc(&y_axis_d,ny*sizeof(double));
            cudaMemcpy(x_axis_d,x_axis.data(),nx*sizeof(double),cudaMemcpyHostToDevice);
            cudaMemcpy(y_axis_d,y_axis.data(),ny*sizeof(double),cudaMemcpyHostToDevice);
            cudaMalloc(&kx_axis_d,nx*sizeof(double));
            cudaMalloc(&ky_axis_d,ny*sizeof(double));
            cudaMemcpy(kx_axis_d,kx_axis.data(),nx*sizeof(double),cudaMemcpyHostToDevice);
            cudaMemcpy(ky_axis_d,ky_axis.data(),ny*sizeof(double),cudaMemcpyHostToDevice);
            r2mod.reinit(nx,ny);
            for(int i = 0; i < nx; ++i)
                for(int j = 0; j < ny; ++j)
                    r2mod(i,j) = std::pow(x(i),2)+std::pow(y(j),2);
            cudaMalloc(&r2mod_d,npoints*sizeof(double));
            cudaMemcpy(r2mod_d,r2mod.data(),npoints*sizeof(double),cudaMemcpyHostToDevice);
        }

        /**
         * @brief Constructor for 3d problems
         */

        GPSolver::GPSolver(Vector<double> &x,
                           Vector<double> &y,
                           Vector<double> &z,
                           Vector<std::complex<double>> &psi_0,
                           Vector<double> &Vext,
                           double scattering_length)
        {

            // Check the order and extent of the Vectors provided
            assert(x.order()==1);
            assert(y.order()==1);
            assert(z.order()==1);
            assert(psi_0.order() == 3);
            assert(Vext.order() == 3);
            nx=x.extent(0);
            ny=y.extent(0);
            nz=z.extent(0);
            assert(psi_0.extent(0) == nx);
            assert(psi_0.extent(1) == ny);
            assert(psi_0.extent(2) == nz);
            assert(Vext.extent(0) == nx);
            assert(Vext.extent(1) == ny);
            assert(Vext.extent(2) == nz);
            problem_is_3d=true;
            npoints=nx*ny*nz;

            // Initialize the thread grid, i.e. choose the number of cuda threads per block and the number of blocks.
            blockSize = 512;
            gridSize = (npoints + blockSize - 1) / blockSize;

            // Allocate memory for all device arrays
            cudaMalloc(&external_potential_d,npoints*sizeof(double));
            cudaMalloc(&kmod2_d,             npoints*sizeof(double));
            cudaMalloc(&density_d,           npoints*sizeof(double));
            cudaMalloc(&wave_function_d,     npoints*sizeof(cuDoubleComplex));
            cudaMalloc(&hpsi_d,              npoints*sizeof(cuDoubleComplex));
            cudaMalloc(&ft_wave_function_d,  npoints*sizeof(cuDoubleComplex));

            // Allocate space for device and managed scalars
            cudaMalloc(&scattering_length_d,sizeof(double));
            cudaMallocManaged(&norm_d,              sizeof(double));
            cudaMallocManaged(&initial_norm_d,      sizeof(double));
            cudaMallocManaged(&chemical_potential_d,sizeof(double));

            // Get the first necessary copies of input data from host to device
            cudaMemcpy(external_potential_d,Vext.data(),       npoints*sizeof(double),         cudaMemcpyHostToDevice);
            cudaMemcpy(wave_function_d,     psi_0.data(),      npoints*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
            cudaMemcpy(scattering_length_d, &scattering_length,1      *sizeof(double),         cudaMemcpyHostToDevice);

            // Initialize the mesh in Fourier space, and copy it to the device
            Vector<double> kx(nx);
            Vector<double> ky(ny);
            Vector<double> kz(nz);
            Vector<double> kmod2(nx,ny,nz);
            create_mesh_in_Fourier_space(x,y,z,kx,ky,kz);
            for (size_t i = 0; i < nx; ++i)
                for (size_t j = 0; j < ny; ++j)
                    for (size_t k = 0; k < nz; ++k)
                        kmod2(i,j,k) = std::pow(kx(i),2)+
                                       std::pow(ky(j),2)+
                                       std::pow(kz(k),2);
            cudaMemcpy(kmod2_d, kmod2.data(),npoints*sizeof(double),cudaMemcpyHostToDevice);

            // Initialize space steps
            dx = x(1)-x(0);
            dy = y(1)-y(0);
            dz = z(1)-z(0);
            dv = dx*dy*dz;

            // Initialize the device reduce kernel
            cub::DeviceReduce::Sum(temporary_storage_d,size_temporary_storage,density_d,norm_d,npoints);
            cudaDeviceSynchronize();

            // Allocate temporary storage memory, required for reduction kernels
            cudaMalloc(&temporary_storage_d,size_temporary_storage);
            cudaDeviceSynchronize();

            // Calculate initial norm
            calculate_density(density_d,wave_function_d,npoints);
            cudaDeviceSynchronize();
            cub::DeviceReduce::Sum(temporary_storage_d,size_temporary_storage,density_d,norm_d,npoints);
            cudaDeviceSynchronize();
            norm_d[0]=norm_d[0]*dv;
            initial_norm_d[0]=norm_d[0];
            std::cout << "Initial norm: " << initial_norm_d[0] << std::endl;

            // Initialize the wave function to return as a result
            result_wave_function.reinit(nx,ny,nz);

            // Initialize the host vectors containing the mesh axis. This can be useful in particular for data output
            x_axis.reinit(nx);
            y_axis.reinit(ny);
            z_axis.reinit(nz);
            x_axis=x;
            y_axis=y;
            z_axis=z;
            kx_axis.reinit(nx);
            ky_axis.reinit(ny);
            kz_axis.reinit(nz);
            kx_axis=kx;
            ky_axis=ky;
            kz_axis=kz;
            cudaMalloc(&x_axis_d,nx*sizeof(double));
            cudaMalloc(&y_axis_d,ny*sizeof(double));
            cudaMalloc(&z_axis_d,nz*sizeof(double));
            cudaMemcpy(x_axis_d,x_axis.data(),nx*sizeof(double),cudaMemcpyHostToDevice);
            cudaMemcpy(y_axis_d,y_axis.data(),ny*sizeof(double),cudaMemcpyHostToDevice);
            cudaMemcpy(z_axis_d,z_axis.data(),nz*sizeof(double),cudaMemcpyHostToDevice);
            cudaMalloc(&kx_axis_d,nx*sizeof(double));
            cudaMalloc(&ky_axis_d,ny*sizeof(double));
            cudaMalloc(&kz_axis_d,nz*sizeof(double));
            cudaMemcpy(kx_axis_d,kx_axis.data(),nx*sizeof(double),cudaMemcpyHostToDevice);
            cudaMemcpy(ky_axis_d,ky_axis.data(),ny*sizeof(double),cudaMemcpyHostToDevice);
            cudaMemcpy(kz_axis_d,kz_axis.data(),nz*sizeof(double),cudaMemcpyHostToDevice);
            r2mod.reinit(nx,ny,nz);
            for(int i = 0; i < nx; ++i)
                for(int j = 0; j < ny; ++j)
                    for(int k = 0; k < nz; ++k)
                        r2mod(i,j,k) = std::pow(x(i),2)+std::pow(y(j),2)+std::pow(z(k),2);
            cudaMalloc(&r2mod_d,npoints*sizeof(double));
            cudaMemcpy(r2mod_d,r2mod.data(),npoints*sizeof(double),cudaMemcpyHostToDevice);

            // Initialize a vector to output the wave function
            wave_function_output.reinit(nx,ny,nz);
        }

        /**
         * @brief Destructor frees device memory
         *
         */

        GPSolver::~GPSolver()
        {
            cudaFree(external_potential_d);
            cudaFree(density_d);
            cudaFree(norm_d);
            cudaFree(initial_norm_d);
            cudaFree(wave_function_d);
            cudaFree(ft_wave_function_d);
            cudaFree(hpsi_d);
            cudaFree(x_axis_d);
            cudaFree(y_axis_d);
            cudaFree(z_axis_d);
            cudaFree(kx_axis_d);
            cudaFree(ky_axis_d);
            cudaFree(kz_axis_d);
            cudaFree(kmod2_d);
            cudaFree(r2mod_d);
            cudaFree(chemical_potential_d);
            cudaFree(scattering_length_d);
            cudaFree(alpha_d);
            cudaFree(beta_d);
            cudaFree(time_step_d);
            cudaFree(temporary_storage_d);

        }

        /**
         * @brief Get a pointer to the wave function stored in the device
         * */

        cuDoubleComplex* GPSolver::get_wave_function_device_pointer()
        {
            return wave_function_d;
        }

        /**
         * @brief Get a pointer to the fourier transform of the wave function stored in the device
         * */

        cuDoubleComplex* GPSolver::get_ft_wave_function_device_pointer()
        {
            return ft_wave_function_d;
        }

        /**
         * @brief Calculate the density profile
         *
         * */

        void GPSolver::calculate_density(double *density, cuDoubleComplex *wave_function,int size)
        {
            SimpleKernels::square_vector<<<gridSize,blockSize>>>(density,wave_function,size);
        }

        /**
         *
         * @brief Run the gradient descent
         *
         * \warning No check of the residual!
         *
         */

        std::tuple<Vector<std::complex<double>>, double>
                GPSolver::run_gradient_descent(int max_num_iter,
                                               double alpha,
                                               double beta,
                                               std::ostream &output_stream,
                                               int write_output_every)

        {
            // Initialize the fft plan required for the calculation of the laplacian
            cufftHandle ft_plan;
            if(problem_is_1d)
                cufftPlan1d(&ft_plan,nx,CUFFT_Z2Z,1);
            else if(problem_is_2d)
                cufftPlan2d(&ft_plan,nx,ny,CUFFT_Z2Z);
            else if(problem_is_3d)
                cufftPlan3d(&ft_plan,nx,ny,nz,CUFFT_Z2Z);

            //--------------------------------------------------//
            //    Here the gradient-descent iterations start    //
            //--------------------------------------------------//

            // Allocate space for some new data on the device
            cudaMalloc(&alpha_d,sizeof(double));
            cudaMemcpy(alpha_d,&alpha,sizeof(double),cudaMemcpyHostToDevice);
            cudaMalloc(&beta_d,sizeof(double));
            cudaMemcpy(beta_d,&beta,sizeof(double),cudaMemcpyHostToDevice);
            cuDoubleComplex* psi_new;
            cuDoubleComplex* psi_old;
            cudaMalloc(&psi_new,npoints*sizeof(cuDoubleComplex));
            cudaMalloc(&psi_old,npoints*sizeof(cuDoubleComplex));

            // Loop starts here
            for (int it = 0; it < max_num_iter; ++it)
            {

                // Calculate the action of the laplacian
                cufftExecZ2Z(ft_plan, wave_function_d, ft_wave_function_d, CUFFT_FORWARD);
                cudaDeviceSynchronize();
                SimpleKernels::vector_multiplication<<<gridSize,blockSize>>>(ft_wave_function_d,kmod2_d,npoints);
                cudaDeviceSynchronize();
                cufftExecZ2Z(ft_plan, ft_wave_function_d, hpsi_d, CUFFT_INVERSE);
                cudaDeviceSynchronize();
                SimpleKernels::rescale<<<gridSize,blockSize>>>(hpsi_d,0.5*pow(TWOPI,2)/npoints,npoints);
                cudaDeviceSynchronize();

                // Calculate the rest of H|psi>
                SolverKernels::step_2_hpsi<<<gridSize,blockSize>>>(hpsi_d,
                                                                    wave_function_d,
                                                                    external_potential_d,
                                                                    scattering_length_d,
                                                                    npoints);
                cudaDeviceSynchronize();

                // Perform a gradient descent (plus heavy-ball) step
                SolverKernels::gradient_descent_step<<<gridSize,blockSize>>>(wave_function_d,
                                                                               hpsi_d,
                                                                               psi_new,
                                                                               psi_old,
                                                                               alpha_d,
                                                                               beta_d,
                                                                               npoints);
                cudaDeviceSynchronize();

                // Normalize the wave function
                SimpleKernels::square_vector<<<gridSize,blockSize>>>(density_d,psi_new,npoints);
                cudaDeviceSynchronize();
                cub::DeviceReduce::Sum(temporary_storage_d,size_temporary_storage,density_d,norm_d,npoints);
                cudaDeviceSynchronize();
                norm_d[0] = norm_d[0]*dv;
                SimpleKernels::rescale<<<gridSize,blockSize>>>(wave_function_d,
                                                                 psi_new,
                                                                 sqrt(initial_norm_d[0]/norm_d[0]),
                                                                 npoints);
                cudaDeviceSynchronize();

                // Calculate the chemical potential
                SimpleKernels::vector_multiplication<<<gridSize,blockSize>>>(density_d,hpsi_d,wave_function_d,npoints);
                cudaDeviceSynchronize();
                cub::DeviceReduce::Sum(temporary_storage_d,size_temporary_storage,density_d,chemical_potential_d,npoints);
                cudaDeviceSynchronize();
                chemical_potential_d[0] = chemical_potential_d[0]*dv/norm_d[0];

                // Eventually write some output
                if(it % write_output_every == 0)
                    write_gradient_descent_output(it);

            }

            // Free the remaining arrays from the device
            cudaFree(psi_new);
            cudaFree(psi_old);

            // Copy out the results
            cudaMemcpy(result_wave_function.data(),
                       wave_function_d,
                       npoints*sizeof(cuDoubleComplex),
                       cudaMemcpyDeviceToHost);
            double result_chemical_potential = chemical_potential_d[0];

            // Return
            return std::make_pair(result_wave_function,result_chemical_potential);

        }

        /**
         *
         * @brief Write gradient descent output
         *
         * */

        void GPSolver::write_gradient_descent_output(int it)
        {
            std::cout << it << " " << chemical_potential_d[0] << std::endl;
        }

        /**
         * @brief Real-time operator splitting
         * */

        void GPSolver::run_operator_splitting(int number_of_time_steps, double time_step, std::ostream &output_stream,
                                              int write_output_every)
        {
            // Copy input data into the device
            cudaMallocManaged(&time_step_d,sizeof(double));
            cudaMemcpy(time_step_d,&time_step,sizeof(double),cudaMemcpyHostToDevice);

            // Initialize the fft plan required for the calculation of the laplacian
            cufftHandle ft_plan;
            if(problem_is_1d)
                cufftPlan1d(&ft_plan,nx,CUFFT_Z2Z,1);
            else if(problem_is_2d)
                cufftPlan2d(&ft_plan,nx,ny,CUFFT_Z2Z);
            else if(problem_is_3d)
                cufftPlan3d(&ft_plan,nx,ny,nz,CUFFT_Z2Z);

            //----------------------------------------------------//
            //    Here the operator-splitting iterations start    //
            //----------------------------------------------------//

            for (size_t it = 0; it < number_of_time_steps; ++it)
            {

                // Write output starting from the very first iteration
                if(it % write_output_every == 0)
                {
                    cudaMemcpy(wave_function_output.data(),
                               wave_function_d,
                               npoints*sizeof(cuDoubleComplex),
                               cudaMemcpyDeviceToHost);
                    write_operator_splitting_output(it,output_stream);
                }

                // Solve step-1 of operator splitting, i.e. the one NOT involving Fourier transforms
                SolverKernels::step_1_operator_splitting<<<gridSize,blockSize>>>(wave_function_d,
                                                                                   external_potential_d,
                                                                                   time_step_d,
                                                                                   scattering_length_d,
                                                                                   npoints);
                cudaDeviceSynchronize();

                // Solve step-2 of operator splitting, i.e. the one actually involving Fourier transforms
                cufftExecZ2Z(ft_plan,wave_function_d,ft_wave_function_d,CUFFT_FORWARD);
                cudaDeviceSynchronize();
                SolverKernels::aux_step_2_operator_splitting<<<gridSize,blockSize>>>(ft_wave_function_d,
                                                                                       kmod2_d,
                                                                                       time_step_d,
                                                                                       npoints);
                cudaDeviceSynchronize();
                cufftExecZ2Z(ft_plan,ft_wave_function_d,wave_function_d,CUFFT_INVERSE);
                cudaDeviceSynchronize();
                SimpleKernels::rescale<<<gridSize,blockSize>>>(wave_function_d,1./npoints,npoints);
                cudaDeviceSynchronize();

            }

        }

        /**
         *
         * @brief Operator splitting output.
         *
         * This function is called after a copy of the current wave function outside of the GPU, is such a way that it
         * can be used for example for data analysis or to write it to a file for visualization. Since each call
         * blocks the real-time evolution on the GPU until the function has finished, it is better to use it with
         * moderation to avoid a big loss of performance.
         *
         */

        void GPSolver::write_operator_splitting_output(int it,std::ostream& output_stream)
        {}

        /**
         *
         * @brief Reinitialize the solver with a new external potential and wave function.
         *
         * */

        void GPSolver::reinit(Vector<std::complex<double>> &psi, Vector<double> &Vext)
        {
            cudaMemcpy(wave_function_d,psi.data(),npoints*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
            cudaMemcpy(external_potential_d,Vext.data(),npoints*sizeof(double),cudaMemcpyHostToDevice);
        }
    }
}
