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

#include <random>
#include "cudaDipolarGPSolver.cuh"
#include "mesh_fourier_space.hpp"
#include "DataWriter.hpp"
#include "cub/cub.cuh"
#include "simple_kernels.cuh"
#include "solver_kernels.cuh"
#include "DFtCalculator.hpp"

#define PI 3.1415926535897932384626433
#define TWOPI (2*PI)

namespace UltraCold
{
    namespace cudaSolvers
    {

        /**
         * @brief Constructor for 2d problems
         */

        DipolarGPSolver::DipolarGPSolver(Vector<double> &x,
                                         Vector<double> &y,
                                         Vector<std::complex<double>> &psi_0,
                                         Vector<double> &Vext,
                                         double scattering_length,
                                         double dipolar_length,
                                         double alpha)
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

            // Initialize the Fourier transform of the dipolar potential

            double epsilon_dd = 0.0;
            if(scattering_length != 0)
                epsilon_dd = dipolar_length/scattering_length;

            Vtilde.reinit(nx,ny);
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                {
                    double qd = TWOPI * kx[i] / sqrt(2);
                    double q  = TWOPI * sqrt(pow(kx[i], 2) + pow(ky[j], 2)) / sqrt(2);
                    double value =
                            sqrt(8*PI) * scattering_length * epsilon_dd *
                            (
                                    (-1 + 3*sqrt(PI) * pow(qd,2)/q * exp(pow(q,2)) * erfc(q)) * pow(sin(alpha),2) +
                                    ( 2 - 3*sqrt(PI) * q * exp(pow(q,2)) * erfc(q)) * pow(cos(alpha),2)
                            );
                    if (isnan(value) && kx(i) == 0 && ky(j) == 0)
                        Vtilde(i, j) = sqrt(8*PI)*scattering_length*epsilon_dd*(3*std::pow(std::cos(alpha),2)-1);
                    else if(isnan(value) && kx(i) != 0 && ky(j) != 0)
                        Vtilde(i,j) = 0.0;
                    else
                        Vtilde(i, j) = value;
                }

            cudaMalloc(&Vtilde_d,npoints*sizeof(cuDoubleComplex));
            cudaMemcpy(Vtilde_d,Vtilde.data(),npoints*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
            cudaMalloc(&Phi_dd_d,npoints*sizeof(cuDoubleComplex));

            // Initialize gamma(\epsilon_dd) for the LHY correction
            cudaMallocManaged(&gamma_epsilondd_d,sizeof(double));
            gamma_epsilondd_d[0] = 0.0;

            // Scattering length is divided by sqrt(2PI) here, since in the propagators it is multiplied by 4PI
            scattering_length *= 1./sqrt(2*PI);
            cudaMemcpy(scattering_length_d, &scattering_length,1*sizeof(double),cudaMemcpyHostToDevice);

            // Initialize the wave function for the output
            wave_function_output.reinit(nx,ny);

        }

        /**
         * @brief Constructor for 3d problems
         */

        DipolarGPSolver::DipolarGPSolver(Vector<double> &x,
                                         Vector<double> &y,
                                         Vector<double> &z,
                                         Vector<std::complex<double>> &psi_0,
                                         Vector<double> &Vext,
                                         double scattering_length,
                                         double dipolar_length,
                                         double theta_mu,
                                         double phi_mu,
                                         bool add_lhy_correction)
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
                        r2mod(i,j,k) = std::pow(x(i),2);
            cudaMalloc(&r2mod_d,npoints*sizeof(double));
            cudaMemcpy(r2mod_d,r2mod.data(),npoints*sizeof(double),cudaMemcpyHostToDevice);

            // Initialize the Fourier transform of the dipolar potential
            cudaMallocManaged(&epsilon_dd_d,sizeof(double));
            epsilon_dd_d[0] = 0.0;
            if(scattering_length != 0)
                epsilon_dd_d[0] = dipolar_length/scattering_length;
            Vtilde.reinit(nx,ny,nz);
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k)
                    {
                        double aux = TWOPI * (
                                kx[i]*sin(theta_mu)*cos(phi_mu) +
                                ky[j]*sin(theta_mu)*sin(phi_mu)+
                                kz[k]*cos(theta_mu));
                        double aux1 = TWOPI * sqrt(pow(kx[i], 2) + pow(ky[j], 2) + pow(kz[k], 2));
                        if (aux1 <= 1.E-6)
                            Vtilde(i,j,k) = -4*PI*scattering_length*epsilon_dd_d[0];
                        else
                            Vtilde(i,j,k) =
                                    12.0 * PI * scattering_length * epsilon_dd_d[0] * (pow(aux/aux1,2)-1.0/3.0);
                    }
            cudaMalloc(&Vtilde_d,npoints*sizeof(cuDoubleComplex));
            cudaMemcpy(Vtilde_d,Vtilde.data(),npoints*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
            cudaMalloc(&Phi_dd_d,npoints*sizeof(cuDoubleComplex));

            // Initialize gamma(\epsilon_dd) for the LHY correction
            cudaMallocManaged(&gamma_epsilondd_d,sizeof(double));
            gamma_epsilondd_d[0] = 0.0;
            if (epsilon_dd_d[0] != 0 && add_lhy_correction)
            {
                gamma_epsilondd_d[0] = 64.0*sqrt(PI)/3*sqrt(pow(scattering_length,5));
                double F_epsilon_dd=0.0;
                int n_theta=1000;
                double d_theta=PI/(n_theta-1);

                std::complex<double> csum;
                std::complex<double> caux;
                csum={0.0,0.0};
                for (int i = 0; i < n_theta; ++i)
                {
                    double theta=i*d_theta;
                    caux = pow(1.0+epsilon_dd_d[0]*(3.0*pow(cos(theta),2)-1.0),5);
                    caux = sqrt(caux);
                    csum += sin(theta)*caux;
                }
                csum *= d_theta;
                F_epsilon_dd = csum.real();
                gamma_epsilondd_d[0] *= F_epsilon_dd;
            }

            // Initialize the wave function for the output
            wave_function_output.reinit(nx,ny,nz);

        }

        /**
         * @brief Destructor frees device memory
         *
         * */

        DipolarGPSolver::~DipolarGPSolver()
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
            cudaFree(Vtilde_d);
            cudaFree(Phi_dd_d);
            cudaFree(epsilon_dd_d);
            cudaFree(gamma_epsilondd_d);
        }

        /**
         * @brief Get a pointer to the wave function stored in the device
         * */

        cuDoubleComplex* DipolarGPSolver::get_wave_function_device_pointer()
        {
            return wave_function_d;
        }

        /**
         * @brief Get a pointer to the wave function stored on the device
         * */

        cuDoubleComplex* DipolarGPSolver::get_ft_wave_function_device_pointer()
        {
            return ft_wave_function_d;
        }

        /**
         * @brief Calculate the density profile
         *
         * */

        void DipolarGPSolver::calculate_density(double *density, cuDoubleComplex *wave_function,int size)
        {
            SimpleKernels::square_vector<<<gridSize,blockSize>>>(density,wave_function,size);
        }

        /**
         *
         * @brief Run the gradient descent
         *
         * \warning No check of the residual!
         *
         * */

        std::tuple<Vector<std::complex<double>>, double>
        DipolarGPSolver::run_gradient_descent(int max_num_iter,
                                              double alpha,
                                              double beta,
                                              std::ostream &output_stream,
                                              int write_output_every)

        {
            // Initialize the fft plan required for the calculation of the laplacian
            cufftHandle ft_plan;
            if(problem_is_2d)
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
            cuDoubleComplex* c_density_d;
            cudaMalloc(&c_density_d,npoints*sizeof(cuDoubleComplex));

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

                // Calculate the dipolar potential
                SimpleKernels::square_vector<<<gridSize,blockSize>>>(c_density_d,wave_function_d,npoints);
                cudaDeviceSynchronize();
                cufftExecZ2Z(ft_plan,c_density_d,ft_wave_function_d,CUFFT_FORWARD);
                cudaDeviceSynchronize();
                SimpleKernels::vector_multiplication<<<gridSize,blockSize>>>(ft_wave_function_d,
                                                                               Vtilde_d,
                                                                               npoints);
                cudaDeviceSynchronize();
                cufftExecZ2Z(ft_plan,ft_wave_function_d,Phi_dd_d,CUFFT_INVERSE);
                cudaDeviceSynchronize();
                SimpleKernels::rescale<<<gridSize,blockSize>>>(Phi_dd_d,1./npoints,npoints);

                // Calculate the rest of H|psi>
                SolverKernels::step_2_dipolar_hpsi<<<gridSize,blockSize>>>(hpsi_d,
                                                                             wave_function_d,
                                                                             external_potential_d,
                                                                             Phi_dd_d,
                                                                             scattering_length_d,
                                                                             gamma_epsilondd_d,
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
                    write_gradient_descent_output(it,output_stream);

            }

            // Free the remaining arrays from the device
            cudaFree(psi_new);
            cudaFree(psi_old);
            cudaFree(c_density_d);

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

        void DipolarGPSolver::write_gradient_descent_output(size_t        iteration_number,
                                                            std::ostream& output_stream)
        {
            output_stream << iteration_number << " " << chemical_potential_d[0] << std::endl;
        }

        /**
         * @brief Real-time operator splitting
         */

        void DipolarGPSolver::run_operator_splitting(int number_of_time_steps,
                                                     double time_step,
                                                     std::ostream &output_stream,
                                                     int write_output_every)
        {
            // Copy input data into the device
            cudaMallocManaged(&time_step_d,sizeof(double));
            cudaMemcpy(time_step_d,&time_step,sizeof(double),cudaMemcpyHostToDevice);

            // Initialize the fft plan required for the calculation of the laplacian
            cufftHandle ft_plan;
            if(problem_is_2d)
                cufftPlan2d(&ft_plan,nx,ny,CUFFT_Z2Z);
            else if(problem_is_3d)
                cufftPlan3d(&ft_plan,nx,ny,nz,CUFFT_Z2Z);
            cuDoubleComplex* c_density_d;
            cudaMalloc(&c_density_d,npoints*sizeof(cuDoubleComplex));

            // Initialize other variables
            this->write_output_every=write_output_every;

            //----------------------------------------------------//
            //    Here the operator-splitting iterations start    //
            //----------------------------------------------------//
            for (size_t it = 0; it < number_of_time_steps; ++it)
            {

                // Write output starting from the very first iteration
                if(it % write_output_every == 0)
                    write_operator_splitting_output(it,output_stream);

                // Calculate the current value of dipolar potential
                SimpleKernels::square_vector<<<gridSize,blockSize>>>(c_density_d,wave_function_d,npoints);
                cudaDeviceSynchronize();
                cufftExecZ2Z(ft_plan,c_density_d,ft_wave_function_d,CUFFT_FORWARD);
                cudaDeviceSynchronize();
                SimpleKernels::vector_multiplication<<<gridSize,blockSize>>>(ft_wave_function_d,
                                                                               Vtilde_d,
                                                                               npoints);
                cudaDeviceSynchronize();
                cufftExecZ2Z(ft_plan,ft_wave_function_d,Phi_dd_d,CUFFT_INVERSE);
                cudaDeviceSynchronize();
                SimpleKernels::rescale<<<gridSize,blockSize>>>(Phi_dd_d,1./npoints,npoints);
                cudaDeviceSynchronize();

                // Solve step-1 of operator splitting, i.e. the one NOT involving Fourier transforms
                SolverKernels::step_1_operator_splitting_dipolars<<<gridSize,blockSize>>>(wave_function_d,
                                                                                            external_potential_d,
                                                                                            Phi_dd_d,
                                                                                            time_step_d,
                                                                                            scattering_length_d,
                                                                                            gamma_epsilondd_d,
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
            cudaFree(c_density_d);
        }

        /**
         * @brief Real-time operator splitting for Faraday patterns
         */

        void DipolarGPSolver::run_operator_splitting_faraday(int number_of_time_steps,
                                                             double time_step,
                                                             double modulation_amplitude,
                                                             double modulation_frequency,
                                                             std::ostream &output_stream,
                                                             int write_output_every)
        {

            // Copy input data into the device
            cudaMallocManaged(&time_step_d,sizeof(double));
            cudaMemcpy(time_step_d,&time_step,sizeof(double),cudaMemcpyHostToDevice);

            // Initialize the fft plan required for the calculation of the laplacian
            cufftHandle ft_plan;
            if(problem_is_2d)
                cufftPlan2d(&ft_plan,nx,ny,CUFFT_Z2Z);
            else if(problem_is_3d)
                cufftPlan3d(&ft_plan,nx,ny,nz,CUFFT_Z2Z);
            cuDoubleComplex* c_density_d;
            cudaMalloc(&c_density_d,npoints*sizeof(cuDoubleComplex));

            // Initialize other variables
            this->write_output_every=write_output_every;
            double initial_scattering_length;
            double current_scattering_length;
            cudaMemcpy(&initial_scattering_length,scattering_length_d,sizeof(double),cudaMemcpyDeviceToHost);

            //----------------------------------------------------//
            //    Here the operator-splitting iterations start    //
            //----------------------------------------------------//
            for (size_t it = 0; it < number_of_time_steps; ++it)
            {

                // Write output starting from the very first iteration
                if(it % write_output_every == 0)
                    write_operator_splitting_output(it,output_stream);

                // Calculate the current value of dipolar potential
                SimpleKernels::square_vector<<<gridSize,blockSize>>>(c_density_d,wave_function_d,npoints);
                cudaDeviceSynchronize();
                cufftExecZ2Z(ft_plan,c_density_d,ft_wave_function_d,CUFFT_FORWARD);
                cudaDeviceSynchronize();
                SimpleKernels::vector_multiplication<<<gridSize,blockSize>>>(ft_wave_function_d,
                        Vtilde_d,
                        npoints);
                cudaDeviceSynchronize();
                cufftExecZ2Z(ft_plan,ft_wave_function_d,Phi_dd_d,CUFFT_INVERSE);
                cudaDeviceSynchronize();
                SimpleKernels::rescale<<<gridSize,blockSize>>>(Phi_dd_d,1./npoints,npoints);
                cudaDeviceSynchronize();

                // Solve step-1 of operator splitting, i.e. the one NOT involving Fourier transforms
                double current_time = it*time_step;
                current_scattering_length =
                        initial_scattering_length*(1.0+modulation_amplitude*std::cos(TWOPI*modulation_frequency*current_time));
                cudaMemcpy(scattering_length_d,&current_scattering_length,sizeof(double),cudaMemcpyHostToDevice);

                SolverKernels::step_1_operator_splitting_dipolars<<<gridSize,blockSize>>>(wave_function_d,
                        external_potential_d,
                        Phi_dd_d,
                        time_step_d,
                        scattering_length_d,
                        gamma_epsilondd_d,
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
            cudaFree(c_density_d);
        }

        /**
         * @brief Operator splitting output.
         *
         * This function is called after a copy of the current wave function outside of the GPU, is such a way that it
         * can be used for example for data analysis or to write it to a file for visualization. Since each call
         * blocks the real-time evolution on the GPU until the function has finished, it is better to use it with
         * moderation to avoid a big loss of performance.
         *
         * */

        void DipolarGPSolver::write_operator_splitting_output(size_t        iteration_number,
                                                              std::ostream& output_stream)
        {}

        /**
         * @brief Copy the wave function out from device to host
         *
         * In derived classes, the wave function will be available as wave_function_output
         *
         * */

        void DipolarGPSolver::copy_out_wave_function()
        {
            cudaMemcpy(wave_function_output.data(),
                       wave_function_d,
                       npoints*sizeof(cuDoubleComplex),
                       cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
        }

        /**
         * @brief Reinitialize the solver with new external potential and wave function
         *
         * */

        void DipolarGPSolver::reinit(Vector<double> &Vext,Vector<std::complex<double>> &psi)
        {
            cudaMemcpy(wave_function_d,psi.data(),npoints*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
            cudaMemcpy(external_potential_d,Vext.data(),npoints*sizeof(double),cudaMemcpyHostToDevice);
        }

        /**
         *
         * @brief Reinitialize the solver with new external potential, wave function and scattering length.
         *
         * */

        void DipolarGPSolver::reinit(Vector<double> &Vext,Vector<std::complex<double>> &psi,double scattering_length)
        {
            cudaMemcpy(wave_function_d,psi.data(),npoints*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
            cudaMemcpy(external_potential_d,Vext.data(),npoints*sizeof(double),cudaMemcpyHostToDevice);
            cudaMemcpy(scattering_length_d,&scattering_length,sizeof(double),cudaMemcpyHostToDevice);
        }

        /**
         * @brief Set initial conditions for a Truncated Wigner run
         *
         */

        void DipolarGPSolver::set_tw_initial_conditions(bool system_is_trapped)
        {

            // Need to get in a copy of the scattering length
            double scattering_length;
            cudaMemcpy(&scattering_length,scattering_length_d,sizeof(double),cudaMemcpyDeviceToHost);

            // First, consider the case in which the system is not trapped

            if (!system_is_trapped)
            {

                // Obtain a random seed from the clock
                std::default_random_engine generator;
                typedef std::chrono::high_resolution_clock clock;
                clock::time_point beginning = clock::now();
                clock::duration d = clock::now() - beginning;
                generator.seed(d.count());
                std::uniform_real_distribution<double> distribution(-1, 1);

                // Generate the alphas
                double u, v, s; // useful additional variables that we use to generate our random numbers
                std::complex<double> alphak;

                // Now refill the initial wave function with the single particle modes
                if (problem_is_2d)
                {

                    Vector<std::complex<double>> psitilde_tw(nx, ny);
                    Vector<std::complex<double>> psi(nx, ny);
                    MKLWrappers::DFtCalculator dft_tw(psi, psitilde_tw);

                    cudaMemcpy(psi.data(), wave_function_d, npoints * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

                    dft_tw.compute_forward();

                    double density = initial_norm_d[0] / (4. * x_axis(nx - 1) * y_axis(ny - 1));
                    double eps_k;
                    std::complex<double> Ek, uk, vk;

                    for (int i = 0; i < nx; ++i)
                        for (int j = 0; j < ny; ++j)
                        {
                            do
                            {
                                u = distribution(generator);
                                v = distribution(generator);
                                s = u * u + v * v;
                            } while (s >= 1.0 || s == 0);

                            s = sqrt((-2.0 * log(s)) / s);
                            u = u * s;
                            v = v * s;
                            alphak.real(std::sqrt(0.25) * u);
                            alphak.imag(std::sqrt(0.25) * v);

                            eps_k = 0.5 * std::sqrt(pow(kx_axis[i], 2) + pow(ky_axis[j], 2));
                            Ek = std::sqrt(eps_k*(eps_k+2*density*(4*PI*scattering_length+Vtilde(i, j))));
                            if (eps_k == 0)
                            {
                                uk = 1;
                                vk = 0;
                            }
                            else
                            {
                                uk = 0.5 * (sqrt(eps_k / Ek) + sqrt(Ek / eps_k));
                                vk = 0.5 * (sqrt(eps_k / Ek) - sqrt(Ek / eps_k));
                            }

                            psitilde_tw(i, j) += (alphak * uk - conj(alphak) * vk);

                        }

                    dft_tw.compute_backward();

                    cudaMemcpy(wave_function_d,psi.data(),npoints*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);

                }
                else if (problem_is_3d)
                {

                    Vector<std::complex<double>> psitilde_tw(nx,ny,nz);
                    Vector<std::complex<double>> psi(nx,ny,nz);
                    MKLWrappers::DFtCalculator dft_tw(psi,psitilde_tw);

                    cudaMemcpy(psi.data(), wave_function_d, npoints * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

                    dft_tw.compute_forward();

                    double density = initial_norm_d[0] / (4. * x_axis(nx - 1) * y_axis(ny - 1) * z_axis(nz-1));
                    double eps_k;
                    std::complex<double> Ek, uk, vk;

                    for (int i = 0; i < nx; ++i)
                        for (int j = 0; j < ny; ++j)
                            for(int k = 0; k < nz; ++k)
                            {
                                do
                                {
                                    u = distribution(generator);
                                    v = distribution(generator);
                                    s = u * u + v * v;
                                } while (s >= 1.0 || s == 0);

                                s = sqrt((-2.0 * log(s)) / s);
                                u = u * s;
                                v = v * s;
                                alphak.real(std::sqrt(0.25) * u);
                                alphak.imag(std::sqrt(0.25) * v);

                                eps_k = 0.5 * std::sqrt(pow(kx_axis[i], 2) + pow(ky_axis[j], 2) + pow(kz_axis[k],2));
                                Ek = std::sqrt(eps_k*(eps_k+2*density*(4*PI*scattering_length+Vtilde(i, j, k))));
                                if (eps_k == 0)
                                {
                                    uk = 1;
                                    vk = 0;
                                }
                                else
                                {
                                    uk = 0.5 * (sqrt(eps_k / Ek) + sqrt(Ek / eps_k));
                                    vk = 0.5 * (sqrt(eps_k / Ek) - sqrt(Ek / eps_k));
                                }

                                psitilde_tw(i,j,k) += (alphak * uk - conj(alphak) * vk);

                            }

                    dft_tw.compute_backward();

                    cudaMemcpy(wave_function_d,psi.data(),npoints*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);

                }

                // We are now ready to perform a new run of the TWA.
            }
        }

        /**
          * @brief Evaluate integrated momentum distributions for a three-dimensional problem.
          *
          * */

        void DipolarGPSolver::evaluate_integrated_occupation_number()
        {
            // If this is the first time that the function is called, it needs to set up accordingly some internal
            // control variable, plus it needs to initialize some fata
            if(first_call_evaluate_integrated_occupation_number)
            {
                first_call_evaluate_integrated_occupation_number = false;
                cufftPlan3d(&dft_handle_integrated_occupation_number,nx,ny,nz,CUFFT_Z2Z);
                integrated_occupation_number.resize(3);
                integrated_occupation_number[0].reinit(nx);
                integrated_occupation_number[1].reinit(ny);
                integrated_occupation_number[2].reinit(nz);
                cudaMalloc(&xlayer_d,ny*nz*sizeof(double));
                cudaMalloc(&ylayer_d,nx*nz*sizeof(double));
                cudaMalloc(&zlayer_d,nx*ny*sizeof(double));
                cudaMallocManaged(&x_integral,sizeof(double));
                cudaMallocManaged(&y_integral,sizeof(double));
                cudaMallocManaged(&z_integral,sizeof(double));
                cudaMalloc(&occupation_number_d,nx*ny*nz*sizeof(double));
            }

            // Execute the FFT
            cufftExecZ2Z(dft_handle_integrated_occupation_number,
                         wave_function_d,
                         ft_wave_function_d,
                         CUFFT_FORWARD);
            cudaDeviceSynchronize();

            calculate_density(occupation_number_d,ft_wave_function_d,npoints);
            cudaDeviceSynchronize();

            // Calculate integrated occupation number along x, i.e. \int dydz |\psi_k(x,y,z)|^2
            dim3 threadsPerBlock_x(16,16);
            dim3 numBlocks_x((ny + threadsPerBlock_x.x -1) / threadsPerBlock_x.x,
                             (nz + threadsPerBlock_x.y -1) / threadsPerBlock_x.y);
            for(int layer_number = 0; layer_number < nx; ++layer_number)
            {
                SimpleKernels::extract_layer_x<<<threadsPerBlock_x,numBlocks_x>>>(xlayer_d,
                                                                                  occupation_number_d,
                                                                                  layer_number,
                                                                                  nx,ny,nz);
                cudaDeviceSynchronize();

                x_integral[0] = 0.0;
                cub::DeviceReduce::Sum(temporary_storage_d,size_temporary_storage,xlayer_d,x_integral,ny*nz);
                cudaDeviceSynchronize();
                integrated_occupation_number[0][layer_number] = x_integral[0];
            }

            // Calculate integrated occupation number along y, i.e. \int dxdz |\psi_k(x,y,z)|^2
            dim3 threadsPerBlock_y(16,16);
            dim3 numBlocks_y((nx + threadsPerBlock_y.x -1) / threadsPerBlock_y.x,
                             (nz + threadsPerBlock_y.y -1) / threadsPerBlock_y.y);
            for(int layer_number = 0; layer_number < ny; ++layer_number)
            {
                SimpleKernels::extract_layer_y<<<threadsPerBlock_y,numBlocks_y>>>(ylayer_d,
                                                                                  occupation_number_d,
                                                                                  layer_number,
                                                                                  nx,ny,nz);
                cudaDeviceSynchronize();

                y_integral[0] = 0.0;
                cub::DeviceReduce::Sum(temporary_storage_d,size_temporary_storage,ylayer_d,y_integral,nx*nz);
                cudaDeviceSynchronize();
                integrated_occupation_number[1][layer_number] = y_integral[0];
            }

            // Calculate integrated occupation number along z, i.e. \int dxdy |\psi_k(x,y,z)|^2
            dim3 threadsPerBlock_z(16,16);
            dim3 numBlocks_z((nx + threadsPerBlock_z.x -1) / threadsPerBlock_z.x,
                             (ny + threadsPerBlock_z.y -1) / threadsPerBlock_z.y);
            for(int layer_number = 0; layer_number < nz; ++layer_number)
            {
                SimpleKernels::extract_layer_z<<<threadsPerBlock_z,numBlocks_z>>>(zlayer_d,
                                                                                  occupation_number_d,
                                                                                  layer_number,
                                                                                  nx,ny,nz);
                cudaDeviceSynchronize();

                z_integral[0] = 0.0;
                cub::DeviceReduce::Sum(temporary_storage_d,size_temporary_storage,zlayer_d,z_integral,nx*ny);
                cudaDeviceSynchronize();
                integrated_occupation_number[2][layer_number] = z_integral[0];
            }

        }

        /**
         * @brief Apply a low-pass filter on the wave function
         *
         * This function is useful for example for visualizing vortex tangles
         *
         * */

        void DipolarGPSolver::apply_momentum_cutoff(double kc)
        {

            if(first_call_low_pass_filter)
            {
                first_call_low_pass_filter = false;
                cufftPlan3d(&dft_quasi_condensate,nx,ny,nz,CUFFT_Z2Z);
                cudaMalloc(&wave_function_quasi_condensate_d,npoints*sizeof(cuDoubleComplex));
                cudaMalloc(&ft_wave_function_quasi_condensate_d,npoints*sizeof(cuDoubleComplex));
                wave_function_quasi_condensate.reinit(nx,ny,nz);
            }

            SimpleKernels::low_pass_filter<<<gridSize,blockSize>>>(ft_wave_function_quasi_condensate_d,
                                                                   ft_wave_function_d,
                                                                   kmod2_d,
                                                                   kc,
                                                                   npoints);
            cudaDeviceSynchronize();
            cufftExecZ2Z(dft_quasi_condensate,
                         ft_wave_function_quasi_condensate_d,
                         wave_function_quasi_condensate_d,
                         CUFFT_INVERSE);
            cudaDeviceSynchronize();

            SimpleKernels::rescale<<<gridSize,blockSize>>>(wave_function_quasi_condensate_d,1./npoints,npoints);

            cudaMemcpy(wave_function_quasi_condensate.data(),
                       wave_function_quasi_condensate_d,
                       npoints*sizeof(cuDoubleComplex),
                       cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();

        }

        /**
         * @brief Apply a density threshold on a wave function.
         *
         * This is useful to eventually isolate a vortex tangle.
         *
         * */

        void DipolarGPSolver::calculate_vortex_tangle_length(double n_threshold)
        {
            if(first_call_calculate_vortex_tangle_length)
            {
                first_call_calculate_vortex_tangle_length=false;
                cudaMalloc(&vortex_tangle_density_d,npoints*sizeof(double));
                vortex_tangle_density.reinit(nx,ny,nz);
                cudaMallocManaged(&vortex_tangle_length,sizeof(double));
                file_vortex_tangle_length.open("vortex_tangle_length.txt",std::ios_base::app);
            }

            // Apply density threshold
            SimpleKernels::density_threshold<<<gridSize,blockSize>>>(wave_function_quasi_condensate_d,
                                                                     vortex_tangle_density_d,
                                                                     n_threshold,
                                                                     npoints);
            cudaDeviceSynchronize();

            // Eventually copy out the full vortex tangle
            cudaMemcpy(vortex_tangle_density.data(),
                       vortex_tangle_density_d,
                       npoints*sizeof(double),
                       cudaMemcpyDeviceToHost);

            // Calculate the total length of the tangle
            vortex_tangle_length[0] = 0.0;
            cub::DeviceReduce::Sum(temporary_storage_d,size_temporary_storage,vortex_tangle_density_d,
                                   vortex_tangle_length,ny*nz*nz);
            cudaDeviceSynchronize();
            file_vortex_tangle_length << vortex_tangle_length_time_iteration << " "
                                      << vortex_tangle_length[0]/npoints << std::endl;
            vortex_tangle_length_time_iteration += 1;
        }
    }
}