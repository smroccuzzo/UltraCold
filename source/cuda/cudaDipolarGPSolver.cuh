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

#ifndef ULTRACOLD_CUDA_DIPOLAR_GP_SOLVER
#define ULTRACOLD_CUDA_DIPOLAR_GP_SOLVER

#include <vector>
#include <cufft.h>
#include <fstream>
#include "cuComplex.h"
#include "Vector.hpp"

namespace UltraCold
{

    namespace cudaSolvers
    {

        /**
         * @brief GPU-accelerated solver for a dipolar Gross-Pitaevskii equation.
         *
         * The functionality of this solver class is practically identical to the one using only CPUs.
         *
         * \note For optimal performance, the total number of points should be a multiple of 32.
         *
         * */

        class DipolarGPSolver
        {
            public:

                // Constructors
                DipolarGPSolver(Vector<double>&               x,
                                Vector<double>&               y,
                                Vector<std::complex<double>>& psi_0,
                                Vector<double>&               Vext,
                                double                        scattering_length,
                                double                        dipolar_length,
                                double                        alpha);  // 2D problems
                DipolarGPSolver(Vector<double>&               x,
                                Vector<double>&               y,
                                Vector<double>&               z,
                                Vector<std::complex<double>>& psi_0,
                                Vector<double>&               Vext,
                                double                        scattering_length,
                                double                        dipolar_length,
                                double                        theta_mu,
                                double                        phi_mu,
                                bool add_lhy_correction);  // 3D problems

                // Destructor
                ~DipolarGPSolver();

                // Re-initializers
                void reinit(Vector<double>& Vext,Vector<std::complex<double>>& psi);
                void reinit(Vector<double>& Vext,Vector<std::complex<double>>& psi,double scattering_length);

                // Calculate a ground-state solution
                std::tuple<Vector<std::complex<double>>,double> run_gradient_descent(int    max_num_iter,
                                                                                     double alpha,
                                                                                     double beta,
                                                                                     std::ostream& output_stream,
                                                                                     int write_output_every);

                // Solve the GPE in real-time using the classical operator-splitting method.
                // The intermediate PDE is solved using simply Fourier transforms
                void run_operator_splitting(int number_of_time_steps,
                                            double time_step,
                                            std::ostream& output_stream,
                                            int write_output_every);

                void run_operator_splitting_faraday(int number_of_time_steps,
                                                    double time_step,
                                                    double modulation_amplitude,
                                                    double modulation_frequency,
                                                    std::ostream& output_stream,
                                                    int write_output_every);

                void set_tw_initial_conditions(bool system_is_trapped);

                // Get a pointer to the wave function stored on the device
                cuDoubleComplex* get_wave_function_device_pointer();
                cuDoubleComplex* get_ft_wave_function_device_pointer();

            protected:

                // Post-processing

                // Write output for gradient descent
                virtual void write_gradient_descent_output(size_t iteration_number,
                                                           std::ostream& output_stream);

                // Write output for operator splitting
                virtual void write_operator_splitting_output(size_t        iteration_number,
                                                             std::ostream& output_stream);

                // Copy the wave function out from device to host
                void copy_out_wave_function();

                // Evaluate the integrated occupation numbers in momentum space for three-dimensional problems
                void evaluate_integrated_occupation_number();

                // Vector data members, needed for the calculations
                // Vectors living in device memory have the _d postfix
                double* external_potential_d;
                double* density_d;
                double* norm_d;
                double* initial_norm_d;
                cuDoubleComplex* wave_function_d;
                cuDoubleComplex* ft_wave_function_d; // This contains the Fourier transform of psi
                cuDoubleComplex* hpsi_d;             // This contains the result of \hat{H}\psi
                double* x_axis_d;
                double* y_axis_d;
                double* z_axis_d;
                double* kx_axis_d;
                double* ky_axis_d;
                double* kz_axis_d;
                double* kmod2_d;
                double* r2mod_d;
                double* chemical_potential_d;
                double* scattering_length_d;
                double* alpha_d;
                double* beta_d;
                double* time_step_d;
                void* temporary_storage_d = nullptr; // This will be initialized by cub. Needed as temporary storage
                size_t size_temporary_storage = 0;   // This will also be initialized by cub with the future size of
                                                     // the problem.
                cuDoubleComplex* Vtilde_d;
                cuDoubleComplex* Phi_dd_d;
                double* epsilon_dd_d;
                double* gamma_epsilondd_d;

                Vector<double> x_axis,y_axis,z_axis,r2mod,kx_axis,ky_axis,kz_axis;
                Vector<std::complex<double>> Vtilde;   // This contains the Fourier transform of the dipolar potential

                // Other cuda necessary data members
                int gridSize;  // Number of cuda blocks to use
                int blockSize; // size of each cuda block, i.e., number of threads per block

                // Member functions calling kernels
                void calculate_density(double* density, cuDoubleComplex* wave_function, int size);

                // Number of points along each space dimension
                int nx,ny,nz;
                int npoints;

                // Other mesh parameter
                double dx = 1.0;
                double dy = 1.0;
                double dz = 1.0;
                double dv = 1.0;

                // Other useful parameters
                bool problem_is_2d=false;
                bool problem_is_3d=false;
                bool add_lhy_correction=false;

                int write_output_every;

                // Wave function for output
                Vector<std::complex<double>> wave_function_output;

                // Wave function to return
                Vector<std::complex<double>> result_wave_function;

                // These are needed for evaluating integrated occupation numbers
                bool first_call_evaluate_integrated_occupation_number=true;
                cufftHandle dft_handle_integrated_occupation_number;
                std::vector<Vector<double>> integrated_occupation_number;
                double* xlayer_d;
                double* ylayer_d;
                double* zlayer_d;
                double* occupation_number_d;
                double* x_integral;
                double* y_integral;
                double* z_integral;

                // These are needed to apply a low-pass filter to the Fourier transform of the wave function,
                // useful for example for visualizing vortex tangles
                void apply_momentum_cutoff(double kc);
                bool first_call_low_pass_filter=true;
                cuDoubleComplex* wave_function_quasi_condensate_d;
                Vector<std::complex<double>> wave_function_quasi_condensate;
                cufftHandle dft_quasi_condensate;
                cuDoubleComplex* ft_wave_function_quasi_condensate_d;

                // These are needed for thresholding the quasi-condensate density, in order to only isolate particles
                // forming a vortex tangle
                void calculate_vortex_tangle_length(double n_threshold);
                bool first_call_calculate_vortex_tangle_length=true;
                double* vortex_tangle_density_d;
                Vector<double> vortex_tangle_density;
                double* vortex_tangle_length;
                std::fstream file_vortex_tangle_length;
                int vortex_tangle_length_time_iteration = 0;
        };
    }
}


#endif
