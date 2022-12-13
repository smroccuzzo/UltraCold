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

#ifndef ULTRACOLD_CUDA_GP_SOLVER
#define ULTRACOLD_CUDA_GP_SOLVER

#include <vector>
#include "cuComplex.h"
#include "Vector.hpp"

namespace UltraCold
{

    /**
     * @brief Namespace containing GPU-accelerated solvers for various flavors of Gross-Pitaevskii equations
     * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
     *
     * This namespace contains several GPU-accelerated solver classes for various flavors of Gross-Pitaevskii equations.
     * Such solvers are implemented in CUDA, and all the computationally intensive tasks are delegated to a GPU.
     * In order to use such classes, it is necessary of course to have an NVIDIA GPU available, as well as CUDA
     * installed.
     *
     * */

    namespace cudaSolvers
    {

        /**
         * @brief GPU-accelerated solver for a simple Gross-Pitaevskii equation
         *
         * The functionality of this solver class is practically identical to the one using only CPUs.
         *
         * \note For optimal performance, the total number of points should be a multiple of 32.
         *
         * */

        class GPSolver
        {
            public:

                // Constructors
                GPSolver(Vector<double>&               x,
                         Vector<std::complex<double>>& psi_0,
                         Vector<double>&               Vext,
                         double                        scattering_length);  // 1D problems
                GPSolver(Vector<double>&               x,
                         Vector<double>&               y,
                         Vector<std::complex<double>>& psi_0,
                         Vector<double>&               Vext,
                         double                        scattering_length);  // 2D problems
                GPSolver(Vector<double>&               x,
                         Vector<double>&               y,
                         Vector<double>&               z,
                         Vector<std::complex<double>>& psi_0,
                         Vector<double>&               Vext,
                         double                        scattering_length);  // 3D problems

                // Destructor
                ~GPSolver();

                // Re-initializers
                void reinit(Vector<std::complex<double>>& psi, Vector<double>& Vext);

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

                // Get a pointer to the wave function stored on the device

                cuDoubleComplex* get_wave_function_device_pointer();
                cuDoubleComplex* get_ft_wave_function_device_pointer();

            protected:

                // Write the output
                void write_gradient_descent_output(int it);
                void write_operator_splitting_output(int it,std::ostream& output_stream);

            private:

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

                Vector<double> x_axis,y_axis,z_axis,r2mod,kx_axis,ky_axis,kz_axis;

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

                int last_iteration_number;
                bool problem_is_1d=false;
                bool problem_is_2d=false;
                bool problem_is_3d=false;

                // Eigenstates of the harmonic oscillator, useful for TWA
                std::vector<Vector<double>> eigenstates_harmonic_oscillator;

                // Wave function to return
                Vector<std::complex<double>> result_wave_function;

                // Wave function to output
                Vector<std::complex<double>> wave_function_output;

        };
    }
}


#endif
