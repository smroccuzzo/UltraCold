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

#include "gp_solvers.hpp"

#define PI 3.1415926535897932384626433
#define TWOPI (2*PI)

namespace UltraCold
{
    namespace GPSolvers
    {
        /**
         * @brief Constructor for a DipolarGPSolver in three space dimensions
         *
         * @param x *Vector<double>* representing the x-axis of the Cartesian frame on which the Gross-Pitaevskii
         * equation in two space dimensions will be solved
         * @param y *Vector<double>* representing the y-axis of the Cartesian frame on which the Gross-Pitaevskii
         * equation in two space dimensions will be solved
         * @param z *Vector<double>* representing the z-axis of the Cartesian frame on which the Gross-Pitaevskii
         * equation in two space dimensions will be solved
         * @param psi_0 *Vector<std::complex<double>>* representing the initial wave function
         * @param Vext *Vector<double>* representing the external potential.
         * @param scattering_length *double* the scattering length in appropriate units
         * @param epsilon_dd *double* ratio between scattering and dipolar length
         *
         */

        DipolarGPSolver::DipolarGPSolver(Vector<double>& x,
                                         Vector<double>& y,
                                         Vector<double>& z,
                                         Vector<std::complex<double>>& psi_0,
                                         Vector<double>& Vext,
                                         double scattering_length,
                                         double dipolar_length)
        {

            // Get the input data.

            this->x    = x;
            this->y    = y;
            this->z    = z;
            this->Vext = Vext;
            this->scattering_length = scattering_length;
            this->dipolar_length    = dipolar_length;
            this->psi  = psi_0;

            // Check the dimensions of the Vectors provided are consistent

            nx = this->x.extent(0);
            ny = this->y.extent(0);
            nz = this->z.extent(0);

            if(psi.extent(0) != nx || this->Vext.extent(0) != nx ||
               psi.extent(1) != ny || this->Vext.extent(1) != ny ||
               psi.extent(2) != nz || this->Vext.extent(2) != nz )
            {
                std::cout
                        << "\n\n"
                        << "*********************************************************************************\n"
                        << "Error found in the constructor of a GPSolver for a Gross-Pitaevskii\n"
                        << "equation in three space dimensions. The dimensions of the Vectors provided as\n"
                        << "input are not consistent.\n"
                        << "Terminating the program now...\n"
                        << "*********************************************************************************\n"
                        << "\n\n"
                        <<
                        std::endl;
                exit(1);
            }

            // Initialize the mesh in Fourier space

            kx.reinit(nx);
            ky.reinit(ny);
            kz.reinit(nz);
            kmod2.reinit(nx,ny,nz);
            create_mesh_in_Fourier_space(this->x,this->y,this->z,kx,ky,kz);

            for (size_t i = 0; i < nx; ++i)
                for (size_t j = 0; j < ny; ++j)
                    for (size_t k = 0; k < nz; ++k)
                        kmod2(i,j,k) = std::pow(kx(i),2) +
                                       std::pow(ky(j),2) +
                                       std::pow(kz(k),2);

            // Initialize hpsi and psitilde

            hpsi.reinit(nx,ny,nz);
            psitilde.reinit(nx,ny,nz);

            // Initialize space steps

            dx = this->x(1)-this->x(0);
            dy = this->y(1)-this->y(0);
            dz = this->z(1)-this->z(0);
            dv = dx*dy*dz;

            // Get the norm of the initial wave function

            initial_norm = 0.0;
            for (size_t i = 0; i < this->psi.size(); ++i)
                initial_norm += std::norm(this->psi[i]);
            initial_norm *= dv;

            // Initialize the Fourier transform of the dipolar potential

            epsilon_dd = 0.0;
            if(scattering_length != 0)
                epsilon_dd = dipolar_length/scattering_length;

            Vtilde.reinit(nx,ny,nz/2+1);
            for (int i = 0; i < nx; ++i)
            {
                double aux = TWOPI * kx[i];
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz/2+1; ++k)
                    {
                        double aux1 = TWOPI * std::sqrt(std::pow(kx[i], 2) + std::pow(ky[j], 2) + std::pow(kz[k], 2));
                        if (aux1 <= 1.E-6)
                            Vtilde(i, j, k) = 0.0;
                        else
                            Vtilde(i, j, k) =
                                    12.0 * PI * scattering_length * epsilon_dd * (std::pow(aux/aux1,2)-1.0/3.0);
                    }
            }

            Phi_dd.reinit(nx,ny,nz);
            Phi_tilde.reinit(nx,ny,nz/2+1);

            // Initialize gamma(\epsilon_dd) for the LHY correction

            gamma_epsilon_dd = 0.0;

            if (epsilon_dd != 0)
            {
                gamma_epsilon_dd = 64.0*std::sqrt(PI)/3*std::sqrt(std::pow(scattering_length,5));
                double F_epsilon_dd=0.0;
                int n_theta=1000;
                double d_theta=PI/(n_theta-1);

                std::complex<double> csum=(0.0,0.0);
                std::complex<double> caux;
                for (int i = 0; i < n_theta; ++i)
                {
                    double theta=i*d_theta;
                    caux = std::pow(1.0+epsilon_dd*(3.0*std::pow(std::cos(theta),2)-1.0),5);
                    caux = std::sqrt(caux);
                    csum += std::sin(theta)*caux;
                }

                csum *= d_theta;
                F_epsilon_dd = csum.real();
                gamma_epsilon_dd *= F_epsilon_dd;
            }

        }

        /**
         * @brief Reinitialize the solver with new external potential and initial condition.
         * @param Vext *Vector<double>* the new external potential.
         * @param psi_0 *Vector<std::complex<double>>* the new initial wave function.
         *
         * \warning This function does not perform any bound check, hence you must be careful
         * to pass Vectors with the same dimensionality and extents as those passed to the constructor.
         *
         */

        void DipolarGPSolver::reinit(Vector<double>& Vext,
                                     Vector<std::complex<double>>& psi_0)
        {
            this->Vext = Vext;
            this->psi  = psi_0;
        }

        /**
         *
         * @brief Calculates a ground-state solution to the stationary, extended Gross-Pitaevskii equation with
         * Lee-Huang-Yang correction.
         *
         * @param max_num_iter *int* the maximum number of gradient descent iterations
         * @param tolerance  *double* the maximum norm of the residual, below which the algorithm is considered as
         * converged
         * @param alpha *double* the step-length of the gradient-descent iterations.
         * @param beta *double* acceleration step for the heavy-ball method.
         * @param output_stream *std::ostream* stream to which eventual text output can be passed
         *
         * \return std::tuple<Vector<std::complex<double>>,double> representing the calculated ground-state wave
         * function and chemical potential. Can be recovered by using std::tie(psi,chemical_potential).
         *
         * This member function allows to solve the extended Gross-Pitaevskii equation describing a dipolar Bose gas,
         * with the addition of the Lee-Huang-Yang correction, using the same gradient-descend algorithm, accelerated
         * via the heavy-ball method, used in the corresponding member function of the class GPSolver. See that
         * documentation for more details.
         *
         * */

        std::tuple<Vector<std::complex<double>>,double>
                DipolarGPSolver::run_gradient_descent(int max_num_iter,
                                                      double tolerance,
                                                      double alpha,
                                                      double beta,
                                                      std::ostream& output_stream)
        {

            // Useful additional variables

            size_t sz = psi.size();          // just a useful shortcut
            Vector<std::complex<double>> psi_old(psi); // value of psi at previous step for heavy ball
            Vector<std::complex<double>> psi_new(psi); // value of psi at next step for heavy ball

            // Since the G.P. equation is solved on a cartesian mesh with periodic boundary conditions, a
            // DFtCalculator is needed to calculate the laplacian of psi

            MKLWrappers::DFtCalculator dft(hpsi,psitilde);

            // For a dipolar system, one actually needs also an additional DFtCalculator for the dipolar potential

            size_t sz2 = Phi_tilde.size(); // Another useful shortcut
            MKLWrappers::DFtCalculator dft_dd(Phi_dd,Phi_tilde);

            //--------------------------------------------------//
            //    Here the gradient-descent iterations start    //
            //--------------------------------------------------//

            size_t it;
            for (it = 0; it < max_num_iter; ++it)
            {

                //-----------------------//
                //    Calculate H\psi    //
                //-----------------------//

                // Copy the data of psi into hpsi

#pragma omp parallel for
                for (size_t i = 0; i < sz; ++i)
                    hpsi[i] = psi[i];

                // Calculate -0.5*\nabla^2\psi

                dft.compute_forward();
#pragma omp parallel for
                for (size_t i = 0; i < sz; ++i)
                    psitilde[i] *= 0.5*std::pow(TWOPI,2)*kmod2[i];
                dft.compute_backward();

                // Calculate the contribution from the dipolar potential and add together with the lhy term

#pragma omp parallel for
                for (int i = 0; i < sz; ++i)
                    Phi_dd[i] = std::norm(psi[i]);
                dft_dd.compute_forward();
#pragma omp parallel for
                for (int i = 0; i < sz2; ++i)
                    Phi_tilde[i] *= Vtilde[i];
                dft_dd.compute_backward();

#pragma omp parallel for
                for (int i = 0; i < sz; ++i)
                            hpsi[i] += ( Vext[i]
                                    + 4.0*PI*scattering_length * std::norm(psi[i])
                                    + Phi_dd[i]
                                    + gamma_epsilon_dd*std::pow(std::abs(psi[i]),3)) * psi[i];

                //------------------------------------//
                //    Calculate chemical potential    //
                //------------------------------------//

                chemical_potential = 0.0;
                norm = 0.0;
#pragma omp for reduction (+: chemical_potential,norm)
                for (size_t i = 0; i < sz; ++i)
                {
                    chemical_potential += std::real(std::conj(psi[i])*hpsi[i]);
                    norm += std::norm(psi[i]);
                }
                chemical_potential = chemical_potential/norm;

                //------------------------------------------------------------------------------------------------//
                //   Calculate the residual and check if convergence has been reached with the requested accuracy //
                //------------------------------------------------------------------------------------------------//

                residual = 0.0;
#pragma omp for reduction (+: residual)
                for (size_t i = 0; i < sz; ++i)
                    residual += std::norm(hpsi[i]-chemical_potential*psi[i]);
                residual *= dv;

                if(residual < tolerance)
                {

                    // Convergence reached, exit the gradient descent loop

                    break;

                }
                else
                {
                    // Perform another gradient descent step

#pragma omp parallel for
                    for (size_t i = 0; i < sz; ++i)
                        psi_new[i] = (1. + beta)*psi[i] - alpha*hpsi[i] - beta*psi_old[i];

                    // Save the current value of psi

#pragma omp parallel for
                    for (size_t i = 0; i < sz; ++i)
                        psi_old[i] = psi[i];

                    // Normalize psi

                    norm = 0.0;
#pragma omp for reduction (+: norm)
                    for (size_t i = 0; i < sz; ++i)
                        norm += std::norm(psi_new[i]);
                    norm *= dv;

#pragma omp parallel for
                    for (size_t i = 0; i < sz; ++i)
                        psi[i] = std::sqrt(initial_norm/norm)*psi_new[i];

                    // Write some output

                    write_gradient_descent_output(it,output_stream);

                }

            }

            //----------------------------------------------------------------------//
            //    Gradient descent iterations finished, writing the final output    //
            //----------------------------------------------------------------------//

            last_iteration_number = it;
            if(last_iteration_number == max_num_iter)
            {
                output_stream
                        << "\n\n"
                        << "---------------------------------------------------------------------\n"
                        << "Maximum number of iterations hit without reaching convergence with \n"
                        << "the requested accuracy. The final value calculated for the\n"
                        << "chemical potential is " << chemical_potential << "\n"
                        << "The norm of the residual is still " << residual << "\n"
                        << "---------------------------------------------------------------------\n"
                        << "\n\n"
                        <<
                        std::endl;
            }
            else
            {
                output_stream
                        << "\n\n"
                        << "---------------------------------------------------------------------\n"
                        << "Convergence reached in " << last_iteration_number << " iterations.\n"
                        << "The final value calculated for the chemical potential is "
                        << chemical_potential << "\n"
                        << "---------------------------------------------------------------------\n"
                        << "\n\n"
                        <<
                        std::endl;
            }

            return std::make_pair(psi,chemical_potential);

        }

        /**
         * @brief Write output at each step of the gradient descent iterations. This function can
         * (and should) be overridden in derived classes.
         * @param iteration_number *size_t* current iteration number
         * @param output_stream *std::ostream* stream to which eventual text output can be passed
         *
         * \warning Writing output data files at each gradient descent step may be useful,
         * but it is also very expensive. Override this member function with care!
         *
         */

        void DipolarGPSolver::write_gradient_descent_output(size_t iteration_number,
                                                            std::ostream & output_stream)
        {
            output_stream
                    << iteration_number << " " << residual << " " << chemical_potential
                    <<
                    std::endl;

        }


        /**
         * @brief Solve step-1 of operator splitting
         *
         */

        void DipolarGPSolver::solve_step_1_operator_splitting(MKLWrappers::DFtCalculator& dft_dd)
        {

            // Calculate the contribution from the dipolar potential and add together with the lhy term

#pragma omp parallel for
            for (int i = 0; i < psi.size(); ++i)
                Phi_dd[i] = std::norm(psi[i]);
            dft_dd.compute_forward();
#pragma omp parallel for
            for (int i = 0; i < Vtilde.size(); ++i)
                Phi_tilde[i] *= Vtilde[i];
            dft_dd.compute_backward();

#pragma omp parallel for
            for (size_t i = 0; i < psi.size(); ++i)
                psi[i] *= std::exp(-ci*time_step*(Vext[i]
                        + 4.0*PI*scattering_length * std::norm(psi[i])
                        + Phi_dd[i]
                        + gamma_epsilon_dd*std::pow(std::abs(psi[i]),3)));

        }

        /**
         * @brief Useful possible overload
         *
        */

        void DipolarGPSolver::solve_step_1_operator_splitting(MKLWrappers::DFtCalculator&,double) {};

        /**
         * @brief Useful possible overload
         *
        */

        void DipolarGPSolver::solve_step_1_operator_splitting(MKLWrappers::DFtCalculator&,double, double) {};

        /**
         * @brief Useful possible overload
         *
        */

        void DipolarGPSolver::solve_step_1_operator_splitting(MKLWrappers::DFtCalculator&,double, double, double) {};

        /**
         * @brief Solve step-2 of operator splitting
         * */

        void DipolarGPSolver::solve_step_2_operator_splitting(MKLWrappers::DFtCalculator& dft_calculator_step_2)
        {

            dft_calculator_step_2.compute_forward();
#pragma omp parallel for
            for (size_t i = 0; i < psi.size(); ++i)
                psitilde[i] *= std::exp(-ci*time_step*0.5*std::pow(TWOPI,2)*kmod2[i]);
            dft_calculator_step_2.compute_backward();

        }

        /**
         *
         * @brief Solve the extended Gross-Pitaevskii equation using simple operator splitting.
         * @param number_of_time_steps *int* The total number of time-steps to be performed
         * @param time_step *double* time step in appropriate units
         * @param output_stream *std::ostream* stream to which eventual text output can be passed
         *
         * This member function solves the (a-dimensional) time-dependent extended Gross-Pitaevskii equation
         * \f[
         * i \frac{\partial \psi}{\partial t} = \left[ \frac{-\nabla^2}{2}+V_{ext}({\bf r})+g|\psi|^2 \right]\psi
         * \f]
         * for the description of a dipolar Bose gas using the classic operator splitting technique. The idea is the
         * same as the corresponding member class of GPSolver. See that documentation for details.
         *
         */

        void DipolarGPSolver::run_operator_splitting(int number_of_time_steps,
                                                     double time_step,
                                                     std::ostream& output_stream)
        {

            this->time_step = time_step;

            // Since the G.P. equation is solved on a cartesian mesh with periodic boundary conditions, a
            // DFtCalculator is needed to calculate the laplacian of psi

            MKLWrappers::DFtCalculator dft_calculator_step_2(psi,psitilde);
            MKLWrappers::DFtCalculator dft_dd(Phi_dd,Phi_tilde);

            //----------------------------------------------------//
            //    Here the operator-splitting iterations start    //
            //----------------------------------------------------//

            for (size_t it = 0; it < number_of_time_steps; ++it)
            {

                // Write outputs starting from the first time step

                write_operator_splitting_output(it, output_stream);

                // Solve the ODE for step 1)

                solve_step_1_operator_splitting(dft_dd);

                // Solve the ODE for step-2

                solve_step_2_operator_splitting(dft_calculator_step_2);

            }
        }

        /**
         * @brief Useful possible overload
         *
        */

        void DipolarGPSolver::run_operator_splitting(int, double, double, std::ostream&){};

        /**
         * @brief Useful possible overload
         *
        */

        void DipolarGPSolver::run_operator_splitting(int, double, double, double, std::ostream&){};

        /**
         * @brief Useful possible overload
         *
        */

        void DipolarGPSolver::run_operator_splitting(int, double, double, double, double, std::ostream&){};

        /**
         * @brief Write some output at each time step
         * @param iteration_number *size_t* current iteration number
         * @param output_stream *std::ostream* stream to which eventual text output can be passed
         *
         */

        void DipolarGPSolver::write_operator_splitting_output(size_t iteration_number,
                                                       std::ostream& output_stream)
        {
            output_stream << iteration_number << " " << time_step*iteration_number << std::endl;
        }

        /**
         * @brief Useful possible overload
         *
         */

        void DipolarGPSolver::write_operator_splitting_output(size_t, double, std::ostream&){};

        /**
         * @brief Useful possible overload
         *
         */

        void DipolarGPSolver::write_operator_splitting_output(size_t, double, double, std::ostream&){};

        /**
         * @brief Useful possible overload
         *
         */

        void DipolarGPSolver::write_operator_splitting_output(size_t, double, double, double, std::ostream&){};
    }
}

