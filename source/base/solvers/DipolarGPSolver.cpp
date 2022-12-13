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

#include "GPSolvers.hpp"

#define PI 3.1415926535897932384626433
#define TWOPI (2*PI)

using namespace std;

namespace UltraCold
{
    namespace GPSolvers
    {

        /**
         * @brief Constructor for a DipolarGPSolver in two space dimensions
         *
         * @param x *Vector<double>* representing the x-axis of the Cartesian frame on which the Gross-Pitaevskii
         * equation in two space dimensions will be solved
         * @param y *Vector<double>* representing the y-axis of the Cartesian frame on which the Gross-Pitaevskii
         * equation in two space dimensions will be solved
         * @param psi_0 *Vector<complex<double>>* representing the initial wave function
         * @param Vext *Vector<double>* representing the external potential.
         * @param scattering_length *double* the scattering length in appropriate units
         * @param dipolar_length *double* the dipolar length in appropriate units
         * @param alpha *double* angle between the the polarization direction and the z-axis
         *
         */

        DipolarGPSolver::DipolarGPSolver(Vector<double>& x,
                                         Vector<double>& y,
                                         Vector<complex<double>>& psi_0,
                                         Vector<double>& Vext,
                                         double scattering_length,
                                         double dipolar_length,
                                         double alpha)
        {

            // Get the input data.

            this->x    = x;
            this->y    = y;
            this->Vext = Vext;
            this->scattering_length = scattering_length;
            this->psi  = psi_0;

            // Check the dimensions of the Vectors provided are consistent

            nx = this->x.extent(0);
            ny = this->y.extent(0);
            assert(x.order()==1);
            assert(y.order()==1);
            assert(psi.order() == 2);
            assert(Vext.order() == 2);
            assert(psi.extent(0) == nx);
            assert(psi.extent(1) == ny);
            assert(Vext.extent(0) == nx);
            assert(Vext.extent(1) == ny);

            // Initialize the mesh in Fourier space

            kx.reinit(nx);
            ky.reinit(ny);
            kmod2.reinit(nx,ny);
            create_mesh_in_Fourier_space(this->x,this->y,kx,ky);

            for (size_t i = 0; i < nx; ++i)
                for (size_t j = 0; j < ny; ++j)
                        kmod2(i,j) = pow(kx(i),2) +
                                     pow(ky(j),2);

            // Initialize hpsi and psitilde

            hpsi.reinit(nx,ny);
            psitilde.reinit(nx,ny);

            // Initialize space steps

            dx = this->x(1)-this->x(0);
            dy = this->y(1)-this->y(0);
            dv = dx*dy;

            // Get the norm of the initial wave function
            initial_norm = 0.0;
            for (size_t i = 0; i < this->psi.size(); ++i)
                initial_norm += std::norm(this->psi[i]);
            initial_norm *= dv;

            //////////////////////////////////////////////////////////////////
            // Initialize the Fourier transform of the dipolar potential
            //////////////////////////////////////////////////////////////////

            epsilon_dd = 0.0;
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

            Phi_dd.reinit(nx,ny);
            Phi_tilde.reinit(nx,ny);

            gamma_epsilon_dd = 0.0; // No LHY in 2D

        // Scattering length is divided by sqrt(2PI) here, since in the propagators it is multiplied by 4PI
        this->scattering_length *= 1./sqrt(2*PI);

        system_is_2d = true;

    }

        /**
         * @brief Constructor for a DipolarGPSolver in three space dimensions
         *
         * @param x *Vector<double>* representing the x-axis of the Cartesian frame on which the Gross-Pitaevskii
         * equation in two space dimensions will be solved
         * @param y *Vector<double>* representing the y-axis of the Cartesian frame on which the Gross-Pitaevskii
         * equation in two space dimensions will be solved
         * @param z *Vector<double>* representing the z-axis of the Cartesian frame on which the Gross-Pitaevskii
         * equation in two space dimensions will be solved
         * @param psi_0 *Vector<complex<double>>* representing the initial wave function
         * @param Vext *Vector<double>* representing the external potential.
         * @param scattering_length *double* the scattering length in appropriate units
         * @param dipolar_length *double* the dipolar length in appropriate units
         * @param theta_mu *double* angle between the the polarization direction and the z-axis
         * @param phi_mu *double* angle between the projection of the polarization direction in the x-y plane and
         * the x-axis
         *
         */

        DipolarGPSolver::DipolarGPSolver(Vector<double>& x,
                                         Vector<double>& y,
                                         Vector<double>& z,
                                         Vector<complex<double>>& psi_0,
                                         Vector<double>& Vext,
                                         double scattering_length,
                                         double dipolar_length,
                                         double theta_mu,
                                         double phi_mu)
        {

            // Get the input data.

            this->x    = x;
            this->y    = y;
            this->z    = z;
            this->Vext = Vext;
            this->scattering_length = scattering_length;
            this->psi  = psi_0;

            // Check the dimensions of the Vectors provided are consistent

            nx = this->x.extent(0);
            ny = this->y.extent(0);
            nz = this->z.extent(0);
            assert(x.order()==1);
            assert(y.order()==1);
            assert(z.order()==1);
            assert(psi.order() == 3);
            assert(Vext.order() == 3);
            assert(psi.extent(0) == nx);
            assert(psi.extent(1) == ny);
            assert(psi.extent(2) == nz);
            assert(Vext.extent(0) == nx);
            assert(Vext.extent(1) == ny);
            assert(Vext.extent(2) == nz);

            // Initialize the mesh in Fourier space

            kx.reinit(nx);
            ky.reinit(ny);
            kz.reinit(nz);
            kmod2.reinit(nx,ny,nz);
            create_mesh_in_Fourier_space(this->x,this->y,this->z,kx,ky,kz);

            for (size_t i = 0; i < nx; ++i)
                for (size_t j = 0; j < ny; ++j)
                    for (size_t k = 0; k < nz; ++k)
                        kmod2(i,j,k) = pow(kx(i),2) +
                                       pow(ky(j),2) +
                                       pow(kz(k),2);

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

            //////////////////////////////////////////////////////////////
            // Initialize the Fourier transform of the dipolar potential
            //////////////////////////////////////////////////////////////

            epsilon_dd = 0.0;
            if(scattering_length != 0)
                epsilon_dd = dipolar_length/scattering_length;

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
                            Vtilde(i, j, k) = -4*PI*scattering_length*epsilon_dd;
                        else
                            Vtilde(i, j, k) =
                                    12.0 * PI * scattering_length * epsilon_dd * (pow(aux/aux1,2)-1.0/3.0);
                    }

            Phi_dd.reinit(nx,ny,nz);
            Phi_tilde.reinit(nx,ny,nz);

            // Initialize gamma(\epsilon_dd) for the LHY correction

            gamma_epsilon_dd = 0.0;

            if (epsilon_dd != 0)
            {
                gamma_epsilon_dd = 64.0*sqrt(PI)/3*sqrt(pow(scattering_length,5));
                double F_epsilon_dd=0.0;
                int n_theta=1000;
                double d_theta=PI/(n_theta-1);

                complex<double> csum=(0.0,0.0);
                complex<double> caux;
                for (int i = 0; i < n_theta; ++i)
                {
                    double theta=i*d_theta;
                    caux = pow(1.0+epsilon_dd*(3.0*pow(cos(theta),2)-1.0),5);
                    caux = sqrt(caux);
                    csum += sin(theta)*caux;
                }

                csum *= d_theta;
                F_epsilon_dd = csum.real();
                gamma_epsilon_dd *= F_epsilon_dd;
            }

            system_is_3d=true;

        }

        /**
         * @brief Reinitialize the solver with new external potential and initial condition.
         * @param Vext *Vector<double>* the new external potential.
         * @param psi_0 *Vector<complex<double>>* the new initial wave function.
         */

        void DipolarGPSolver::reinit(Vector<double>& Vext,
                                     Vector<complex<double>>& psi_0)
        {
            assert(psi_0.order() == Vext.order());
            assert(psi_0.extent(0) == Vext.extent(0));
            assert(psi_0.extent(1) == Vext.extent(1));
            assert(psi_0.extent(2) == Vext.extent(2));
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
         * @param output_stream *ostream* stream to which eventual text output can be passed
         * @param write_output_every *int* number of iterations after which some output must be written
         *
         * \return tuple<Vector<complex<double>>,double> representing the calculated ground-state wave
         * function and chemical potential. Can be recovered by using tie(psi,chemical_potential).
         *
         * This member function allows to solve the extended Gross-Pitaevskii equation describing a dipolar Bose gas,
         * with the addition of the Lee-Huang-Yang correction, using the same gradient-descend algorithm, accelerated
         * via the heavy-ball method, used in the corresponding member function of the class GPSolver. See that
         * documentation for more details.
         *
         * */

        tuple<Vector<complex<double>>,double>
                DipolarGPSolver::run_gradient_descent(int max_num_iter,
                                                      double tolerance,
                                                      double alpha,
                                                      double beta,
                                                      ostream& output_stream,
                                                      int write_output_every)
        {

            // Useful additional variables
            size_t sz = psi.size();               // just a useful shortcut
            Vector<complex<double>> psi_old(psi); // value of psi at previous step for heavy ball
            Vector<complex<double>> psi_new(psi); // value of psi at next step for heavy ball
            this->write_output_every = write_output_every;

            // Since the G.P. equation is solved on a cartesian mesh with periodic boundary conditions, a
            // DFtCalculator is needed to calculate the laplacian of psi
            MKLWrappers::DFtCalculator dft(hpsi,psitilde);

            // For a dipolar system, one actually needs also an additional DFtCalculator for the dipolar potential
            size_t sz2 = Phi_tilde.size(); // Another useful shortcut
            MKLWrappers::DFtCalculator dft_dd(Phi_dd,Phi_tilde);

            //    Here the gradient-descent iterations start
            size_t it;

            for (it = 0; it < max_num_iter; ++it)
            {

                // Copy the data of psi into hpsi
#pragma omp parallel for
                for (size_t i = 0; i < sz; ++i)
                    hpsi[i] = psi[i];

                // Calculate -0.5*\nabla^2\psi
                dft.compute_forward();
#pragma omp parallel for
                for (size_t i = 0; i < sz; ++i)
                    psitilde[i] *= 0.5*pow(TWOPI,2)*kmod2[i];
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
                                    + gamma_epsilon_dd*pow(abs(psi[i]),3)
                                    ) * psi[i];

                // Calculate chemical potential
                chemical_potential = 0.0;
                norm = 0.0;
#pragma omp for reduction (+: chemical_potential,norm)
                for (size_t i = 0; i < sz; ++i)
                {
                    chemical_potential += real(conj(psi[i]*hpsi[i]));
                    norm += std::norm(psi[i]);
                }
                chemical_potential = chemical_potential/norm;

                // Calculate the residual and check if convergence has been reached with the requested accuracy
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
                        psi_new[i] = (1. + beta )*psi[i] - alpha*hpsi[i] - beta*psi_old[i];

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
                        psi[i] = sqrt(initial_norm/norm)*psi_new[i];

                    // Write some output
                    if(it%write_output_every==0)
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
                        endl;
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
                        endl;
            }

            return make_pair(psi,chemical_potential);

        }

        /**
         * @brief Write output at each step of the gradient descent iterations. This function can
         * (and should) be overridden in derived classes.
         * @param iteration_number *size_t* current iteration number
         * @param output_stream *ostream* stream to which eventual text output can be passed
         *
         * \warning Writing output data files at each gradient descent step may be useful,
         * but it is also very expensive. Override this member function with care!
         *
         */

        void DipolarGPSolver::write_gradient_descent_output(size_t iteration_number,
                                                            ostream & output_stream)
        {
            output_stream
                    << iteration_number << " " << residual << " " << chemical_potential
                    <<
                    endl;
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
                psi[i] *= exp(-ci*time_step*(Vext[i]
                                            + 4.0*PI*scattering_length * std::norm(psi[i])
                                            + Phi_dd[i]
                                            + gamma_epsilon_dd*pow(abs(psi[i]),3))
                                );
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
                psitilde[i] *= exp(-ci*time_step*0.5*pow(TWOPI,2)*kmod2[i]);
            dft_calculator_step_2.compute_backward();

        }

        /**
         *
         * @brief Solve the extended Gross-Pitaevskii equation using simple operator splitting.
         * @param number_of_time_steps *int* The total number of time-steps to be performed
         * @param time_step *double* time step in appropriate units
         * @param output_stream *ostream* stream to which eventual text output can be passed
         * @param write_output_every *int* number of iterations after which some output must be written
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
                                                     ostream& output_stream,
                                                     int write_output_every)
        {

            this->time_step = time_step;
            this-> write_output_every = write_output_every;

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

                if(it%write_output_every == 0)
                    write_operator_splitting_output(it, output_stream);

                // Solve the ODE for step-1

                solve_step_1_operator_splitting(dft_dd);

                // Solve the ODE for step-2

                solve_step_2_operator_splitting(dft_calculator_step_2);

            }
        }

        /**
         * @brief Useful possible overload
         *
        */

        void DipolarGPSolver::run_operator_splitting(int, double, double, ostream&,int){};

        /**
         * @brief Useful possible overload
         *
        */

        void DipolarGPSolver::run_operator_splitting(int, double, double, double, ostream&,int){};

        /**
         * @brief Useful possible overload
         *
        */

        void DipolarGPSolver::run_operator_splitting(int, double, double, double, double, ostream&,int){};

        /**
         * @brief Write some output at each time step
         * @param iteration_number *size_t* current iteration number
         * @param output_stream *ostream* stream to which eventual text output can be passed
         *
         */

        void DipolarGPSolver::write_operator_splitting_output(size_t iteration_number,
                                                              ostream& output_stream)
        {
            output_stream << iteration_number << " " << time_step*iteration_number << endl;
        }

        /**
         * @brief Useful possible overload
         *
         */

        void DipolarGPSolver::write_operator_splitting_output(size_t, double, ostream&){};

        /**
         * @brief Useful possible overload
         *
         */

        void DipolarGPSolver::write_operator_splitting_output(size_t, double, double, ostream&){};

        /**
         * @brief Useful possible overload
         *
         */

        void DipolarGPSolver::write_operator_splitting_output(size_t, double, double, double, ostream&){};

        /**
         * @brief set initial conditions for a truncated wigner run
         * */

        void DipolarGPSolver::set_tw_initial_conditions(bool system_is_trapped)
        {

            // First, consider the case in which the system is not trapped

            if(!system_is_trapped)
            {

                // Obtain a random seed from the clock
                std::default_random_engine generator;
                typedef std::chrono::high_resolution_clock clock;
                clock::time_point                          beginning = clock::now();
                clock::duration                            d         = clock::now() - beginning;
                generator.seed(d.count());
                std::uniform_real_distribution<double> distribution(-1,1);

                // Generate the alphas
                double u, v, s; // useful additional variables that we use to generate our random numbers
                std::complex<double> alphak;

                // Now refill the initial wave function with the single particle modes
                if (system_is_2d)
                {
                    Vector<std::complex<double>> psitilde_tw(nx,ny);
                    MKLWrappers::DFtCalculator dft_tw(psi,psitilde_tw);
                    dft_tw.compute_forward();

                    double density = initial_norm/(4.*x(nx-1)*y(ny-1));
                    double eps_k;
                    std::complex<double> Ek,uk,vk;

                    for(int i = 0; i < nx; ++i)
                        for(int j = 0; j < ny; ++j)
                    {
                        do
                        {
                            u = distribution(generator);
                            v = distribution(generator);
                            s = u * u + v * v;
                        }
                        while (s >= 1.0 || s == 0);

                        s = sqrt((-2.0 * log(s)) / s);
                        u = u * s;
                        v = v * s;
                        alphak.real(std::sqrt(0.25)*u);
                        alphak.imag(std::sqrt(0.25)*v);

                        eps_k = 0.5*std::sqrt(pow(kx[i],2)+pow(ky[j],2));
                        Ek = std::sqrt(eps_k*(eps_k+2*density*(4*PI*scattering_length+Vtilde(i,j))));
                        if(eps_k == 0)
                        {
                            uk = 1;
                            vk = 0;
                        }
                        else
                        {
                            uk = 0.5*(sqrt(eps_k/Ek)+sqrt(Ek/eps_k));
                            vk = 0.5*(sqrt(eps_k/Ek)-sqrt(Ek/eps_k));
                        }

                        psitilde_tw(i,j) += (alphak*uk - conj(alphak)*vk);

                    }

                    dft_tw.compute_backward();

                }
                else if (system_is_3d)
                {
                    Vector<std::complex<double>> psitilde_tw(nx,ny,nz);
                    MKLWrappers::DFtCalculator dft_tw(psi,psitilde_tw);
                    dft_tw.compute_forward();

                    double density = initial_norm/(8.*x(nx-1)*y(ny-1)*z(nz-1));
                    double eps_k;
                    std::complex<double> Ek,uk,vk;

                    for(int i = 0; i < nx; ++i)
                        for(int j = 0; j < ny; ++j)
                            for(int k = 0; k < nz; ++k)
                            {
                                do
                                {
                                    u = distribution(generator);
                                    v = distribution(generator);
                                    s = u * u + v * v;
                                }
                                while (s >= 1.0 || s == 0);

                                s = sqrt((-2.0 * log(s)) / s);
                                u = u * s;
                                v = v * s;
                                alphak.real(std::sqrt(0.25)*u);
                                alphak.imag(std::sqrt(0.25)*v);

                                eps_k = 0.5*std::sqrt(pow(kx[i],2)+pow(ky[j],2)+pow(kz[k],2));
                                Ek = std::sqrt(eps_k*(eps_k+2*density*(4*PI*scattering_length+Vtilde(i,j,k))));
                                if(eps_k == 0)
                                {
                                    uk = 1;
                                    vk = 0;
                                }
                                else
                                {
                                    uk = 0.5*(sqrt(eps_k/Ek)+sqrt(Ek/eps_k));
                                    vk = 0.5*(sqrt(eps_k/Ek)-sqrt(Ek/eps_k));
                                }

                                psitilde_tw(i,j,k) += ( alphak*uk - conj(alphak)*vk);

                            }

                    dft_tw.compute_backward();

                }

                // We are now ready to perform a new run of the TWA.
            }

        }
    }

}

