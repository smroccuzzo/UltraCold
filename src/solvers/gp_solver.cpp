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
         * @brief Constructor for a GPSolver in one space dimension
         *
         * @param x *Vector<double>* representing the cartesian axis on which the Gross-Pitaevskii
         * equation in one space dimension will be solved
         * @param psi_0 *Vector<std::complex<double>>* representing the initial wave function
         * @param Vext *Vector<double>* representing the external potential.
         * @param scattering_length *double* the scattering length in appropriate units
         *
         */

        GPSolver::GPSolver(Vector<double>& x,
                           Vector<std::complex<double>>& psi_0,
                           Vector<double>& Vext,
                           double scattering_length)
        {

            // Get necessary copies of input data

            this->x    = x;
            this->Vext = Vext;
            this->scattering_length = scattering_length;
            this->psi  = psi_0;


            // Check the dimensions of the Vectors provided are consistent

            nx = this->x.extent(0);
            ny = 0;
            nz = 0;
            if(psi.extent(0) != nx || this->Vext.extent(0) != nx)
            {
                std::cout
                    << "\n\n"
                    << "**************************************************************************\n"
                    << "Error found in the constructor of a GPSolver for a Gross-Pitaevskii\n"
                    << "equation in one space dimension. The dimensions of the Vectors provided as\n"
                    << "input are not consistent.\n"
                    << "Terminating the program now...\n"
                    << "**************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }

            // Initialize the mesh in Fourier space

            kx.reinit(nx);
            kmod2.reinit(nx);

            create_mesh_in_Fourier_space(this->x,kx);

            for (size_t i = 0; i < nx; ++i)
                kmod2(i) = std::pow(kx(i),2);

            // Initialize hpsi and psitilde

            hpsi.reinit(nx);
            psitilde.reinit(nx);

            // Initialize space steps

            dx = this->x(1)-this->x(0);
            dv = dx;

            // Get the norm of the initial wave function

            initial_norm = 0.0;
            for (size_t i = 0; i < this->psi.size(); ++i)
                initial_norm += std::norm(this->psi[i]);
            initial_norm *= dv;

        }

        /**
         * @brief Constructor for a GPSolver in two space dimensions
         *
         * @param x *Vector<double>* representing the x-axis of the Cartesian frame on which the Gross-Pitaevskii
         * equation in two space dimensions will be solved
         * @param y *Vector<double>* representing the y-axis of the Cartesian frame on which the Gross-Pitaevskii
         * equation in two space dimensions will be solved
         * @param psi_0 *Vector<std::complex<double>>* representing the initial wave function
         * @param Vext *Vector<double>* representing the external potential.
         * @param scattering_length *double* the scattering length in appropriate units
         *
         */

        GPSolver::GPSolver(Vector<double>& x,
                           Vector<double>& y,
                           Vector<std::complex<double>>& psi_0,
                           Vector<double>& Vext,
                           double scattering_length)
        {

            // Get necessary copies of input data

            this->x    = x;
            this->y    = y;
            this->Vext = Vext;
            this->scattering_length    = scattering_length;
            this->psi  = psi_0;

            // Check the dimensions of the Vectors provided are consistent

            nx = this->x.extent(0);
            ny = this->y.extent(0);
            nz = 0;
            if(psi.extent(0) != nx || this->Vext.extent(0) != nx ||
               psi.extent(1) != ny || this->Vext.extent(1) != ny)
            {
                std::cout
                    << "\n\n"
                    << "***************************************************************************\n"
                    << "Error found in the constructor of a GPSolver for a Gross-Pitaevskii\n"
                    << "equation in two space dimensions. The dimensions of the Vectors provided as\n"
                    << "input are not consistent.\n"
                    << "Terminating the program now...\n"
                    << "***************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }

            // Initialize the mesh in Fourier space

            kx.reinit(nx);
            ky.reinit(ny);
            kmod2.reinit(nx,ny);
            create_mesh_in_Fourier_space(this->x,this->y,kx,ky);

            for (size_t i = 0; i < nx; ++i)
                for (size_t j = 0; j < ny; ++j)
                    kmod2(i,j) = std::pow(kx(i),2) +
                                 std::pow(ky(j),2);

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

        }

        /**
         * @brief Constructor for a GPSolver in three space dimensions
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
         *
         */

        GPSolver::GPSolver(Vector<double>& x,
                           Vector<double>& y,
                           Vector<double>& z,
                           Vector<std::complex<double>>& psi_0,
                           Vector<double>& Vext,
                           double scattering_length)
        {

            // Get the input data.

            this->x    = x;
            this->y    = y;
            this->z    = z;
            this->Vext = Vext;
            this->scattering_length = scattering_length;
            this->psi = psi_0;

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

        void GPSolver::reinit(Vector<double>& Vext,
                              Vector<std::complex<double>>& psi_0)
        {
            this->Vext = Vext;
            this->psi  = psi_0;
        }

        /**
         * @brief Calculates a ground-state solution to the stationary Gross-Pitaevskii equation.
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
         * The Gross-Pitaevskii equation allows to obtain information on the ground-state properties of an
         * ultra-cold bosonic system by searching for stationary solutions of the form
         * \f$ \psi({\bf r},t)=\psi_0({\bf r})e^{-i\mu t/\hbar} \f$, obtaining the stationary Gross-Pitaevskii
         * equation
         * \f[
         * \mu \psi_0 = \left[ \frac{-\hbar^2\nabla^2}{2m}+V_{ext}({\bf r})+ g|\psi_0|^2 \right]\psi_0
         * \f]
         * where \f$ g \f$ is related to the s-wave scattering length \f$ a \f$ by \f$ g=\frac{4\pi\hbar^2 a}{2m} \f$.
         * Solving this for the smallest eigenvalue \f$ \mu \f$, which represents the chemical potential,
         * gives access to the ground-state configuration of the system. \n
         * In order to solve this equation, one can notice that it can be obtained from a constrained minimization
         * formulation of the problem, in particular by requiring that the ground-state order parameter of
         * the system is the exact minimizer of the mean-field energy functional
         * \f[
         * E[\psi] = \int d{\bf r} \left[ \psi^*({\bf r})\left(\frac{-\hbar^2\nabla^2}{2m}+
         * V_{ext}({\bf r})\right)\psi({\bf r}) \right] +\frac{g}{2}\int d{\bf r} |\psi({\bf r})|^4
         * \f]
         * under the constraint of a fixed number of particles \f$ \int d{\bf r}|\psi({\bf r})|^2=N \f$. This allows
         * also to introduce the chemical potential \f$ \mu \f$ as the Lagrange multiplier fixing the total number
         * of particles. One can hence find the ground state order parameter and chemical potential using, for
         * example, a <a href="https://en.wikipedia.org/wiki/Gradient_descent"> gradient descent </a> method, which
         * is the one implemented in this function. \n
         * The idea is to start from a guess solution \f$\psi_0\f$ and generate a sequence of iterates
         * \f$ \{\psi_n\}_{\{n=0,...,\infty\}} \f$ that terminates when one is sufficiently confident to have
         * reached a (hopefully global) minimizer of the mean-field energy functional with good accuracy. In
         * particular, a good stopping criterion consists in fixing a tolerance threshold \f$\epsilon\f$
         * (which, for this function, is a user-defined input parameter) for the norm of the residual, i.e.
         * \f$ || \hat{H}\psi_n − \mu_n \psi_n ||^2 \leq \epsilon \f$, where \f$\hat{H}\f$ is the mean-field
         * Hamiltonian of the system and the estimate \f$ \mu_n \f$ of the chemical potential \f$ \mu \f$ at
         * iteration \f$ n \f$ can be calculated as
         * \f$ \mu_n = \langle \psi_n | \hat{H} | \psi_n \rangle /  \langle \psi_n | \psi_n \rangle \f$. \n
         * In deciding how to move from one iterate \f$ \psi_n \f$ to the next \f$ \psi_{n+1} \f$, line search
         * algorithms like the gradient descent method use information about the functional
         * \f$ E[\psi] \f$ at \f$ \psi_n \f$, and possibly also from earlier iterates
         * \f$ \psi_0, \psi_1, \dots, \psi_{n-1} \f$. The update criterion should be that the energy functional is
         * smaller in \f$ \psi_{n+1} \f$ then in \f$ \psi_n \f$. One thus generates a sequence
         * \f$ \psi_{n+1} = \psi_n + \alpha \chi_n \f$ such that \f$ E[\psi_{n+1}]
         * < E[\psi_n] < \dots < E[ \psi_0 ] \f$ until the stopping criterion is satisfied. The update ”direction”
         * \f$ \chi_n \f$ must thus be chosen to be a descent direction, i.e. a direction along which the functional
         * \f$ E[\psi] \f$ decreases. The step-length \f$ \alpha \f$
         * should instead be (ideally) chosen in such a way that the decrease in the energy functional is, at
         * each iteration step, the maximum possible. Since this is not, in general, an easy task, in this
         * function we accept the compromise to choose the step length \f$ \alpha \f$ empirically as an input
         * parameter at the beginning of the iteration procedure. This also implies that the user may need to run
         * the function a few times before finding an optimal value for \f$ \alpha \f$. \n
         * Coming back to the choice of \f$ \chi_n \f$, the gradient descent method consists in choosing such
         * descent direction as the opposite of the gradient of the functional \f$ E[\psi] \f$ calculated in
         * \f$ \psi_n \f$ . So, in this case, the descent direction is (minus) the functional derivative of the
         * energy functional with respect to \f$ \psi^* \f$ evaluated at \f$ \psi_n \f$ , i.e.
         * \f[
         *
         * \chi_n = - \frac{\delta E[\psi_n]}{\delta \psi^*} =
         * - \left[ \frac{-\hbar^2\nabla^2}{2m}+V_{ext}({\bf r})+g|\psi_n|^2 \right] \psi_n
         *
         * \f]
         * The ground state wave function obtained from this algorithm is normalized in such a way to preserve the
         * initial norm, i.e. the norm of the initial wave function provided. Such normalization condition, fixing
         * the \f$ L^2 \f$ norm of the ground-state wave-function \f$ \psi \f$ should be included in the iteration
         * procedure by introducing a Lagrange multiplier and minimizing the corresponding lagrangian functional.
         * In practice, it is much cheaper to normalize ”by hand” each \f$ \psi_n \f$ obtained via the gradient
         * descent iteration, by fixing
         * \f[
         *  \psi_{n+1}^{(1)} = \psi_n + \alpha \chi_n
         * \f]
         * \f[
         * \psi_{n+1} = \sqrt{\frac{N}{\int d{\bf r}|\psi_{n+1}^{(1)}|^2}}\psi_{n+1}^{(1)}
         * \f]
         * where \f$ N \f$ is the initial norm. Notice that, physically, this represents the total number of
         * particles in the system, and as such must be an integer. Hence, the initial wave function must be
         * properly normalized before the gradient-descent iterations start. \n
         * Finally, this function employs an acceleration algorithm, known as the
         * <a href="https://www.worldscientific.com/doi/abs/10.1142/S0219199700000025"> heavy ball method </a>, in
         * order to speed up the convergence of the sequence \f$ {\psi_n} \f$. The method consists in adding a
         * ”momentum” term into the gradient-descent iterations, in order to make larger steps if the descent
         * direction does not change very much, and smaller steps if it changes a lot. In practice, what one does is
         * just to modify the gradient-descent expression in
         * \f$ \psi_{n+1} = \psi_n + \alpha \chi_n + \beta(\psi_n-\psi_{n-1}) \f$
         * where, again, \f$ \beta \f$ is a parameter chosen empirically and given to this function as input.
         *
         */

        std::tuple<Vector<std::complex<double>>,double>
                GPSolver::run_gradient_descent(int max_num_iter,
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

                // Now hpsi contains the term -0.5*\nabla^2\psi. Let's add the contribution
                // from the external potential and the two-body contact interaction

#pragma omp parallel for
                for (size_t i = 0; i < sz; ++i)
                    hpsi[i] += (Vext[i] + 4.0*PI*scattering_length * std::norm(psi[i]) ) * psi[i];

                //------------------------------------//
                //    Calculate chemical potential    //
                //------------------------------------//

                chemical_potential = 0.0;
                norm = 0.0;
#pragma omp for reduction(+:chemical_potential,norm)
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
#pragma omp for reduction(+:residual)
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
#pragma omp for reduction(+:norm)
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

        void GPSolver::write_gradient_descent_output(size_t iteration_number,
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

        void GPSolver::solve_step_1_operator_splitting()
        {
#pragma omp parallel for
                for (size_t i = 0; i < psi.size(); ++i)
                    psi[i] *= std::exp(-ci*time_step*(Vext[i] + 4.0*PI*scattering_length * std::norm(psi[i])));

        }

        /**
         * @brief Useful possible overload
         *
        */
        void GPSolver::solve_step_1_operator_splitting(double) {};

        /**
         * @brief Useful possible overload
         *
        */

        void GPSolver::solve_step_1_operator_splitting(double, double) {};

        /**
         * @brief Useful possible overload
         *
        */

        void GPSolver::solve_step_1_operator_splitting(double, double, double) {};

        /**
         * @brief Solve step-2 of operator splitting
         * */

        void GPSolver::solve_step_2_operator_splitting(MKLWrappers::DFtCalculator& dft_calculator_step_2)
        {

            dft_calculator_step_2.compute_forward();
#pragma omp parallel for
            for (size_t i = 0; i < psi.size(); ++i)
                psitilde[i] *= std::exp(-ci*time_step*0.5*std::pow(TWOPI,2)*kmod2[i]);
            dft_calculator_step_2.compute_backward();

        }

        /**
         * @brief Solve the Gross-Pitaevskii equation using simple operator splitting.
         * @param number_of_time_steps *int* The total number of time-steps to be performed
         * @param time_step *double* time step in appropriate units
         * @param output_stream *std::ostream* stream to which eventual text output can be passed
         *
         * This member function solves the (a-dimensional) time-dependent Gross-Pitaevskii equation
         * \f[
         * i \frac{\partial \psi}{\partial t} = \left[ \frac{-\nabla^2}{2}+V_{ext}({\bf r})+g|\psi|^2 \right]\psi
         * \f]
         * using the classic operator splitting technique.\n
         * The general idea of operator-splitting methods is to consider an initial value problem
         * \f[
         * y' = \text{A} y + \text{B} y
         * \f]
         * where \f$ A \f$ and \f$ B \f$ are differential operators, and solve the equation considering the actions
         * of the two operators separately. There is a vast literature on operator-splitting approaches for the
         * solution of differential equations, and different methods with a higher or lower level of accuracy. In the
         * case of the Gross-Pitaevskii equation, things are even more simplified by the fact that part of the
         * method implies steps that can be solved *exactly*. \n
         * After choosing a time-step \f$ \Delta t \f$, the operator spltting scheme consists in the following steps
         * \f[
         * \begin{align}
         * & \text{For}   \quad n=0,1,\dots,\text{number_of_time_steps} \nonumber \\
         * & \text{Step 1: solve} \quad y'_1 = \text{A}y_1 \quad \text{in} \quad [t_n,t_n+\Delta t]
         * \quad \text{for} \quad y_1(0)=y(t_n) \nonumber \\
         * & \text{Step 2: solve} \quad y'_2 = \text{B}y_2 \quad \text{in} \quad [t_n,t_n+\Delta t]
         * \quad \text{for} \quad y_2(0)=y_1(t_n+\Delta t) \nonumber \\
         * & \text{Set} \quad y(t_{n+1}) = y_2(t_n+\Delta t) \nonumber \\
         * \end{align}
         * \f]
         * The steps imply the solution of an ordinary differential equation. In the case of the
         * Gross-Pitaevskii equation, a good choice of the operators \f$ A \f$ and \f$ B \f$ is the following
         * \f[
         * \begin{align}
         *  & A = V_{ext} + g |\psi(t)|^2 \nonumber \\
         *  & B = \frac{-\nabla^2}{2} \nonumber \\
         * \end{align}
         * \f]
         *
         * In this case, step 1, although involve the solution of a non-linear differential equation,
         * can be solved **analytically**. In fact, it is easy to show that the equation
         * \f[
         *  i\frac{d}{dt}\psi(t) = V_{ext} + g |\psi(t)|^2
         * \f]
         * preserves the norm \f$ \psi \f$ in time, and hence an explicit solution of this equation is
         * \f[
         * \psi(t+\Delta t) = e^{-i\Delta t\left( V_{ext} + g |\psi(t)|^2 \right)}\psi(t)
         * \f]
         * Moreover, also the second step of the method can be solved analytically, provided that one imposes
         * *periodic boundary conditions*. In fact, given the equation
         * \f[
         * i\frac{d}{dt}\psi(t) = \frac{-\nabla^2}{2}\psi(t)
         * \f]
         * and taking the Fourier transform (in space) at both sides, one finds
         * \f[
         * i\frac{d}{dt}\tilde{\psi}(t) = \frac{k^2}{2}\tilde{\psi}(t)
         * \f]
         * which can again be solved exactly as
         * \f[
         * \tilde{\psi}(t+\Delta t) = e^{-i\Delta t\frac{k^2}{2}}\tilde{\psi}(t)
         * \f]
         * Finally, an inverse Fourier transform allows to recover the solution back in real space, ready to take
         * another time step. \n
         * To summarize, this member function solves the Gross-Pitaevskii equation using classic operator
         * splitting via the following steps
         *  - Step 1) Set \f$ \psi(t_n + \Delta t) =
         *  e^{-i\Delta t\left( V_{ext} + g |\psi(t_n)|^2 \right)}\psi(t_n) \f$;
         *  - Step 2)
         *      - 2.1) Take the Fourier transform of \f$ \psi(t_n + \Delta t) \f$, and call
         *      it \f$ \tilde{\psi}(t_n) \f$;
         *      - 2.2) Set \f$ \tilde{\psi}(t_n+\Delta t) = e^{-i\Delta t\frac{k^2}{2}}\tilde{\psi}(t_n) \f$
         *      - 2.3) Take the inverse Fourier transform of \f$ \tilde{\psi}(t_n+\Delta t) \f$, and obtain
         *      \f$ \psi(t_n+\Delta t) \f$
         *
         * Notice that all these steps can be computed *locally*, meaning that we don't need additional vectors to
         * store intermediate values of \f$ \psi \f$.
         *
         */

        void GPSolver::run_operator_splitting(int number_of_time_steps,
                                              double time_step,
                                              std::ostream& output_stream)
        {

            this->time_step = time_step;

            // Since the G.P. equation is solved on a cartesian mesh with periodic boundary conditions, a
            // DFtCalculator is needed to calculate the laplacian of psi

            MKLWrappers::DFtCalculator dft_calculator_step_2(psi,psitilde);

            //----------------------------------------------------//
            //    Here the operator-splitting iterations start    //
            //----------------------------------------------------//

            for (size_t it = 0; it < number_of_time_steps; ++it)
            {

                // Write outputs starting from the first time step

                write_operator_splitting_output(it, output_stream);

                // Solve the ODE for step 1

                solve_step_1_operator_splitting();

                // Solve the ODE for step-2

                solve_step_2_operator_splitting(dft_calculator_step_2);

            }
        }

        /**
         * @brief Useful possible overload
         *
        */

        void GPSolver::run_operator_splitting(int, double, double, std::ostream&){};

        /**
         * @brief Useful possible overload
         *
        */

        void GPSolver::run_operator_splitting(int, double, double, double, std::ostream&){};

        /**
         * @brief Useful possible overload
         *
        */

        void GPSolver::run_operator_splitting(int, double, double, double, double, std::ostream&){};

        /**
         * @brief Write some output at each time step
         * @param iteration_number *size_t* current iteration number
         * @param output_stream *std::ostream* stream to which eventual text output can be passed
         *
         */

        void GPSolver::write_operator_splitting_output(size_t iteration_number,
                                                       std::ostream& output_stream)
        {
            output_stream << iteration_number << " " << time_step*iteration_number << std::endl;
        }

        /**
         * @brief Useful possible overload
         *
         */

        void GPSolver::write_operator_splitting_output(size_t, double, std::ostream&){};

        /**
         * @brief Useful possible overload
         *
         */

        void GPSolver::write_operator_splitting_output(size_t, double, double, std::ostream&){};

        /**
         * @brief Useful possible overload
         *
         */

        void GPSolver::write_operator_splitting_output(size_t, double, double, double, std::ostream&){};

    } // namespace GP_Solvers
} // namespace UltraCold


