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

#ifndef ULTRACOLD_GP_SOLVERS
#define ULTRACOLD_GP_SOLVERS

#include <utility>
#include <tuple>
#include <iostream>
#include <cmath>
#include <omp.h>
#include <complex>

#include "vector.hpp"
#include "dft.hpp"
#include "data_out.hpp"
#include "mesh_fourier_space.hpp"

namespace UltraCold
{

    /**
     * @brief This namespace contains several solver classes for various flavors of Non-Linear Schrodinger equations.
     *
     * Ultra-cold bosonic systems are very often described in terms of the so-called
     * <a href="https://en.wikipedia.org/wiki/Gross%E2%80%93Pitaevskii_equation"> Gross-Pitaevskii equation </a>, which
     * is a non-linear Schrodinger equation arising from a mean-field description of the system, valid when the
     * particles are very weakly interacting and the temperature is close to the absolute zero. In these conditions,
     * an ultra-cold gas of weakly interacting bosonic atoms is described by a complex order parameter \f$ \psi \f$,
     * whose square modulus gives the local density of atoms, and satisfying the following Gross-Pitaevskii equation
     * \f[
     *  i\hbar\frac{\partial \psi}{\partial t} = \left[ \frac{-\hbar^2\nabla^2}{2m}+V_{ext}({\bf r})+g|\psi|^2
     * \right]\psi
     * \f]
     * Here \f$ V_{ext} \f$ is some external trapping potential, and \f$ g \f$ is related to the s-wave scattering
     * length \f$ a \f$ by \f$ g = \frac{4\pi\hbar^2a}{m} \f$. This model is valid as long as the two-body interaction
     * between the atoms can be modeled by a contact interaction of the form
     * \f$ V({\bf r}-{\bf r'})=g\delta({\bf r}-{\bf r'}) \f$. Other models (including, for example, a dipole-dipole
     * interaction and the effects of quantum fluctuations) have a similar form, as it is discussed in other solver
     * classes belonging to this namespace. \n
     * There are two kinds of information we can derive from the solution of the Gross-Pitaevskii equation:
     *  - *Ground-state properties*: these are obtained by searching for stationary solutions of the form
     * \f$ \psi({\bf r},t)=\psi_0({\bf r})e^{-i\mu t/\hbar} \f$, obtaining the stationary Gross-Pitaevskii equation
     *
     * \f[
     *  \mu \psi_0 = \left[ \frac{-\hbar^2\nabla^2}{2m}+V_{ext}({\bf r})+g|\psi_0|^2 \right]\psi_0
     * \f]
     *
     * Solving this for the smallest eigenvalue \f$ \mu \f$, which represents the chemical potential,
     * gives access to the ground-state configuration of the system.
     *  - *Dynamics*: solving the Gross-Pitaevskii equation for appropriate initial conditions allows to simulate the
     * dynamical behavior of the system and to compare the results of the model with experiments.
     *
     * In order to solve these equations, there are several possibilities. \n
     * For what concern the stationary Gross-Pitaevskii equation, one possibility comes from noticing that such equation
     * can be obtained from a constrained minimization formulation of the problem, in particular by requiring that the
     * ground-state order parameter of the system is the exact minimizer of the mean-field energy functional
     *
     * \f[
     *  E[\psi] = \int d{\bf r} \left[ \psi^*({\bf r})\left(\frac{-\hbar^2\nabla^2}{2m}+V_{ext}({\bf r})\right)
     *  \psi({\bf r}) \right]+\frac{g}{2}\int d{\bf r} |\psi({\bf r})|^4
     * \f]
     *
     * under the constraint of a fixed number of particles \f$ \int d{\bf r}|\psi({\bf r})|^2=N \f$. This allows also to
     * introduce the chemical potential \f$ \mu \f$ as the Lagrange multiplier fixing the total number of particles. One
     * can hence find the ground state order parameter and chemical potential using, for example, a
     * <a href="https://en.wikipedia.org/wiki/Gradient_descent"> gradient descent </a> method, as implemented in the
     * class GPSolver (see the related documentation). \n
     * For what concern instead the full Gross-Pitaevskii equation, this is solved using a classical
     * operator-splitting method. See again the full documentation of the GPSolver class for more details.
     *
     */

    namespace GPSolvers
    {

        /**
         * @brief Class to solve the Gross-Pitaevskii equation.
         * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
         *
         * This class allows to solve the Gross-Pitaevskii equation
         *
         * \f[
         *  i\hbar\frac{\partial \psi}{\partial t} = \left[ \frac{-\hbar^2\nabla^2}{2m}+V_{ext}({\bf r})+g|\psi|^2
         * \right]\psi
         * \f]
         *
         * both for ground-state configurations as well as for the dynamics, on a cartesian mesh with periodic
         * boundary conditions. The basic usage of the class is as follows:
         *
         * \code {.cpp}
         *
         *  GPSolver gp_solver(x,y,z,psi0,Vext,scattering_length);
         *
         *  // for ground state calculations
         *
         *  std::tie(psi,chemical_potential) = gp_solver.run_gradient_descent(maximum_number_of_iterations,
         *                                                                    tolerance,
         *                                                                    alpha,
         *                                                                    beta,
         *                                                                    std::cout);
         *  // or
         *
         *  gp_solver.run_operator_splitting(number_of_time_steps,time_step); // for the dynamics
         *
         * \endcode
         *
         * In the constructor, you need to provide the (cartesian) axis along which you want to solve the equation (of
         * course, 1, 2 or 3), the initial wave function (all the methods are iterative) and the external potential. All
         * of these must be properly initialized before passing them to the solver class. \n
         * For ground-state calculations, it is possible to generate some output at each step of the gradient-descent
         * iterations. By default, the GPSolver::run_gradient_descent() just writes, in the standard output, the current
         * iteration number, chemical potential and norm of the residual. This default behavior can however be
         * customized by overriding the GPSolver::write_gradient_descent_output() member function. \n
         * A similar behavior is the default also for dynamics calculations, i.e. for the member function
         * GPSolver::run_operator_splitting(). In this case, the default is just to write the current time step and current
         * time. For examples of usage and customization via overloading, see the example 1 in the
         * <code>examples</code> folder.
         * \note Several member functions of this class take advantage of
         * <a href="https://www.openmp.org/"> OpenMP </a> parallelization. For optimal performance, be sure to run
         * \code
         * $ export OMP_NUM_THREADS=<number of physical cores on the machine used>
         * \endcode
         * on the shell before the program execution.
         * \warning The Gross-Pitaevskii equation is solved in its a-dimensional form. In the case of a harmonic
         * potential in one space dimension it has the form
         *
         * \f[
         *  i\frac{\partial \psi}{\partial t} = \left[ -\frac{1}{2}\frac{\partial^2}{\partial x^2}+\frac{1}{2}x^2
         *  +4\pi a|\psi|^2 \right]\psi
         * \f]
         * where \f$ a \f$ is the scattering length, which here must be given in harmonic units. Be sure to
         * provide the axis, initial wave function, and external potential initialized properly. See the examples in
         * the examples folder.
         *
         */

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

                // Re-initializer

                void reinit(Vector<double>&               Vext,
                            Vector<std::complex<double>>& psi_0);

                // Destructor

                ~GPSolver()=default;

                // Calculate a ground-state solution

                std::tuple<Vector<std::complex<double>>,double> run_gradient_descent(int    max_num_iter,
                                                                                     double tolerance,
                                                                                     double alpha,
                                                                                     double beta,
                                                                                     std::ostream& output_stream);

                // Solve for the dynamics.
                // May be overloaded if necessary (e.g. for time varying scattering length or external potential).
                // Some useful possible overloads are provided, with the function optionally taking a different number
                // of external arguments. Unfortunately, it is currently not possible to declare virtual member
                // functions as templates (so cannot merge a different number of input parameters in the
                // parameter pack of a variadic template)

                virtual void run_operator_splitting(int    number_of_time_steps,
                                                    double time_step,
                                                    std::ostream& output_stream);

                virtual void run_operator_splitting(int,double,double,std::ostream &);
                virtual void run_operator_splitting(int,double,double,double,std::ostream &);
                virtual void run_operator_splitting(int,double,double,double,double,std::ostream &);

            protected:

                // Write the output.

                virtual void write_gradient_descent_output(size_t iteration_number,
                                                           std::ostream& output_stream);

                virtual void write_operator_splitting_output(size_t iteration_number,
                                                             std::ostream& output_stream);
                virtual void write_operator_splitting_output(size_t, double, std::ostream&);
                virtual void write_operator_splitting_output(size_t, double, double, std::ostream&);
                virtual void write_operator_splitting_output(size_t, double, double, double, std::ostream&);

                // Member functions to implement operator splitting, plus some overloads

                virtual void solve_step_1_operator_splitting();
                virtual void solve_step_1_operator_splitting(double);
                virtual void solve_step_1_operator_splitting(double, double);
                virtual void solve_step_1_operator_splitting(double, double, double);
                void solve_step_2_operator_splitting(MKLWrappers::DFtCalculator&);

                // Vector data members, needed for the calculations

                Vector<std::complex<double>> psi;
                Vector<double>    Vext;
                Vector<double>    x;
                Vector<double>    y;
                Vector<double>    z;
                Vector<double>    kx;
                Vector<double>    ky;
                Vector<double>    kz;
                Vector<double>    kmod2;               // This contains the squared moduli of the k-vectors
                Vector<std::complex<double>> psitilde; // This contains the Fourier transform of psi
                Vector<std::complex<double>> hpsi;     // This contains the result of \hat{H}\psi

                // Number of points along each space dimension

                int nx,ny,nz;

                // Other mesh parameter

                double dx = 1.0;
                double dy = 1.0;
                double dz = 1.0;
                double dv = 1.0;

                // Other useful parameters

                double chemical_potential;
                double scattering_length;
                double residual;
                double initial_norm;
                double norm;
                double time_step;
                std::complex<double> ci={0.0,1.0};
                int last_iteration_number;

        };

        /**
         * @brief Class to solve the Gross-Pitaevskii equation for a dipolar Bose gas with the Lee-Huang-Yang
         * correction.
         * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
         *
         * This class allows to solve the extended Gross-Pitaevskii equation for a dipolar Bose gas
         *
         * \f[
         * \begin{align}
         *  &i\hbar\frac{\partial}{\partial t}  \Psi({\bf r},t)= \mathcal{H}({\bf r})\,\Psi({\bf r},t)\,,
         * \end{align}
         * \f]
         *
         * where the Hamiltonian \f$ H \f$ is
         *
         * \f[
         * \begin{align}
         *  \mathcal{H}({\bf r})=-&\frac{\hbar^2}{2m}\nabla^2+V_{\rm ext}({\bf r})+g|\Psi({\bf r},t)|^2+\gamma
         *  (\varepsilon_{dd})|\Psi({\bf r},t)|^3\nonumber\\
         *  +& \int d{\bf r'}V_{dd}({\bf r}-{\bf r'})|\Psi({\bf r'},t)|^2\,,
         * \end{align}
         * \f]
         *
         * with \f$ g=4\pi\hbar^2a/m \f$ the coupling constant fixed by the \f$ s \f$-wave scattering length \f$ a \f$
         * and \f$ V_{dd}({\bf r}_{i}-{\bf r}_{j})=\frac{\mu_0\mu^2}{4\pi}\frac{1-3\cos^2\theta}{|{\bf r}_{i}-{\bf
         * r}_{j}|^3} \f$ the dipole-dipole potential, being \f$ \mu_0 \f$ the magnetic permeability in vacuum, \f$
         * \mu \f$ the magnetic dipole moment and \f$ \theta \f$ the angle between the vector distance between
         * dipoles and the polarization direction, which we choose as the \f$ x \f$-axis. In the absence of trapping,
         * the system can be fully characterised by the single parameter
         * \f$ \varepsilon_{dd}=\mu_0\mu^2/(3g)=a_{dd}/a \f$, i.e., the ratio between the strength of the dipolar and
         * the contact interaction, eventually written in terms of the dipolar length \f$ a_{dd} \f$ and the
         * scattering length \f$ a \f$. The third term of the Hamiltonian corresponds to the local density
         * approximation of the beyond-mean-field
         * <a href = "https://journals.aps.org/pra/abstract/10.1103/PhysRevA.86.063609"> Lee-Huang-Yang </a> (LHY)
         * correction with
         *
         * \f[
         * \gamma(\varepsilon_{dd})=\frac{16}{3\sqrt{\pi}} ga^{\frac{3}{2}}\,\mbox{Re}\bigg[\!\int_0^{\pi}
         * \!\!\!\!d\theta\sin\theta [1+\varepsilon_{dd}(3\cos^2\theta-1)]^{\frac{5}{2}}\bigg]\,.
         * \f]
         *
         * Experimental measurements and microscopic Monte Carlo calculations have confirmed that the LHY term is an
         * accurate correction to the mean-field theory given by the Gross-Pitaevskii equation in dipolar gases. \n
         * This class can be used to solve the extended Gross-Pitaevskii equation both for ground-state configurations
         * as well as for the dynamics, on a cartesian mesh with periodic boundary conditions. The basic usage of the
         * class is as follows:
         *
         * \code {.cpp}
         *
         *  DipolarGPSolver dipolar_gp_solver(x,
         *                                    y,
         *                                    z,
         *                                    psi0,
         *                                    Vext,
         *                                    scattering_length,
         *                                    dipolar_length);
         *
         *  // for ground state calculations
         *
         *  std::tie(psi,chemical_potential) = dipolar_gp_solver.run_gradient_descent(maximum_number_of_iterations,
         *                                                                            tolerance,
         *                                                                            alpha,
         *                                                                            beta,
         *                                                                            std::cout);
         *
         *  // or
         *
         *  dipolar_gp_solver.run_operator_splitting(number_of_time_steps,time_step); // for the dynamics
         *
         * \endcode
         *
         * In the constructor, you need to provide the three cartesian axis defining the mesh on which the equation
         * is going to be solved, the initial wave function (all the methods are iterative), the external potential,
         * the s-wave scattering length \f$ a \f$ and the dipolar length \f$ a_{dd} \f$, both in appropriate units. \n
         * The algorithms used and their working principles are the same as those used for a non-dipolar Bose gas,
         * described in the documentation of the GPSolver class. See that documentation for details.
         * \note Several member functions of this class take advantage of
         * <a href="https://www.openmp.org/"> OpenMP </a> parallelization. For optimal performance, be sure to run
         * \code
         * $ export OMP_NUM_THREADS=<number of physical cores on the machine used>
         * \endcode
         * on the shell before the program execution.
         * \warning The extended Gross-Pitaevskii equation is solved in its a-dimensional form, which, in the case
         * of, for example, harmonic trapping, is given by
         *
         * \f[
         *  i\frac{\partial \psi({\bf r},t)}{\partial t} = \left[ -\frac{\nabla^2}{2}
         *  + \left( \frac{1}{2} \left( \frac{x}{a_x^2} \right)^2 + \frac{1}{2} \left( \frac{y}{a_y^2} \right)^2 +
         *  \frac{1}{2} \left( \frac{z}{a_z^2} \right)^2  \right)
         *  + 4\pi a |\psi({\bf r},t)|^2
         *  + 3 a_{dd} \int d{\bf r'}\frac{1-3\cos^2\theta}{|{\bf r}-{\bf r}'|^3}|\psi({\bf r'},t)|^2
         *  + \frac{64\sqrt{\pi}}{3}a^{5/2} \mbox{Re}\bigg[\!\int_0^{\pi}
         * \!\!\!\!d\theta\sin\theta [1+\varepsilon_{dd}(3\cos^2\theta-1)]^{\frac{5}{2}}\bigg]\, |\psi({\bf r},t)|^3
         *  \right]\psi({\bf r},t)
         * \f]
         *
         * where \f$ a_{x,y,z}=\sqrt{\frac{\hbar}{m\omega_{x,y,z}}} \f$ and all the lengths (including the s-wave and
         * the dipolar length) are in units of the harmonic oscillator length \f$ a_{ho} = (a_xa_ya_z)^{1/3} \f$
         *
         */

        class DipolarGPSolver
        {

            public:

                // Constructors

                DipolarGPSolver(Vector<double>& x,
                                Vector<double>& y,
                                Vector<double>& z,
                                Vector<std::complex<double>>& psi0,
                                Vector<double>& Vext,
                                double scattering_length,
                                double dipolar_length);             // 3D problems

                // Re-initializer

                void reinit(Vector<double>&               Vext,
                            Vector<std::complex<double>>& psi);

                // Destructor

                ~DipolarGPSolver()=default;

                // GP_Solvers

                // Calculate a ground-state solution

                std::tuple<Vector<std::complex<double>>,double> run_gradient_descent(int max_num_iter,
                                                                                     double tolerance,
                                                                                     double alpha,
                                                                                     double beta,
                                                                                     std::ostream& output_stream);
                // Solve for the dynamics.
                // May be overloaded if necessary (e.g. for time varying scattering length or external potential).
                // Some useful possible overloads are provided, with the function optionally taking a different number
                // of external arguments. Unfortunately, it is currently not possible to declare virtual member
                // functions as templates (so cannot merge a different number of input parameters in the
                // parameter pack of a variadic template)

                virtual void run_operator_splitting(int    number_of_time_steps,
                                                    double time_step,
                                                    std::ostream& output_stream);

                virtual void run_operator_splitting(int,double,double,std::ostream &);
                virtual void run_operator_splitting(int,double,double,double,std::ostream &);
                virtual void run_operator_splitting(int,double,double,double,double,std::ostream &);

            protected:

                // Write the output.

                virtual void write_gradient_descent_output(size_t iteration_number,
                                                           std::ostream& output_stream);

                virtual void write_operator_splitting_output(size_t        iteration_number,
                                                             std::ostream& output_stream);
                virtual void write_operator_splitting_output(size_t, double, std::ostream&);
                virtual void write_operator_splitting_output(size_t, double, double, std::ostream&);
                virtual void write_operator_splitting_output(size_t, double, double, double, std::ostream&);

                // Member functions to implement operator splitting, plus some overloads

                virtual void solve_step_1_operator_splitting(MKLWrappers::DFtCalculator&);
                virtual void solve_step_1_operator_splitting(MKLWrappers::DFtCalculator&, double);
                virtual void solve_step_1_operator_splitting(MKLWrappers::DFtCalculator&,double, double);
                virtual void solve_step_1_operator_splitting(MKLWrappers::DFtCalculator&,double, double, double);
                void solve_step_2_operator_splitting(MKLWrappers::DFtCalculator&);

                // Vector data members, needed for the calculations

                Vector<std::complex<double>> psi;
                Vector<double>    Vext;
                Vector<double>    x;
                Vector<double>    y;
                Vector<double>    z;
                Vector<double>    kx;
                Vector<double>    ky;
                Vector<double>    kz;
                Vector<double>    kmod2;               // This contains the squared moduli of the k-vectors
                Vector<std::complex<double>> psitilde; // This contains the Fourier transform of psi
                Vector<std::complex<double>> hpsi;     // This contains the result of \hat{H}\psi
                Vector<double> Vtilde;   // This contains the Fourier transform of the dipolar potential
                Vector<double> Phi_dd;
                Vector<std::complex<double>> Phi_tilde;

                // Number of points along each space dimension

                int nx,ny,nz;

                // Other mesh parameter

                double dx = 1.0;
                double dy = 1.0;
                double dz = 1.0;
                double dv = 1.0;

                // Other useful parameters

                double chemical_potential;
                double scattering_length;
                double dipolar_length;
                double epsilon_dd;
                double gamma_epsilon_dd;
                double residual;
                double initial_norm;
                double norm;
                double time_step;
                std::complex<double> ci={0.0,1.0};
                int last_iteration_number;

        };
    } // namespace GP_Solvers
} // namespace UltraCold

#endif // ULTRACOLD_GP_SOLVERS