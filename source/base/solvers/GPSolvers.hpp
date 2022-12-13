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
#include <vector>
#include <random>
#include <chrono>
#include <assert.h>

#include "Vector.hpp"
#include "DFtCalculator.hpp"
#include "DataWriter.hpp"
#include "mesh_fourier_space.hpp"

namespace UltraCold
{

    /**
     * @brief Solver classes for various flavors of Gross-Pitaevskii equations.
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
         *                                                                    output_stream,
         *                                                                    write_output_every);
         *  // or
         *
         *  gp_solver.run_operator_splitting(number_of_time_steps,
         *                                   time_step,
         *                                   output_stream,
         *                                   write_output_every); // for the dynamics
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

                // Re-initializers

                void reinit(Vector<double>&               Vext,
                            Vector<std::complex<double>>& psi_0);
                void reinit(Vector<double>&               Vext,
                            Vector<std::complex<double>>& psi_0,
                            double scattering_length);

                // Destructor

                ~GPSolver()=default;

                // Calculate a ground-state solution

                std::tuple<Vector<std::complex<double>>,double> run_gradient_descent(int    max_num_iter,
                                                                                     double tolerance,
                                                                                     double alpha,
                                                                                     double beta,
                                                                                     std::ostream& output_stream,
                                                                                     int write_output_every);

                // Solve for the dynamics.
                // May be overloaded if necessary (e.g. for time varying scattering length or external potential).
                // Some useful possible overloads are provided, with the function optionally taking a different number
                // of external arguments. Unfortunately, it is currently not possible to declare virtual member
                // functions as templates (so cannot merge a different number of input parameters in the
                // parameter pack of a variadic template)

                virtual void run_operator_splitting(int    number_of_time_steps,
                                                    double time_step,
                                                    std::ostream& output_stream,
                                                    int write_output_every);

                virtual void run_operator_splitting(int,double,double,std::ostream &,int);
                virtual void run_operator_splitting(int,double,double,double,std::ostream &,int);
                virtual void run_operator_splitting(int,double,double,double,double,std::ostream &,int);

                // Set initial conditions for a run in the context of the Truncated Wigner approximation

                void set_tw_initial_conditions(bool system_is_in_harmonic_trap,
                                               bool first_call,
                                               int number_of_modes);

            protected:

                // Write the output.

                virtual void write_gradient_descent_output(size_t iteration_number,
                                                           std::ostream& output_stream);

                virtual void write_operator_splitting_output(size_t iteration_number,
                                                             std::ostream& output_stream);
                virtual void write_operator_splitting_output(size_t, size_t, std::ostream&);
                virtual void write_operator_splitting_output(size_t, double, std::ostream&);
                virtual void write_operator_splitting_output(size_t, double, double, std::ostream&);
                virtual void write_operator_splitting_output(size_t, double, double, double, std::ostream&);

                // Member functions to implement operator splitting, plus some overloads

                virtual void solve_step_1_operator_splitting();
                virtual void solve_step_1_operator_splitting(int);
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
                int write_output_every;

                // Eigenstates of the harmonic oscillator, useful for TWA

                std::vector<Vector<double>> eigenstates_harmonic_oscillator;

        };

        /**
         * @brief Class to solve the Gross-Pitaevskii equation for a dipolar Bose gas in two or three space dimensions.
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
         * with \f$ g=4\pi\hbar^2a/m \f$ the coupling constant fixed by the \f$ s \f$-wave scattering length \f$ a \f$.\n
         * In **three space dimensions**, the dipole-dipole interaction, for dipoles aligned along the same axis,
         * eventually tilted about the z-axis by the angles \f$(\theta,\phi)\f$ (\f$\theta \f$ is the angle between the
         * polarization axis and the z-axis, \f$\phi\f$ is the angle beteen the projection of the magnetic moment into
         * the \f$x\f$-\f$y\f$ plane and the \f$x\f$-axis) is given by the following expression
         * \f[
         * V_{dd}({\bf r})=
         * \frac{\mu_0\mu^2}{4\pi}\frac{1-3\hat{\bf e}\cdot\frac{{\bf r}}{|{\bf r}|}}{|{\bf r}|^3}
         * \f]
         * where \f$ \mu_0 \f$ is the magnetic permeability in vacuum, \f$ \mu \f$ the magnetic dipole moment,
         * \f$ \hat{{\bf e}}=\left(\sin(\theta)\cos(\phi),\sin(\theta)\sin(\phi),\cos(\theta)\right) \f$ the unit vector
         * specifying the polarization direction, and \f${\bf r}\f$ the relative position of two interacting dipoles.
         * The corresponding non-local term in the Hamiltonian is evaluated, in the solver class, by means of Fast
         * Fourier Transforms, exploting the exact knowledge of the dipolar potential in momentum space (see, e.g.,
         * <a href="https://arxiv.org/pdf/0905.0386.pdf"> this </a> review). \n
         *
         * The third term of the Hamiltonian corresponds to the local density approximation of the beyond-mean-field
         * <a href = "https://journals.aps.org/pra/abstract/10.1103/PhysRevA.86.063609"> Lee-Huang-Yang </a> (LHY)
         * correction, with
         *
         * \f[
         * \gamma(\varepsilon_{dd})=\frac{16}{3\sqrt{\pi}} ga^{\frac{3}{2}}\,\mbox{Re}\bigg[\!\int_0^{\pi}
         * \!\!\!\!d\theta\sin\theta [1+\varepsilon_{dd}(3\cos^2\theta-1)]^{\frac{5}{2}}\bigg]\,.
         * \f]
         *
         * where \f$ \varepsilon_{dd}=\mu_0\mu^2/(3g)=a_{dd}/a \f$ gives the ratio between the strength of
         * the dipolar and the contact interaction, eventually written in terms of the dipolar length \f$ a_{dd} \f$
         * and the scattering length \f$ a \f$. \n
         *
         * Experimental measurements and microscopic Monte Carlo calculations have confirmed that the LHY term is an
         * accurate correction to the mean-field theory given by the Gross-Pitaevskii equation in dipolar gases, and it
         * is **fundamental** in the description of the currently available experimental phenomenology. \n
         *
         * In **two space dimensions**, the dipole-dipole interaction assumes a different form. A typical situation of
         * both theoretical and experimental interest is the one in which the system is strongly confined along the
         * \f$z\f$-axis with a tight harmonic trap, such that \f$\hbar\omega_z >> \mu\f$, where \f$ \mu \f$ is the
         * chemical potential. In this condition, the condensate wave-function can be factorized as
         *
         * \f[
         *  \Psi(x,y,z,t) = \chi(z)\psi(x,y,t)
         * \f]
         *
         * where \f$\chi(z)\f$ is the ground-state wave function of the harmonic oscillator along the z-direction.
         * Inserting this ansatz into the mean-field Gross-Pitaevskii equation and integrating out the z-direction,
         * one finds that the contact interaction parameter \f$g\f$ is renormalized by the factor \f$1/\sqrt{2\pi}l_z\f$,
         * where \f$l_z=\sqrt{\frac{\hbar}{m\omega_z}}\f$ is the harmonic oscillator length along the \f$z\f$-axis,
         * while the dipole-dipole interaction potential can be rewritten as (for more details, see for example
         * <a href="https://arxiv.org/pdf/1312.4129.pdf"> here </a>)
         *
         * \f[
         *  \Phi_{dd}({\bf r},t) =
         *  \mathrm{FT}^{-1}\left[\tilde{n}({\bf k},t)\tilde{U}_{dd}^{2D}\left(\frac{{\bf k}l_z}{\sqrt{2}}\right)\right]
         * \f]
         * where \f$\mathrm{FT}^{-1}\f$ indicates the inverse Fourier transform, \f$\tilde{n}({\bf k},t)\f$ is the
         * Fourier transform of the density profile at time \f$t\f$, and the dipole-dipole interaction in
         * this quasi-two-dimensional setup is given, in momentum space, by
         *
         * \f[
         *  U_{dd}^{2D}({\bf q}) = \frac{4\pi}{3}g_{dd}
         *  \left[ F_{\parallel}({\bf q})\sin^2(\alpha) + F_{\perp}({\bf q})\cos^2(\alpha)  \right]
         * \f]
         *
         * with \f$ g_{dd} = \frac{\mu_0\mu^2}{4\pi\sqrt{2\pi}l_z} \f$,
         * \f$ F_{\parallel}({\bf q}) =  -1+3\sqrt{\pi}\frac{q_x^2}{q}\mathrm{e}^{q^2}\mathrm{erfc}(q)\f$,
         * \f$ F_{\perp}({\bf q}) =  2-3\sqrt{\pi}q\mathrm{e}^{q^2}\mathrm{erfc}(q)\f$, where \f$\mathrm{erfc}\f$ is the
         * <a href = "https://en.wikipedia.org/wiki/Error_function"> complementary error function </a>, and \f$\alpha\f$
         * is the angle between the polarization direction and the z-axis (assuming that the projection of the
         * polarization direction in the \f$x\f$-\f$y\f$ plane is along the \f$x\f$-axis).
         *
         * \note
         * In the two-dimensional solver, at the time of writing, the Lee-Huang-Yang correction is **not implemented**.
         * This is because there are no available experiments at this time, and it is not clear whether or not such
         * correction is necessary to describe the physics of the dipolar Bose gas in quasi-2D.
         *
         * This class can be used to solve the extended Gross-Pitaevskii equation both for ground-state configurations
         * as well as for the dynamics, on a cartesian mesh with periodic boundary conditions. The basic usage of the
         * class for **three dimensional problems** is is as follows:
         *
         * \code {.cpp}
         *
         *  DipolarGPSolver dipolar_gp_solver(x,
         *                                    y,
         *                                    z,
         *                                    psi0,
         *                                    Vext,
         *                                    scattering_length,
         *                                    dipolar_length,
         *                                    phi_mu,
         *                                    theta_mu);
         *
         *  // for ground state calculations
         *
         *  std::tie(psi,chemical_potential) = dipolar_gp_solver.run_gradient_descent(maximum_number_of_iterations,
         *                                                                            tolerance,
         *                                                                            alpha,
         *                                                                            beta,
         *                                                                            output_stream,
         *                                                                            write_output_every);
         *
         *  // or
         *
         *  dipolar_gp_solver.run_operator_splitting(number_of_time_steps,time_step); // for the dynamics
         *
         * \endcode
         *
         * while for **two-dimensional problems**, the constructor is the following
         *
         * \code{.cpp}
         *
         *  DipolarGPSolver dipolar_gp_solver(x,
         *                                    y,
         *                                    psi0,
         *                                    Vext,
         *                                    scattering_length,
         *                                    dipolar_length,
         *                                    alpha);
         *
         * \endcode
         *
         * In the constructor, you need to provide the cartesian axis defining the mesh on which the equation
         * is going to be solved, the initial wave function (all the methods are iterative), the external potential,
         * the s-wave scattering length \f$ a \f$ and the dipolar length \f$ a_{dd} \f$, both in appropriate units, as
         * well as the tilting angles. \n
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
                                Vector<std::complex<double>>& psi0,
                                Vector<double>& Vext,
                                double scattering_length,
                                double dipolar_length,
                                double alpha);             // 2D problems

                DipolarGPSolver(Vector<double>& x,
                                Vector<double>& y,
                                Vector<double>& z,
                                Vector<std::complex<double>>& psi0,
                                Vector<double>& Vext,
                                double scattering_length,
                                double dipolar_length,
                                double theta_mu,
                                double phi_mu);             // 3D problems

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
                                                                                     std::ostream& output_stream,
                                                                                     int write_output_every);
                // Solve for the dynamics.
                // May be overloaded if necessary (e.g. for time varying scattering length or external potential).
                // Some useful possible overloads are provided, with the function optionally taking a different number
                // of external arguments. Unfortunately, it is currently not possible to declare virtual member
                // functions as templates (so cannot merge a different number of input parameters in the
                // parameter pack of a variadic template)
                virtual void run_operator_splitting(int    number_of_time_steps,
                                                    double time_step,
                                                    std::ostream& output_stream,
                                                    int write_output_every);
                virtual void run_operator_splitting(int,double,double,std::ostream &,int);
                virtual void run_operator_splitting(int,double,double,double,std::ostream &,int);
                virtual void run_operator_splitting(int,double,double,double,double,std::ostream &,int);

                // Set initial conditions for a run in the context of the Truncated Wigner approximation
                void set_tw_initial_conditions(bool system_is_trapped);

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
                Vector<std::complex<double>> Vtilde;   // This contains the Fourier transform of the dipolar potential
                Vector<std::complex<double>> Phi_dd;
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
                double epsilon_dd;
                double gamma_epsilon_dd;
                double residual;
                double initial_norm;
                double norm;
                double time_step;
                std::complex<double> ci={0.0,1.0};
                int last_iteration_number;
                int write_output_every;
                bool system_is_2d=false;
                bool system_is_3d=false;

        };

    } // namespace GP_Solvers
} // namespace UltraCold

#endif // ULTRACOLD_GP_SOLVERS