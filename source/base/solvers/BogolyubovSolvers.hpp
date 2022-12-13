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


#ifndef ULTRACOLD_BOGOLYUBOV_SOLVERS
#define ULTRACOLD_BOGOLYUBOV_SOLVERS

#include <utility>
#include <tuple>
#include <iostream>
#include <cmath>
#include <omp.h>
#include <complex>
#include <vector>
#include <assert.h>

#include "Vector.hpp"
#include "DFtCalculator.hpp"
#include "DataWriter.hpp"
#include "mesh_fourier_space.hpp"
#include "GPSolvers.hpp"

namespace UltraCold
{

    /**
     * @brief Solver classes for several flavours of Bogolyubov equations.
     *
     * While the GPSolver class and related classes allows to study the ground state and the dynamics of
     * Bose-Einstein condensates, the elementary excitations on top of a certain state can be studied by solving the
     * so-called *Bogolyubov equations*. The idea is to consider a certain configuration, described a condensate
     * wave-function \f$ \psi_0 \f$, and with chemical potential \f$ \mu \f$, and to study small oscillations on top
     * of it by searching for solutions of the Gross-Pitaevskii equation of the form
     *
     * \f[
     *  \psi({\bf r},t) = e^{-i\frac{\mu}{\hbar}t}\left[ \psi_0({\bf r})
     *  + \sum_{n=0}^{\infty} \left( u_n({\bf r})e^{-i\omega_n t} + v^*_n({\bf r})e^{i\omega_n t} \right) \right]
     * \f]
     *
     * Keeping only terms linear the functions \f$ u \f$ and \f$ v \f$, one finds that the *quasi-particle
     * amplitudes* \f$ u_n \f$ and \f$ v_n \f$ and their energies \f$ \hbar\omega_n \f$ are given by the solutions of a
     * certain eigenvalue problem. The form depends on the "flavour" of Gross-Pitaevskii equation one is considering
     * (ordinary, dipolar, mixtures...) and is described in the documentation of the specific solver classes.
     *
     * */

    namespace BogolyubovSolvers
    {

        /**
         * @brief Class to solve the Bogolyubov equations for a trapped Bose gas.
         * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
         *
         * Studying the elementary excitations of a trapped Bose gas in the condensate phase on top of a certain
         * ground (or stationary) state, means searching for solutions of the time-dependent Gross-Pitaevskii
         * equation of the form
         *
         * \f[
         *  \psi({\bf r},t) = e^{-i\frac{\mu}{\hbar}t}\left[ \psi_0({\bf r})
         *  + \sum_{n=0}^{\infty} \left( u_n({\bf r})e^{-i\omega_n t} + v^*_n({\bf r})e^{i\omega_n t} \right) \right]
         * \f]
         *
         * and solving the eigenvalue problem that comes out by keeping only terms linear in the quasi-particle
         * amplitudes \f$ u \f$ and \f$ v \f$. In the case of the ordinary Bose gas (i.e., an ensemble of bosonic
         * particles at very low temperature interacting only via a contact interaction), this amounts to solving the
         * following eigenvalue problem
         *
         * \f[
         *
         *  \begin{bmatrix}
         *      u \\
         *      v
         *  \end{bmatrix}
         *      =
         *  \begin{bmatrix}
         *      -\frac{\hbar^2}{2m}\nabla^2 + V_{ext}({\bf r}) + g |\psi_0|^2 - \mu & g \psi_0^2 \\
         *      - g (\psi_0^*)^2 &  -\left(-\frac{\hbar^2}{2m}\nabla^2 + V_{ext}({\bf r}) + g |\psi_0|^2-\mu\right)
         *  \end{bmatrix}
         *  \begin{bmatrix}
         *      u \\
         *      v
         *  \end{bmatrix}
         *
         * \f]
         *
         * In the case in which the condensate wave function is real (e.g., in absence of vortices, solitons...) the
         * problem can be recast in a more convenient form. In fact, taking the sum and the difference between the
         * two equations, one easily finds
         *
         * \f[
         *  \begin{align}
         *      & \hat{H}\hat{X} (u+v) = (\hbar\omega)^2 (u+v) \nonumber \\
         *      & \hat{X}\hat{H} (u-v) = (\hbar\omega)^2 (u-v) \nonumber \\
         *  \end{align}
         * \f]
         *
         * with
         *
         * \f[
         *  \begin{align}
         *      & \hat{H} = -\frac{\hbar^2}{2m}\nabla^2 + V_{ext}({\bf r}) + g |\psi_0|^2 - \mu\nonumber \\
         *      & \hat{X} = -\frac{\hbar^2}{2m}\nabla^2 + V_{ext}({\bf r}) + 3 g |\psi_0|^2 - \mu \nonumber \\
         *  \end{align}
         * \f]
         *
         * Now, both equations allow to find the (square) of the energy of the Bogolyubov modes, but solving a system
         * of half the dimensionality of the original problem. This typically allows a great saving of computational
         * time. The eigenvectors of the two problems correspond to \f$(u+v)\f$ and \f$(u-v)\f$ respectively, so that
         * if one is interested in finding the Bogolyubov quasi-particle amplitudes \f$ u \f$ and \f$ v \f$, one also
         * needs to solve the second problem, and then set \f$ u = 0.5 \left( (u+v) + (u-v) \right)\f$ and
         * \f$ v = 0.5 \left( (u+v) - (u-v) \right)\f$
         *
         * This class solves the eigenvalue problem using the matrix-free routines provided as part of the package
         * <a href="https://github.com/opencollab/arpack-ng"> arpack-ng </a>, which is distributed as a bundled
         * package with UltraCold.
         *
         * */

        class TrappedBogolyubovSolver
        {

            public:
    
                // Constructors for real problems
    
                TrappedBogolyubovSolver(Vector<double>& x,
                                        Vector<double>& psi0,
                                        Vector<double>& Vext,
                                        double scattering_length,
                                        double chemical_potential,
                                        int number_of_modes,
                                        double tolerance,
                                        int maximum_number_arnoldi_iterations,
                                        bool eigenvectors_required);             // 1D problems
                TrappedBogolyubovSolver(Vector<double>& x,
                                        Vector<double>& y,
                                        Vector<double>& psi0,
                                        Vector<double>& Vext,
                                        double scattering_length,
                                        double chemical_potential,
                                        int number_of_modes,
                                        double tolerance,
                                        int maximum_number_arnoldi_iterations,
                                        bool eigenvectors_required);             // 2D problems
                TrappedBogolyubovSolver(Vector<double>& x,
                                        Vector<double>& y,
                                        Vector<double>& z,
                                        Vector<double>& psi0,
                                        Vector<double>& Vext,
                                        double scattering_length,
                                        double chemical_potential,
                                        int number_of_modes,
                                        double tolerance,
                                        int maximum_number_arnoldi_iterations,
                                        bool eigenvectors_required);             // 3D problems

                // Constructors for complex problems

                TrappedBogolyubovSolver(Vector<double>& x,
                                        Vector<std::complex<double>>& psi0,
                                        Vector<double>& Vext,
                                        double scattering_length,
                                        double chemical_potential,
                                        int number_of_modes,
                                        double tolerance,
                                        int maximum_number_arnoldi_iterations,
                                        bool eigenvectors_required);             // 1D problems
                TrappedBogolyubovSolver(Vector<double>& x,
                                        Vector<double>& y,
                                        Vector<std::complex<double>>& psi0,
                                        Vector<double>& Vext,
                                        double scattering_length,
                                        double chemical_potential,
                                        int number_of_modes,
                                        double tolerance,
                                        int maximum_number_arnoldi_iterations,
                                        bool eigenvectors_required);             // 2D problems
                TrappedBogolyubovSolver(Vector<double>& x,
                                        Vector<double>& y,
                                        Vector<double>& z,
                                        Vector<std::complex<double>>& psi0,
                                        Vector<double>& Vext,
                                        double scattering_length,
                                        double chemical_potential,
                                        int number_of_modes,
                                        double tolerance,
                                        int maximum_number_arnoldi_iterations,
                                        bool eigenvectors_required);             // 3D problems


                // Destructor

                ~TrappedBogolyubovSolver()=default;

                // Runner

                std::tuple< std::vector<std::complex<double>>,
                            std::vector< Vector<std::complex<double>> >,
                            std::vector< Vector<std::complex<double>> > > run();

            protected:

                // Vector data members, needed for the calculations

                Vector<double> Vext;
                Vector<double> x;
                Vector<double> y;
                Vector<double> z;
                Vector<double> kx;
                Vector<double> ky;
                Vector<double> kz;
                Vector<double> kmod2;
                Vector<std::complex<double>> temp;
                Vector<std::complex<double>> temp2;
                Vector<std::complex<double>> temp_tilde;
                Vector<std::complex<double>> temp2_tilde;
                Vector<std::complex<double>> psi0;

                std::vector<Vector<std::complex<double>>> u,v;

                // Number of points along each space dimension

                int nx,ny,nz;

                // Other mesh parameter

                double dx = 1.0;
                double dy = 1.0;
                double dz = 1.0;
                double dv = 1.0;

                // Other useful parameters

                int number_of_modes;
                double chemical_potential;
                double scattering_length;
                bool problem_is_1d=false;
                bool problem_is_2d=false;
                bool problem_is_3d=false;
                bool eigenvectors_required=false;
                bool problem_is_real=false;
                bool problem_is_complex=false;
                double tolerance;
                int maximum_number_of_arnoldi_iterations;

        };

        /**
         * @brief Class to solve the Bogolyubov equations for a trapped **dipolar** Bose gas.
         * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
         *
         * Studying the elementary excitations of a trapped **dipolar** Bose gas in the condensate phase on top of a
         * certain ground (or stationary) state, means searching for solutions of the time-dependent
         * extended Gross-Pitaevskii equation of the form
         *
         * \f[
         *  \psi({\bf r},t) = e^{-i\frac{\mu}{\hbar}t}\left[ \psi_0({\bf r})
         *  + \sum_{n=0}^{\infty} \left( u_n({\bf r})e^{-i\omega_n t} + v^*_n({\bf r})e^{i\omega_n t} \right) \right]
         * \f]
         *
         * and solving the eigenvalue problem that comes out by keeping only terms linear in the quasi-particle
         * amplitudes \f$ u \f$ and \f$ v \f$. In the case of a **dipolar** Bose gas, taking also into account the
         * effects of quantum fluctuations via the
         * <a href = "https://journals.aps.org/pra/abstract/10.1103/PhysRevA.86.063609"> Lee-Huang-Yang </a> (LHY)
         * correction, this amounts to solving the following eigenvalue problem
         *
         * \f[
         *
         *  \begin{bmatrix}
         *      u \\
         *      v
         *  \end{bmatrix}
         *      =
         *  \begin{bmatrix}
         *      \hat{H}-\mu+\hat{X} & \hat{X}^\dagger \\
         *      -\hat{X} & -(\hat{H}-\mu+\hat{X}^\dagger)
         *  \end{bmatrix}
         *  \begin{bmatrix}
         *      u \\
         *      v
         *  \end{bmatrix}
         *
         * \f]
         *
         * with
         *
         * \f[
         *  \begin{align}
         *   \mathcal{H}({\bf r})=-&\frac{\hbar^2}{2m}\nabla^2+V_{\rm ext}({\bf r})+g|\Psi({\bf r},t)|^2+\gamma
         *   (\varepsilon_{dd})|\Psi({\bf r},t)|^3\nonumber\\
         *   +& \int d{\bf r'}V_{dd}({\bf r}-{\bf r'})|\Psi({\bf r'},t)|^2\,,
         *  \end{align}
         * \f]
         *
         * and
         *
         * \f[
         *
         * \hat{X}f({\bf r}) = \psi_0({\bf r})\int d{\bf r}' V_{dd}({\bf r}-{\bf r}')f({\bf r}')\psi_0^*({\bf r}') +
         * \frac{3}{2}\gamma(\varepsilon_{dd})|\psi_0({\bf r})|^3f({\bf r})
         *
         * \f]
         *
         * and finally
         *
         * \f[
         * \gamma(\varepsilon_{dd})=\frac{16}{3\sqrt{\pi}} ga^{\frac{3}{2}}\,\mbox{Re}\bigg[\!\int_0^{\pi}
         * \!\!\!\!d\theta\sin\theta [1+\varepsilon_{dd}(3\cos^2\theta-1)]^{\frac{5}{2}}\bigg]\,.
         * \f]
         *
         * with \f$ g=4\pi\hbar^2a/m \f$ the coupling constant fixed by the \f$ s \f$-wave scattering length \f$ a
         * \f$, \f$ V_{dd}({\bf r}_{i}-{\bf r}_{j})=\frac{\mu_0\mu^2}{4\pi}\frac{1-3\cos^2\theta}{|{\bf r}_{i}-{\bf
         * r}_{j}|^3} \f$ the dipole-dipole potential, being \f$ \mu_0 \f$ the magnetic permeability in vacuum, \f$
         * \mu \f$ the magnetic dipole moment and \f$ \theta \f$ the angle between the vector distance between
         * dipoles and the polarization direction, which we choose as the \f$ x \f$-axis, and
         * \f$ \varepsilon_{dd}=\mu_0\mu^2/(3g)=a_{dd}/a \f$ the ratio between the strength of the dipolar and
         * the contact interaction, eventually written in terms of the dipolar length \f$ a_{dd} \f$ and the
         * scattering length \f$ a \f$.
         *
         * In the case in which the condensate wave function is real (e.g., in absence of vortices, solitons...) the
         * problem can be recast in a more convenient form. In fact, taking the sum and the difference between the
         * two equations, one easily finds
         *
         * \f[
         *  \begin{align}
         *      & (\hat{H}-\mu)(\hat{H}-\mu+2\hat{X}) (u+v) = (\hbar\omega)^2 (u+v) \nonumber \\
         *      & (\hat{H}-\mu+2\hat{X})(\hat{H}-\mu) (u-v) = (\hbar\omega)^2 (u-v) \nonumber \\
         *  \end{align}
         * \f]
         *
         * Now, both equations allow to find the (square) of the energy of the Bogolyubov modes, but solving a system
         * of half the dimensionality of the original problem. This typically allows a great saving of computational
         * time. The eigenvectors of the two problems correspond to \f$(u+v)\f$ and \f$(u-v)\f$ respectively, so that
         * if one is interested in finding the Bogolyubov quasi-particle amplitudes \f$ u \f$ and \f$ v \f$, one also
         * needs to solve the second problem, and then set \f$ u = 0.5 \left( (u+v) + (u-v) \right)\f$ and
         * \f$ v = 0.5 \left( (u+v) - (u-v) \right)\f$
         *
         * This class solves the eigenvalue problem using the matrix-free routines provided as part of the package
         * <a href="https://github.com/opencollab/arpack-ng"> arpack-ng </a>, which is distributed as a bundled
         * package with UltraCold.
         *
         * */

        class TrappedDipolarBogolyubovSolver
        {

            public:

                // Constructors for real problems

//                TrappedDipolarBogolyubovSolver(Vector<double>& x,
//                                               Vector<double>& y,
//                                               Vector<double>& psi0,
//                                               Vector<double>& Vext,
//                                               double scattering_length,
//                                               double dipolar_length,
//                                               double chemical_potential,
//                                               int number_of_modes,
//                                               double tolerance,
//                                               int maximum_number_arnoldi_iterations,
//                                               bool eigenvectors_required); // 2D problems

                TrappedDipolarBogolyubovSolver(Vector<double>& x,
                                               Vector<double>& y,
                                               Vector<double>& z,
                                               Vector<double>& psi0,
                                               Vector<double>& Vext,
                                               double scattering_length,
                                               double dipolar_length,
                                               double theta_mu,
                                               double phi_mu,
                                               double chemical_potential,
                                               int number_of_modes,
                                               double tolerance,
                                               int maximum_number_arnoldi_iterations,
                                               bool eigenvectors_required); // 3D problems

                // Destructor

                ~TrappedDipolarBogolyubovSolver()=default;

                // Runner

                std::tuple<std::vector<std::complex<double>>,
                        std::vector<Vector<std::complex<double>>>,
                        std::vector<Vector<std::complex<double>>>> run();

            protected:

                // Vector data members, needed for the calculations

                Vector<double> Vext;
                Vector<double> x;
                Vector<double> y;
                Vector<double> z;
                Vector<double> kx;
                Vector<double> ky;
                Vector<double> kz;
                Vector<double> kmod2;
                Vector<std::complex<double>> temp;
                Vector<std::complex<double>> temp2;
                Vector<std::complex<double>> temp_tilde;
                Vector<std::complex<double>> temp2_tilde;
                Vector<std::complex<double>> psi0;
                Vector<std::complex<double>> Vtilde;   // This contains the Fourier transform of the dipolar potential
                Vector<std::complex<double>> Phi_dd;   // This will contain the full dipolar potential

                std::vector<Vector<std::complex<double>>> u,v;

                // Number of points along each space dimension

                int nx,ny,nz;

                // Other mesh parameter

                double dx = 1.0;
                double dy = 1.0;
                double dz = 1.0;
                double dv = 1.0;

                // Other useful parameters

                int number_of_modes;
                double chemical_potential;
                double scattering_length;
                double epsilon_dd;
                double gamma_epsilon_dd;
                bool eigenvectors_required=false;
                bool problem_is_real=false;
                bool problem_is_complex=false;

                double tolerance;
                int maximum_number_of_arnoldi_iterations;

        };

    } // namespace BogolyubovSolvers
} // namespace UltraCold

#endif //ULTRACOLD_BOGOLYUBOV_SOLVERS
