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

#include "BogolyubovSolvers.hpp"
#include "arpack.hpp"

#define PI 3.1415926535897932384626433
#define TWOPI (2*PI)

namespace UltraCold
{
    namespace BogolyubovSolvers
    {

        // Constructors for real problems

        /**
         * @brief Constructor for a TrappedBogolyubovSolver in one space dimension and real stationary condensate
         * wave function
         *
         * @param x *Vector<double>* representing the cartesian axis on which the Bogolyubov
         * equations in one space dimension will be solved
         * @param psi_0 *Vector<double>* representing a stationary condensate wave function
         * @param Vext *Vector<double>* representing the external potential.
         * @param scattering_length *double* the scattering length in appropriate units
         * @param chemical_potential *double* the ground-state chemical potential
         * @param eigenvectors_required *bool* specifies if also the eigenvectors are required. Default is false.
         *
         */

        TrappedBogolyubovSolver::TrappedBogolyubovSolver(Vector<double>& x,
                                                         Vector<double>& psi0,
                                                         Vector<double>& Vext,
                                                         double scattering_length,
                                                         double chemical_potential,
                                                         int number_of_modes,
                                                         double tolerance,
                                                         int maximum_number_arnoldi_iterations,
                                                         bool eigenvectors_required)
        {

            // Get necessary copies of input data
            this->x    = x;
            this->Vext = Vext;
            this->scattering_length = scattering_length;
            this->number_of_modes = number_of_modes;
            this->chemical_potential = chemical_potential;
            this->tolerance = tolerance;
            this->maximum_number_of_arnoldi_iterations = maximum_number_arnoldi_iterations;
            this->eigenvectors_required=eigenvectors_required;

            // Check the shape and dimensions of the Vectors provided are consistent
            nx = this->x.extent(0);
            ny = 0;
            nz = 0;
            assert(x.order()==1);
            assert(psi0.order() == 1);
            assert(Vext.order() == 1);
            assert(psi0.extent(0) == nx);
            assert(Vext.extent(0) == nx);

            // Initialize the mesh in Fourier space
            kx.reinit(nx);
            kmod2.reinit(nx);
            create_mesh_in_Fourier_space(this->x,kx);
            for (size_t i = 0; i < nx; ++i)
                kmod2(i) = std::pow(kx(i),2);

            // Initialize additional workspace vectors
            temp.reinit(nx);
            temp2.reinit(nx);
            temp_tilde.reinit(nx);
            temp2_tilde.reinit(nx);

            // Initialize space steps
            dx = this->x(1)-this->x(0);
            dv = dx;

            // Initialize boolean
            problem_is_1d = true;
            problem_is_real = true;

            // Eventually resize u and v
            if (this->eigenvectors_required)
            {
                u.resize(this->number_of_modes);
                v.resize(this->number_of_modes);
                for (int i = 0; i < this->number_of_modes; ++i)
                {
                    u[i].reinit(nx);
                    v[i].reinit(nx);
                }
            }

            // Reinitialize psi0
            this->psi0.reinit(nx);
            for (int i = 0; i < nx; ++i) this->psi0[i] = psi0[i];

        }

        /**
         * @brief Constructor for a TrappedBogolyubovSolver in two space dimensions and real stationary condensate
         * wave function
         *
         * @param x *Vector<double>* representing the x-axis on which the Bogolyubov
         * equations in one space dimension will be solved
         * @param y *Vector<double>* representing the y-axis on which the Bogolyubov
         * equations in one space dimension will be solved
         * @param psi_0 *Vector<double>* representing a stationary condensate wave function
         * @param Vext *Vector<double>* representing the external potential.
         * @param scattering_length *double* the scattering length in appropriate units
         * @param chemical_potential *double* the ground-state chemical potential
         * @param eigenvectors_required *bool* specifies if also the eigenvectors are required. Default is false.
         *
         */

        TrappedBogolyubovSolver::TrappedBogolyubovSolver(Vector<double>& x,
                                                         Vector<double>& y,
                                                         Vector<double>& psi0,
                                                         Vector<double>& Vext,
                                                         double scattering_length,
                                                         double chemical_potential,
                                                         int number_of_modes,
                                                         double tolerance,
                                                         int maximum_number_arnoldi_iterations,
                                                         bool eigenvectors_required)
        {

            // Get necessary copies of input data
            this->x    = x;
            this->y    = y;
            this->Vext = Vext;
            this->scattering_length = scattering_length;
            this->number_of_modes = number_of_modes;
            this->chemical_potential = chemical_potential;
            this->tolerance = tolerance;
            this->maximum_number_of_arnoldi_iterations = maximum_number_arnoldi_iterations;
            this->eigenvectors_required=eigenvectors_required;

            // Check that the shape and dimensions of the Vectors provided are consistent
            nx = this->x.extent(0);
            ny = this->y.extent(0);
            nz = 0;
            assert(x.order()==1);
            assert(y.order()==1);
            assert(psi0.order() == 2);
            assert(Vext.order() == 2);
            assert(psi0.extent(0) == nx);
            assert(psi0.extent(1) == ny);

            // Initialize the mesh in Fourier space
            kx.reinit(nx);
            ky.reinit(ny);
            kmod2.reinit(nx,ny);
            create_mesh_in_Fourier_space(this->x,this->y,kx,ky);

            for (size_t i = 0; i < nx; ++i)
                for (size_t j = 0; j < ny; ++j)
                    kmod2(i,j) = std::pow(kx(i),2) +
                                 std::pow(ky(j),2);

            // Initialize additional workspace vectors
            temp.reinit(nx,ny);
            temp2.reinit(nx,ny);
            temp_tilde.reinit(nx,ny);
            temp2_tilde.reinit(nx,ny);

            // Initialize space steps
            dx = this->x(1)-this->x(0);
            dy = this->y(1)-this->y(0);
            dv = dx*dy;

            // Initialize boolean
            problem_is_2d = true;
            problem_is_real = true;

            // Eventually resize u and v
            if (this->eigenvectors_required)
            {
                u.resize(this->number_of_modes);
                v.resize(this->number_of_modes);
                for (int i = 0; i < this->number_of_modes; ++i)
                {
                    u[i].reinit(nx,ny);
                    v[i].reinit(nx,ny);
                }
            }

            // Reinitialize psi0
            this->psi0.reinit(nx,ny);
            for (int i = 0; i < nx*ny; ++i) this->psi0[i] = psi0[i];

        }

        /**
         * @brief Constructor for a TrappedBogolyubovSolver in three space dimensions and real stationary condensate
         * wave function
         *
         * @param x *Vector<double>* representing the x-axis on which the Bogolyubov
         * equations in one space dimension will be solved
         * @param y *Vector<double>* representing the y-axis on which the Bogolyubov
         * equations in one space dimension will be solved
         * @param z *Vector<double>* representing the z-axis on which the Bogolyubov
         * equations in one space dimension will be solved
         * @param psi_0 *Vector<double>* representing a stationary condensate wave function
         * @param Vext *Vector<double>* representing the external potential.
         * @param scattering_length *double* the scattering length in appropriate units
         * @param chemical_potential *double* the ground-state chemical potential
         * @param eigenvectors_required *bool* specifies if also the eigenvectors are required. Default is false.
         *
         */

        TrappedBogolyubovSolver::TrappedBogolyubovSolver(Vector<double>& x,
                                                         Vector<double>& y,
                                                         Vector<double>& z,
                                                         Vector<double>& psi0,
                                                         Vector<double>& Vext,
                                                         double scattering_length,
                                                         double chemical_potential,
                                                         int number_of_modes,
                                                         double tolerance,
                                                         int maximum_number_arnoldi_iterations,
                                                         bool eigenvectors_required)
        {

            // Get necessary copies of input data
            this->x    = x;
            this->y    = y;
            this->z    = z;
            this->Vext = Vext;
            this->scattering_length = scattering_length;
            this->number_of_modes = number_of_modes;
            this->chemical_potential = chemical_potential;
            this->tolerance = tolerance;
            this->maximum_number_of_arnoldi_iterations = maximum_number_arnoldi_iterations;
            this->eigenvectors_required=eigenvectors_required;

            // Check the dimensions of the Vectors provided are consistent
            nx = this->x.extent(0);
            ny = this->y.extent(0);
            nz = this->z.extent(0);
            assert(x.order()==1);
            assert(y.order()==1);
            assert(z.order()==1);
            assert(psi0.order() == 3);
            assert(Vext.order() == 3);
            assert(psi0.extent(0) == nx);
            assert(psi0.extent(1) == ny);
            assert(psi0.extent(2) == nz);
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
                        kmod2(i,j,k) = std::pow(kx(i),2) +
                                       std::pow(ky(j),2) +
                                       std::pow(kz[k],2);

            // Initialize hpsi and psitilde
            temp.reinit(nx,ny,nz);
            temp2.reinit(nx,ny,nz);
            temp_tilde.reinit(nx,ny,nz);
            temp2_tilde.reinit(nx,ny,nz);

            // Initialize space steps
            dx = this->x(1)-this->x(0);
            dy = this->y(1)-this->y(0);
            dy = this->z(1)-this->z(0);
            dv = dx*dy*dz;

            // Initialize boolean
            problem_is_3d = true;
            problem_is_real = true;

            // Eventually resize u and v
            if (this->eigenvectors_required)
            {
                u.resize(this->number_of_modes);
                v.resize(this->number_of_modes);
                for (int i = 0; i < this->number_of_modes; ++i)
                {
                    u[i].reinit(nx,ny,nz);
                    v[i].reinit(nx,ny,nz);
                }
            }

            // Reinitialize psi0
            this->psi0.reinit(nx,ny,nz);
            for (int i = 0; i < nx*ny*nz; ++i) this->psi0[i] = psi0[i];

        }

        /**
         * @brief Constructor for a TrappedBogolyubovSolver in one space dimension and complex stationary condensate
         * wave function
         *
         * @param x *Vector<double>* representing the cartesian axis on which the Bogolyubov
         * equations in one space dimension will be solved
         * @param psi_0 *Vector<std::complex<double>>* representing a stationary condensate wave function
         * @param Vext *Vector<double>* representing the external potential.
         * @param scattering_length *double* the scattering length in appropriate units
         * @param chemical_potential *double* the ground-state chemical potential
         * @param eigenvectors_required *bool* specifies if also the eigenvectors are required. Default is false.
         *
         */

        TrappedBogolyubovSolver::TrappedBogolyubovSolver(Vector<double>& x,
                                                         Vector<std::complex<double>>& psi0,
                                                         Vector<double>& Vext,
                                                         double scattering_length,
                                                         double chemical_potential,
                                                         int number_of_modes,
                                                         double tolerance,
                                                         int maximum_number_arnoldi_iterations,
                                                         bool eigenvectors_required)
        {

            // Get necessary copies of input data
            this->x    = x;
            this->Vext = Vext;
            this->scattering_length = scattering_length;
            this->psi0 = psi0;
            this->number_of_modes = number_of_modes;
            this->chemical_potential = chemical_potential;
            this->tolerance = tolerance;
            this->maximum_number_of_arnoldi_iterations = maximum_number_arnoldi_iterations;
            this->eigenvectors_required=eigenvectors_required;

            // Check the dimensions of the Vectors provided are consistent
            nx = this->x.extent(0);
            ny = 0;
            nz = 0;
            assert(x.order()==1);
            assert(psi0.order() == 1);
            assert(Vext.order() == 1);
            assert(psi0.extent(0) == nx);
            assert(Vext.extent(0) == nx);

            // Initialize the mesh in Fourier space
            kx.reinit(nx);
            kmod2.reinit(nx);
            create_mesh_in_Fourier_space(this->x,kx);
            for (size_t i = 0; i < nx; ++i)
                kmod2(i) = std::pow(kx(i),2);

            // Initialize additional workspace vectors
            temp.reinit(nx);
            temp2.reinit(nx);
            temp_tilde.reinit(nx);
            temp2_tilde.reinit(nx);

            // Initialize space steps
            dx = this->x(1)-this->x(0);
            dv = dx;

            // Initialize boolean
            problem_is_1d = true;
            problem_is_complex = true;

            // Eventually resize u and v
            if (this->eigenvectors_required)
            {
                u.resize(this->number_of_modes);
                v.resize(this->number_of_modes);
                for (int i = 0; i < this->number_of_modes; ++i)
                {
                    u[i].reinit(nx);
                    v[i].reinit(nx);
                }
            }

        }

        /**
         * @brief Constructor for a TrappedBogolyubovSolver in two space dimensions and complex stationary condensate
         * wave function
         *
         * @param x *Vector<double>* representing the x-axis on which the Bogolyubov
         * equations in one space dimension will be solved
         * @param y *Vector<double>* representing the y-axis on which the Bogolyubov
         * equations in one space dimension will be solved
         * @param psi_0 *Vector<std::complex<double>>* representing a stationary condensate wave function
         * @param Vext *Vector<double>* representing the external potential.
         * @param scattering_length *double* the scattering length in appropriate units
         * @param chemical_potential *double* the ground-state chemical potential
         * @param eigenvectors_required *bool* specifies if also the eigenvectors are required. Default is false.
         *
         */

        TrappedBogolyubovSolver::TrappedBogolyubovSolver(Vector<double>& x,
                                                         Vector<double>& y,
                                                         Vector<std::complex<double>>& psi0,
                                                         Vector<double>& Vext,
                                                         double scattering_length,
                                                         double chemical_potential,
                                                         int number_of_modes,
                                                         double tolerance,
                                                         int maximum_number_arnoldi_iterations,
                                                         bool eigenvectors_required)
        {

            // Get necessary copies of input data
            this->x    = x;
            this->y    = y;
            this->Vext = Vext;
            this->scattering_length = scattering_length;
            this->psi0 = psi0;
            this->number_of_modes = number_of_modes;
            this->chemical_potential = chemical_potential;
            this->tolerance = tolerance;
            this->maximum_number_of_arnoldi_iterations = maximum_number_arnoldi_iterations;
            this->eigenvectors_required=eigenvectors_required;

            // Check the dimensions of the Vectors provided are consistent
            nx = this->x.extent(0);
            ny = this->y.extent(0);
            nz = 0;
            assert(x.order()==1);
            assert(y.order()==1);
            assert(psi0.order() == 2);
            assert(Vext.order() == 2);
            assert(psi0.extent(0) == nx);
            assert(psi0.extent(1) == ny);
            assert(Vext.extent(0) == nx);
            assert(Vext.extent(1) == ny);

            // Initialize the mesh in Fourier space
            kx.reinit(nx);
            ky.reinit(ny);
            kmod2.reinit(nx,ny);
            create_mesh_in_Fourier_space(this->x,this->y,kx,ky);
            for (size_t i = 0; i < nx; ++i)
                for (size_t j = 0; j < ny; ++j)
                    kmod2(i,j) = std::pow(kx(i),2) +
                                 std::pow(ky(j),2);

            // Initialize additional workspace vectors
            temp.reinit(nx,ny);
            temp2.reinit(nx,ny);
            temp_tilde.reinit(nx,ny);
            temp2_tilde.reinit(nx,ny);

            // Initialize space steps
            dx = this->x(1)-this->x(0);
            dy = this->y(1)-this->y(0);
            dv = dx*dy;

            // Initialize boolean
            problem_is_2d = true;
            problem_is_complex = true;

            // Eventually resize u and v
            if (this->eigenvectors_required)
            {
                u.resize(this->number_of_modes);
                v.resize(this->number_of_modes);
                for (int i = 0; i < this->number_of_modes; ++i)
                {
                    u[i].reinit(nx,ny);
                    v[i].reinit(nx,ny);
                }
            }

        }

        /**
         * @brief Constructor for a TrappedBogolyubovSolver in three space dimensions and complex stationary condensate
         * wave function
         *
         * @param x *Vector<double>* representing the x-axis on which the Bogolyubov
         * equations in one space dimension will be solved
         * @param y *Vector<double>* representing the y-axis on which the Bogolyubov
         * equations in one space dimension will be solved
         * @param z *Vector<double>* representing the z-axis on which the Bogolyubov
         * equations in one space dimension will be solved
         * @param psi_0 *Vector<std::complex<double>>* representing a stationary condensate wave function
         * @param Vext *Vector<double>* representing the external potential.
         * @param scattering_length *double* the scattering length in appropriate units
         * @param chemical_potential *double* the ground-state chemical potential
         * @param eigenvectors_required *bool* specifies if also the eigenvectors are required. Default is false.
         *
         */

        TrappedBogolyubovSolver::TrappedBogolyubovSolver(Vector<double>& x,
                                                         Vector<double>& y,
                                                         Vector<double>& z,
                                                         Vector<std::complex<double>>& psi0,
                                                         Vector<double>& Vext,
                                                         double scattering_length,
                                                         double chemical_potential,
                                                         int number_of_modes,
                                                         double tolerance,
                                                         int maximum_number_arnoldi_iterations,
                                                         bool eigenvectors_required)
        {

            // Get necessary copies of input data
            this->x    = x;
            this->y    = y;
            this->z    = z;
            this->Vext = Vext;
            this->scattering_length = scattering_length;
            this->psi0 = psi0;
            this->number_of_modes = number_of_modes;
            this->chemical_potential = chemical_potential;
            this->tolerance = tolerance;
            this->maximum_number_of_arnoldi_iterations = maximum_number_arnoldi_iterations;
            this->eigenvectors_required=eigenvectors_required;

            // Check the dimensions of the Vectors provided are consistent
            nx = this->x.extent(0);
            ny = this->y.extent(0);
            nz = this->z.extent(0);
            assert(x.order()==1);
            assert(y.order()==1);
            assert(z.order()==1);
            assert(psi0.order() == 3);
            assert(Vext.order() == 3);
            assert(psi0.extent(0) == nx);
            assert(psi0.extent(2) == ny);
            assert(psi0.extent(2) == nz);
            assert(Vext.extent(0) == nx);
            assert(Vext.extent(2) == ny);
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
                        kmod2(i,j,k) = std::pow(kx(i),2) +
                                       std::pow(ky(j),2) +
                                       std::pow(kz[k],2);

            // Initialize hpsi and psitilde
            temp.reinit(nx,ny,nz);
            temp2.reinit(nx,ny,nz);
            temp_tilde.reinit(nx,ny,nz);
            temp2_tilde.reinit(nx,ny,nz);

            // Initialize space steps
            dx = this->x(1)-this->x(0);
            dy = this->y(1)-this->y(0);
            dy = this->z(1)-this->z(0);
            dv = dx*dy*dz;

            // Initialize boolean
            problem_is_3d = true;
            problem_is_complex = true;

            // Eventually resize u and v
            if (this->eigenvectors_required)
            {
                u.resize(this->number_of_modes);
                v.resize(this->number_of_modes);
                for (int i = 0; i < this->number_of_modes; ++i)
                {
                    u[i].reinit(nx,ny,nz);
                    v[i].reinit(nx,ny,nz);
                }
            }
        }


        /**
         * @brief Solve the Bogolyubov equations.
         *
         * This member function actually solves the Bogolyubov equations for a trapped condensate. In the case in which
         * a real wave-function was passed to the constructor, the function will use the following useful recast of the
         * problem in one of halved dimensionality. Starting from
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
         *      - g (\psi_0)^2 &  -\left(-\frac{\hbar^2}{2m}\nabla^2 + V_{ext}({\bf r}) + g |\psi_0|^2-\mu\right)
         *  \end{bmatrix}
         *  \begin{bmatrix}
         *      u \\
         *      v
         *  \end{bmatrix}
         *
         * \f]
         *
         * and taking the sum and the difference between two equations, one obtain the following equivalent formulation
         * of the problem
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
         *
         *  If the eigenvectors are not required, this function will solve only the first equation, calculating the
         *  square of the energies of the eigen-modes, but *returning the energies* (i.e., it calculates the square
         *  root before returning). If the eigenvectors are also requested, the function will solve also the second
         *  eigenvalue problem, obtaining the Bogolyubov quasi-particle amplitudes \f$ u \f$ and \f$ v \f$ by taking
         *  (half) the sum and the difference between the eigenvectors of the two problems.
         *
         *  In the case in which, instead, the wave-function passed to the constructor is complex, the function will
         *  solve, by brute force, the complete problem.
         *
         * */

        std::tuple<std::vector<std::complex<double>>,
                std::vector<Vector<std::complex<double>>>,
                std::vector<Vector<std::complex<double>>>> TrappedBogolyubovSolver::run()
        {

            // Initialize the variables
            a_int N;
            if(problem_is_1d)      N=nx;
            else if(problem_is_2d) N=nx*ny;
            else if(problem_is_3d) N=nx*ny*nz;

            a_int nev = number_of_modes;
            a_int ncv = 2 * nev + 1;
            a_int ldv;
            if(problem_is_real)
                ldv = N;
            else if (problem_is_complex)
                ldv=2*N;
            a_int ldz = ldv + 1;
            a_int lworkl = 3 * (ncv * ncv) + 6 * ncv;
            a_int rvec = 0;
            if(eigenvectors_required)
                rvec = 1;

            std::array<a_int, 11> iparam{};
            iparam[0] = 1;
            iparam[2] = maximum_number_of_arnoldi_iterations;
            iparam[3] = 1;
            iparam[6] = 1;
            std::array<a_int, 14> ipntr{};

            a_int info = 0, ido = 0;

            std::vector<std::complex<double>> eigenvalues(nev);

            if(problem_is_complex)
            {


                MKLWrappers::DFtCalculator dft_calculator1(temp,temp_tilde);
                MKLWrappers::DFtCalculator dft_calculator2(temp2,temp2_tilde);

                // Declare vector data
                std::vector<std::complex<double>> resid(ldv);
                std::vector<std::complex<double>> V(ncv * ldv);
                std::vector<std::complex<double>> workd(3 * ldv);
                std::vector<std::complex<double>> workl(lworkl);
                std::vector<std::complex<double>> d(nev + 1);
                std::vector<std::complex<double>> workev(2 * ncv);
                std::vector<std::complex<double>> Z((ldv+1)*(nev+1));
                std::vector<double> rwork(ncv);
                std::complex<double> const sigma(0.0, 0.0);


                /////////////////////////////
                // Arpack's RCI start here
                /////////////////////////////

                while (ido != 99) {

                    arpack::naupd(ido,
                                  arpack::bmat::identity,
                                  2*N,
                                  arpack::which::smallest_magnitude,
                                  nev,
                                  tolerance,
                                  resid.data(),
                                  ncv,
                                  V.data(),
                                  ldv,
                                  iparam.data(),
                                  ipntr.data(),
                                  workd.data(),
                                  workl.data(),
                                  lworkl,
                                  rwork.data(),
                                  info);

                    // New u
#pragma omp parallel for
                    for (int i = 0; i < N; ++i) temp[i] = workd[ipntr[0]-1+i];
                    dft_calculator1.compute_forward();
#pragma omp parallel for
                    for (int i = 0; i < N; ++i) temp_tilde[i] *= 0.5 * std::pow(TWOPI,2) * kmod2[i];
                    dft_calculator1.compute_backward();
#pragma omp parallel for
                    for (int i = 0; i < N; ++i)
                        temp[i] += ( ( Vext[i] + 8*PI*scattering_length*std::norm(psi0[i]) - chemical_potential) *
                                workd[ipntr[0]-1+i]
                                + 4.0*PI*scattering_length * std::pow(psi0[i],2) * workd[ipntr[0]-1+N+i]);

                    // New v
#pragma omp parallel for
                    for (int i = 0; i < N; ++i) temp2[i] = workd[ipntr[0]-1+N+i];
                    dft_calculator2.compute_forward();
#pragma omp parallel for
                    for (int i = 0; i < N; ++i) temp2_tilde[i] *= 0.5 * std::pow(TWOPI,2) * kmod2[i];
                    dft_calculator2.compute_backward();
#pragma omp parallel for
                    for (int i = 0; i < N; ++i)
                        temp2[i] += ( ( Vext[i] + 8.0*PI*scattering_length * std::norm(psi0[i]) - chemical_potential ) *
                                workd[ipntr[0]-1+N+i]
                                      + 4.0*PI*scattering_length * std::pow(conj(psi0[i]),2) * workd[ipntr[0]-1+i]);


                    // Write the results into workd and give control back to arpack
#pragma omp parallel for
                    for (int i = 0; i < N; ++i)
                    {
                        workd[ipntr[1]-1+i] = temp[i];
                        workd[ipntr[1]-1+N+i] = (-temp2[i]);
                    }

                }

                // Recover the eigenvalues and, optionally, the eigenvectors
                std::vector<a_int> select(ncv);
                for (int i = 0; i < ncv; i++) select[i] = 1;
                arpack::neupd(rvec,
                              arpack::howmny::ritz_vectors,
                              select.data(),
                              d.data(),
                              Z.data(),
                              ldz,
                              sigma,
                              workev.data(),
                              arpack::bmat::identity,
                              2*N,
                              arpack::which::smallest_magnitude,
                              nev,
                              tolerance,
                              resid.data(),
                              ncv,
                              V.data(),
                              ldv,
                              iparam.data(),
                              ipntr.data(),
                              workd.data(),
                              workl.data(),
                              lworkl,
                              rwork.data(),
                              info);

                std::sort(d.begin(),d.end(),
                          [&] (std::complex<double> c1, std::complex<double> c2) {return c1.real() < c2.real();});

                for (int i = 0; i < nev; ++i) eigenvalues[i] = d[i];

                // TODO Eventually calculate the eigenvectors

            }
            else if (problem_is_real)
            {

                // Declare vector data
                std::vector<double> resid(ldv);
                std::vector<double> V(ncv * ldv);
                std::vector<double> workd(3 * ldv);
                std::vector<double> workl(lworkl);
                std::vector<double> d_r(nev + 1);
                std::vector<double> d_i(nev + 1);
                std::vector<double> workev(2 * ncv);
                std::vector<double> Z((ldv+1)*(nev+1));
                std::complex<double> const sigma(0.0, 0.0);

                std::vector<double> upv((N+1)*(nev+1));
                std::vector<double> umv((N+1)*(nev+1));
                MKLWrappers::DFtCalculator dft_calculator(temp,temp_tilde);
                std::vector<std::complex<double>> eigenvalues_upv(nev);
                std::vector<std::complex<double>> eigenvalues_umv(nev);

                /////////////////////////////////////////////////////////////////////////////////////////////
                // First use of Arpack RCI.
                // Solve HX(u+v)=lambda^2(u+v)
                // This will find the eigenvalues and, if requested, also the first step of the eigenvectors,
                // corresponding to (u+v)
                /////////////////////////////////////////////////////////////////////////////////////////////

                while (ido != 99)
                {

                    arpack::naupd(ido,
                                  arpack::bmat::identity,
                                  N,
                                  arpack::which::smallest_magnitude,
                                  nev,
                                  tolerance,
                                  resid.data(),
                                  ncv,
                                  V.data(),
                                  ldv,
                                  iparam.data(),
                                  ipntr.data(),
                                  workd.data(),
                                  workl.data(),
                                  lworkl,
                                  info);

                    // Apply X
                    for (int i = 0; i < N; ++i) temp[i] = workd[ipntr[0]-1+i];
                    dft_calculator.compute_forward();
#pragma omp parallel for
                    for (int i = 0; i < N; ++i) temp_tilde[i] *= 0.5 * std::pow(TWOPI,2) * kmod2[i];
                    dft_calculator.compute_backward();
#pragma omp parallel for
                    for (int i = 0; i < N; ++i)
                        temp[i] += ( Vext[i]
                                     + 12.0*PI*scattering_length * std::norm(psi0[i])
                                     - chemical_potential ) * workd[ipntr[0]-1+i];

                    // Apply H
#pragma omp parallel for
                    for (int i = 0; i < N; ++i) temp2[i] = temp[i];
                    dft_calculator.compute_forward();
#pragma omp parallel for
                    for (int i = 0; i < N; ++i) temp_tilde[i] *= 0.5 * std::pow(TWOPI,2) * kmod2[i];
                    dft_calculator.compute_backward();
#pragma omp parallel for
                    for (int i = 0; i < N; ++i)
                        temp[i] += ( Vext[i]
                                     + 4.0*PI*scattering_length * std::norm(psi0[i])
                                     - chemical_potential ) * temp2[i];

                    // Write the result into workd and give control back to arpack
                    for (int i = 0; i < N; ++i) workd[ipntr[1]-1+i] = temp[i].real();

                }

                // Recover the eigenvalues and, optionally, the eigenvectors
                std::vector<a_int> select(ncv);
                for (int i = 0; i < ncv; i++) select[i] = 1;
                arpack::neupd(rvec,
                              arpack::howmny::ritz_vectors,
                              select.data(),
                              d_r.data(),
                              d_i.data(),
                              upv.data(),
                              ldz,
                              sigma.real(),
                              sigma.imag(),
                              workev.data(),
                              arpack::bmat::identity,
                              N,
                              arpack::which::smallest_magnitude,
                              nev,
                              tolerance,
                              resid.data(),
                              ncv,
                              V.data(),
                              ldv,
                              iparam.data(),
                              ipntr.data(),
                              workd.data(),
                              workl.data(),
                              lworkl,
                              info);

                // Extract the energies of the Bogolyubov modes from the eigenvalues of HX and sort them in
                // ascending order (by real part)
                std::vector<std::pair<std::complex<double>,int>> pairs_upv;
                for (int i = 0; i < nev; ++i)
                {
                    eigenvalues_upv[i] = std::sqrt(std::complex<double>(d_r[i],d_i[i]));
                    pairs_upv.emplace_back(eigenvalues_upv[i],i);
                }

                std::sort(pairs_upv.begin(),pairs_upv.end(),
                          [&] (std::pair<std::complex<double>,int> const p1,
                               std::pair<std::complex<double>,int> const p2)
                          { return ( (p1.first.real() < p2.first.real()) ); }
                );

                for (int i = 0; i < nev; ++i) eigenvalues_upv[i] = pairs_upv[i].first;

                // If eigenvectors are requested, we need to solve an additional eigenvalue problem in order to find
                // (u-v).

                if(eigenvectors_required)
                {

                    /////////////////////////////////////////////////////////////////////////////////////////////
                    // Second use of Arpack RCI.
                    // Solve XH(u-v)=lambda^2(u-v)
                    /////////////////////////////////////////////////////////////////////////////////////////////

                    // Reset the control variables
                    ido = 0;

                    // Restarts RCI
                    while (ido != 99)
                    {

                        arpack::naupd(ido,
                                      arpack::bmat::identity,
                                      N,
                                      arpack::which::smallest_magnitude,
                                      nev,
                                      tolerance,
                                      resid.data(),
                                      ncv,
                                      V.data(),
                                      ldv,
                                      iparam.data(),
                                      ipntr.data(),
                                      workd.data(),
                                      workl.data(),
                                      lworkl,
                                      info);

                        // Apply H
                        for (int i = 0; i < N; ++i) temp[i] = workd[ipntr[0]-1+i];
                        dft_calculator.compute_forward();
#pragma omp parallel for
                        for (int i = 0; i < N; ++i) temp_tilde[i] *= 0.5 * std::pow(TWOPI,2) * kmod2[i];
                        dft_calculator.compute_backward();

#pragma omp parallel for
                        for (int i = 0; i < N; ++i)
                            temp[i] += ( Vext[i]
                                         + 4.0*PI*scattering_length * std::norm(psi0[i])
                                         - chemical_potential ) * workd[ipntr[0]-1+i];

                        // Apply X
#pragma omp parallel for
                        for (int i = 0; i < N; ++i) temp2[i] = temp[i];
                        dft_calculator.compute_forward();
#pragma omp parallel for
                        for (int i = 0; i < N; ++i) temp_tilde[i] *= 0.5 * std::pow(TWOPI,2) * kmod2[i];
                        dft_calculator.compute_backward();
#pragma omp parallel for
                        for (int i = 0; i < N; ++i)
                            temp[i] += ( Vext[i]
                                         + 12.0*PI*scattering_length * std::norm(psi0[i])
                                         - chemical_potential ) * temp2[i];

                        // Write the result into workd and give control back to arpack
                        for (int i = 0; i < N; ++i) workd[ipntr[1]-1+i] = temp[i].real();

                    }

                    // Recover the eigenvalues and eigenvectors
                    for (int i = 0; i < ncv; i++) select[i] = 1;
                    arpack::neupd(rvec,
                                  arpack::howmny::ritz_vectors,
                                  select.data(),
                                  d_r.data(),
                                  d_i.data(),
                                  umv.data(),
                                  ldz,
                                  sigma.real(),
                                  sigma.imag(),
                                  workev.data(),
                                  arpack::bmat::identity,
                                  N,
                                  arpack::which::smallest_magnitude,
                                  nev,
                                  tolerance,
                                  resid.data(),
                                  ncv,
                                  V.data(),
                                  ldv,
                                  iparam.data(),
                                  ipntr.data(),
                                  workd.data(),
                                  workl.data(),
                                  lworkl,
                                  info);

                    // Extract the energies of the Bogolyubov modes from the eigenvalues of HX and sort them in
                    // ascending order (by real part)
                    std::vector<std::pair<std::complex<double>,int>> pairs_umv;
                    for (int i = 0; i < nev; ++i)
                    {
                        if(d_r[i] <= 0)
                            d_r[i] = 0.0;
                        eigenvalues_umv[i] = std::sqrt(std::complex<double> (d_r[i],d_i[i]));
                        pairs_umv.emplace_back(eigenvalues_umv[i],i);

                    }

                    std::sort(pairs_umv.begin(),pairs_umv.end(),
                              [&] (std::pair<std::complex<double>,int> const p1,
                                   std::pair<std::complex<double>,int> const p2)
                              { return ( (p1.first.real() < p2.first.real()) ); }
                    );
                    for (int i = 0; i < nev; ++i) eigenvalues_umv[i] = pairs_umv[i].first;

                    // Vectors are re-initialized, now we can recover u and v
                    for (int i = 0; i < nev; ++i)
                    {

                        // Find the index corresponding to a certain eigenvalue
                        int k_upv = pairs_upv[i].second;
                        int k_umv = pairs_umv[i].second;

                        // Recover u and v

                        for (int j = 0; j < N; ++j)
                        {
                            u[i](j) = 0.5 * (upv[(N+1)*k_upv+j]+umv[(N+1)*k_umv+j]);
                            v[i](j) = 0.5 * (upv[(N+1)*k_upv+j]-umv[(N+1)*k_umv+j]);
                        }

                        // Normalize
                        double bogolyubov_norm = 0.0;
                        for (int j = 0; j < N; ++j)
                            bogolyubov_norm += (std::norm(u[i](j)) - std::norm(v[i](j)));
                        bogolyubov_norm *= dv;
                        for (int j = 0; j < N; ++j)
                        {
                            u[i](j) = u[i](j)/std::sqrt(std::abs(bogolyubov_norm));
                            v[i](j) = v[i](j)/std::sqrt(std::abs(bogolyubov_norm));
                        }

                    }

                }

                for (int i = 0; i < nev; ++i)
                {
                    eigenvalues[i] = eigenvalues_upv[i];
                }

            }

            return std::make_tuple(std::move(eigenvalues),std::move(u),std::move(v));

        }
    }
}