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

#include "arpack.hpp"
#include "BogolyubovSolvers.hpp"

#define PI 3.1415926535897932384626433
#define TWOPI (2*PI)

namespace UltraCold
{
    namespace BogolyubovSolvers
    {

        // Constructors for real problems

        /**
         * @brief Constructor for a TrappedDipolarBogolyubovSolver in three space dimensions and real stationary
         * condensate wave function
         *
         * @param x *Vector<double>* representing the x-axis on which the Bogolyubov equations will be solved
         * @param y *Vector<double>* representing the y-axis on which the Bogolyubov equations will be solved
         * @param z *Vector<double>* representing the z-axis on which the Bogolyubov equations will be solved
         * @param psi_0 *Vector<double>* representing a stationary condensate wave function
         * @param Vext *Vector<double>* representing the external potential.
         * @param scattering_length *double* the scattering length in appropriate units
         * @param dipolar_length *double* the dipolar length in appropriate units
         * @param theta_mu *double* angle between the the polarization direction and the z-axis
         * @param phi_mu *double* angle between the projection of the polarization direction in the x-y plane and
         * the x-axis
         * @param chemical_potential *double* the ground-state chemical potential
         * @param eigenvectors_required *bool* specifies if also the eigenvectors are required. Default is false.
         *
         */

        TrappedDipolarBogolyubovSolver::TrappedDipolarBogolyubovSolver(Vector<double>& x,
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
            this->eigenvectors_required=eigenvectors_required;
            this->tolerance = tolerance;
            this->maximum_number_of_arnoldi_iterations = maximum_number_arnoldi_iterations;

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
            for (int i = 0; i < nx*ny*nz; ++i)
            {
                this->psi0[i] = psi0[i];
            }

            // Initialize the Fourier transform of the dipolar potential

            epsilon_dd = 0.0;
            if(scattering_length != 0)
                epsilon_dd = dipolar_length/scattering_length;

            Vtilde.reinit(nx,ny,nz);
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k)
                    {
                        double aux = TWOPI * (
                                kx[i]*std::sin(theta_mu)*std::cos(phi_mu) +
                                ky[j]*std::sin(theta_mu)*std::sin(phi_mu)+
                                kz[k]*std::cos(theta_mu));
                        double aux1 = TWOPI * std::sqrt(std::pow(kx[i], 2) + std::pow(ky[j], 2) + std::pow(kz[k], 2));
                        if (aux1 <= 1.E-6)
                            Vtilde(i, j, k) = 0.0;
                        else
                            Vtilde(i, j, k) =
                                    12.0 * PI * scattering_length * epsilon_dd * (std::pow(aux/aux1,2)-1.0/3.0);
                    }

            // Initialize dipole-dipole potential

            Phi_dd.reinit(nx,ny,nz);
            Vector<std::complex<double>> Phi_tilde(nx,ny,nz);
            MKLWrappers::DFtCalculator dft_dd(Phi_dd,Phi_tilde);

#pragma omp parallel for
            for (int i = 0; i < nx*ny*nz; ++i)
                Phi_dd[i] = std::norm(this->psi0[i]);

            dft_dd.compute_forward();
#pragma omp parallel for
            for (int i = 0; i < nx*ny*nz; ++i)
                Phi_tilde[i] *= Vtilde[i];
            dft_dd.compute_backward();

            // Initialize gamma(\epsilon_dd) for the LHY correction

            gamma_epsilon_dd = 0.0;
            if (epsilon_dd != 0)
            {
                gamma_epsilon_dd = 64.0*std::sqrt(PI)/3*std::sqrt(std::pow(scattering_length,5));
                double F_epsilon_dd=0.0;
                int n_theta=1000;
                double d_theta=PI/(n_theta-1);

                std::complex<double> csum(0.0,0.0);
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
         * @brief Solve the Bogolyubov equations.
         *
         * This member function actually solves the Bogolyubov equations for a trapped dipolar condensate. In the case
         * in which a real wave-function was passed to the constructor, the function will use the following useful
         * recast of the problem in one of halved dimensionality. Starting from
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
         * and taking the sum and the difference between two equations, one obtain the following equivalent formulation
         * of the problem
         *
         * \f[
         *  \begin{align}
         *      & (\hat{H}-\mu)(\hat{H}-\mu+2\hat{X}) (u+v) = (\hbar\omega)^2 (u+v) \nonumber \\
         *      & (\hat{H}-\mu+2\hat{X})(\hat{H}-\mu) (u-v) = (\hbar\omega)^2 (u-v) \nonumber \\
         *  \end{align}
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
         * If the eigenvectors are not required, this function will solve only the first equation, calculating the
         * square of the energies of the eigen-modes, but *returning the energies* (i.e., it calculates the square
         * root before returning). If the eigenvectors are also requested, the function will solve also the second
         * eigenvalue problem, obtaining the Bogolyubov quasi-particle amplitudes \f$ u \f$ and \f$ v \f$ by taking
         * (half) the sum and the difference between the eigenvectors of the two problems.
         *
         * In the case in which, instead, the wave-function passed to the constructor is complex, the function will
         * solve, by brute force, the complete problem.
         *
         * */

        std::tuple<std::vector<std::complex<double>>,
                std::vector<Vector<std::complex<double>>>,
                std::vector<Vector<std::complex<double>>>> TrappedDipolarBogolyubovSolver::run()
        {

            ///////////////////////////////////
            // Initialize the variables
            ///////////////////////////////////

            a_int N;
            N=nx*ny*nz;

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

            /////////////////////////////
            // Arpack's RCI start here
            /////////////////////////////

            if (problem_is_real)
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

                MKLWrappers::DFtCalculator dft_calculator(temp, temp_tilde);
                MKLWrappers::DFtCalculator dft_calculator2(temp2, temp2_tilde);

                std::vector<std::complex<double>> eigenvalues_upv(nev);
                std::vector<std::complex<double>> eigenvalues_umv(nev);

                /////////////////////////////////////////////////////////////////////////////////////////////
                // First use of Arpack RCI.
                // Solve (H-mu)(H-mu+2X)(u+v)=lambda^2(u+v)
                // This will find the eigenvalues and, if requested, also the first step of the eigenvectors,
                // corresponding to (u+v)
                /////////////////////////////////////////////////////////////////////////////////////////////

                while (ido != 99) {

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

                    // Apply H-mu+2X

#pragma omp parallel for
                    for (int i = 0; i < N; ++i) temp2[i] = psi0[i]*workd[ipntr[0]-1+i];
                    dft_calculator2.compute_forward();
#pragma omp parallel for
                    for (int i = 0; i < N; ++i) temp2_tilde[i] *= Vtilde[i];
                    dft_calculator2.compute_backward();
#pragma omp parallel for
                    for (int i = 0; i < N; ++i) temp2[i] *= (2.0*psi0[i]);

                    for (int i = 0; i < N; ++i) temp[i] = workd[ipntr[0]-1+i];
                    dft_calculator.compute_forward();
#pragma omp parallel for
                    for (int i = 0; i < N; ++i) temp_tilde[i] *= 0.5 * std::pow(TWOPI,2) * kmod2[i];
                    dft_calculator.compute_backward();

#pragma omp parallel for
                    for (int i = 0; i < N; ++i)
                        temp[i] += ( Vext[i]
                                    + 12.0*PI*scattering_length * std::norm(psi0[i])
                                    + 4.0*gamma_epsilon_dd*std::pow(psi0[i],3)
                                    + Phi_dd[i]
                                    - chemical_potential ) *
                                   workd[ipntr[0]-1+i] + temp2[i];

                    // Apply H-mu

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
                                     + gamma_epsilon_dd*std::pow(std::abs(psi0[i]),3)
                                     + Phi_dd[i]
                                     - chemical_potential ) *
                                   temp2[i];

                    // Write the result into workd and give control back to arpack

                    for (int i = 0; i < N; ++i)
                        workd[ipntr[1]-1+i] = temp[i].real();

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

                // Transform the eigenvalues into the Bogolyubov modes, plus sort them in ascending order (by real
                // part)

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
                    // Solve (H-mu+2X)(H-mu)(u-v)=lambda^2(u-v)
                    /////////////////////////////////////////////////////////////////////////////////////////////

                    // Reset the control variables

                    ido = 0;

                    // Need another auxiliary vector

                    Vector<std::complex<double>> temp3(nx,ny,nz);

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

                        // Apply H-mu

                        for (int i = 0; i < N; ++i) temp[i] = workd[ipntr[0]-1+i];
                        dft_calculator.compute_forward();
#pragma omp parallel for
                        for (int i = 0; i < N; ++i) temp_tilde[i] *= 0.5 * std::pow(TWOPI,2) * kmod2[i];
                        dft_calculator.compute_backward();

#pragma omp parallel for
                        for (int i = 0; i < N; ++i)
                            temp[i] += ( Vext[i]
                                         + 4.0*PI*scattering_length * std::norm(psi0[i])
                                         + gamma_epsilon_dd*std::pow(std::abs(psi0[i]),3)
                                         + Phi_dd[i]
                                         - chemical_potential ) * workd[ipntr[0]-1+i];

                        // Apply H-mu+2X

#pragma omp parallel for
                        for (int i = 0; i < N; ++i) temp2[i] = psi0[i]*temp[i];
                        dft_calculator2.compute_forward();
#pragma omp parallel for
                        for (int i = 0; i < N; ++i) temp2_tilde[i] *= Vtilde[i];
                        dft_calculator2.compute_backward();
#pragma omp parallel for
                        for (int i = 0; i < N; ++i) temp2[i] *= (2.0*psi0[i]);

#pragma omp parallel for
                        for (int i = 0; i < N; ++i) temp3[i] = temp[i];

                        dft_calculator.compute_forward();
#pragma omp parallel for
                        for (int i = 0; i < N; ++i) temp_tilde[i] *= 0.5 * std::pow(TWOPI,2) * kmod2[i];
                        dft_calculator.compute_backward();
#pragma omp parallel for
                        for (int i = 0; i < N; ++i)
                            temp[i] += ( Vext[i]
                                         + 12.0*PI*scattering_length * std::norm(psi0[i])
                                         + 4.0*gamma_epsilon_dd*std::pow(psi0[i],3)
                                         + Phi_dd[i]
                                         - chemical_potential ) * temp3[i] + temp2[i];

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
                    eigenvalues[i] = eigenvalues_upv[i];
            }

            return std::make_tuple(std::move(eigenvalues),std::move(u),std::move(v));

        }
    }
}