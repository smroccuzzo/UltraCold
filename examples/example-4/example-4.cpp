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

/**
 *
 * @file example-4.cpp. Full documentation at \subpage example-4.
 * @brief Excitation spectrum of a trapped dipolar Bose-Einstein condensate across the superfluid-supersolid phase transition.
 * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
 *
 */

#include "UltraCold.hpp"
#include <random>

using namespace UltraCold;

int main()
{

    Tools::InputParser ip("example-4.prm");

    ip.read_input_file();

    double xmax = ip.retrieve_double("xmax");
    double ymax = ip.retrieve_double("ymax");
    double zmax = ip.retrieve_double("zmax");

    const int nx = ip.retrieve_int("nx");
    const int ny = ip.retrieve_int("ny");
    const int nz = ip.retrieve_int("nz");

    double scattering_length         = ip.retrieve_double("scattering length");
    double dipolar_length            = ip.retrieve_double("dipolar_length");
    const int    number_of_particles = ip.retrieve_int("number of particles");
    const double atomic_mass         = ip.retrieve_double("atomic mass");
    double omegax                    = ip.retrieve_double("omegax");
    double omegay                    = ip.retrieve_double("omegay");
    double omegaz                    = ip.retrieve_double("omegaz");
    double theta_mu                  = ip.retrieve_double("theta");
    double phi_mu                    = ip.retrieve_double("phi");

    const int    number_of_gradient_descent_steps = ip.retrieve_int("number of gradient descent steps");
    const double residual                         = ip.retrieve_double("residual");
    const double alpha                            = ip.retrieve_double("alpha");
    const double beta                             = ip.retrieve_double("beta");

    const int number_of_modes = ip.retrieve_int("number of modes");
    const int maximum_number_arnoldi_iterations = ip.retrieve_int("maximum number of arnoldi iterations");
    const double tolerance = ip.retrieve_double("tolerance");
    const bool calculate_eigenvectors = ip.retrieve_bool("calculate eigenvectors");

    const double hbar        = 0.6347*1.E5;
    const double bohr_radius = 5.292E-5;

    omegax *= TWOPI;
    omegay *= TWOPI;
    omegaz *= TWOPI;

    const double omega_ho = std::cbrt(omegax*omegay*omegaz);

    omegax = omegax/omega_ho;
    omegay = omegay/omega_ho;
    omegaz = omegaz/omega_ho;

    const double a_ho = std::sqrt(hbar/(atomic_mass*omega_ho));

    scattering_length *= bohr_radius/a_ho;
    dipolar_length *= bohr_radius/a_ho;

    xmax = xmax/a_ho;
    ymax = ymax/a_ho;
    zmax = zmax/a_ho;

    double dx = 2 * xmax / nx;
    double dy = 2 * ymax / ny;
    double dz = 2 * zmax / nz;

    Vector<double> x(nx), y(ny), z(nz), kx(nx), ky(ny), kz(nz);

    for (int i = 0; i < nx; ++i) x[i] = -xmax + i * dx;
    for (int i = 0; i < ny; ++i) y[i] = -ymax + i * dy;
    for (int i = 0; i < nz; ++i) z[i] = -zmax + i * dz;
    create_mesh_in_Fourier_space(x, y, z, kx, ky, kz);

    Vector<std::complex<double>> psi(nx, ny, nz);
    Vector<double> Vext(nx, ny, nz);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0,1);

    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
            {
                double random_number = distribution(generator);
                psi(i,j,k)  = (1.0+0.1*random_number)*
                              std::exp(-0.1*(pow(x(i),2) +
                                             pow(y(j),2) +
                                             pow(z(k),2)) );

                Vext(i,j,k) = 0.5*( std::pow(omegax,2)*pow(x(i),2) +
                                    std::pow(omegay,2)*pow(y(j),2) +
                                    std::pow(omegaz,2)*pow(z(k),2) );
            }

    double norm = 0.0;
    for (size_t i = 0; i < psi.size(); ++i) norm += std::norm(psi[i]);
    norm *= (dx * dy * dz);
    for (size_t i = 0; i < psi.size(); ++i) psi[i] *= std::sqrt(number_of_particles / norm);

    GraphicOutput::DataWriter psi_out;
    psi_out.set_output_name("initial_wave_function");
    psi_out.write_slice2d_vtk(x,y,psi,"xy","initial_wave_function","ASCII");

    GPSolvers::DipolarGPSolver dipolar_gp_solver(x,
                                                 y,
                                                 z,
                                                 psi,
                                                 Vext,
                                                 scattering_length,
                                                 dipolar_length,
                                                 theta_mu,
                                                 phi_mu);
    double chemical_potential;
    std::tie(psi, chemical_potential) = dipolar_gp_solver.run_gradient_descent(number_of_gradient_descent_steps,
                                                                               residual,
                                                                               alpha,
                                                                               beta,
                                                                               std::cout,
                                                                               10);

    psi_out.set_output_name("ground_state_wave_function");
    psi_out.write_vtk(x,y,z,psi,"ground_state_wave_function","BINARY");

    Vector<double> psi_real(nx,ny,nz);
    for (int i = 0; i < nx * ny * nz; ++i)
        psi_real[i] = psi[i].real();

    BogolyubovSolvers::TrappedDipolarBogolyubovSolver bogolyubov_solver(x,
                                                                        y,
                                                                        z,
                                                                        psi_real,
                                                                        Vext,
                                                                        scattering_length,
                                                                        dipolar_length,
                                                                        theta_mu,
                                                                        phi_mu,
                                                                        chemical_potential,
                                                                        number_of_modes,
                                                                        tolerance,
                                                                        maximum_number_arnoldi_iterations,
                                                                        calculate_eigenvectors);

    std::vector<std::complex<double>> eigenvalues(number_of_modes);
    std::vector<Vector<std::complex<double>>> u(number_of_modes),v(number_of_modes);
    std::tie(eigenvalues,u,v) = bogolyubov_solver.run();

    GraphicOutput::DataWriter                    output_fluctuations;
    std::vector<Vector< std::complex<double> >> density_fluctuations(number_of_modes);
    std::vector<Vector< std::complex<double> >> phase_fluctuations(number_of_modes);

    std::cout << "---------------------------" << std::endl;
    std::cout << "expected sloshing frequencies: " << std::endl;
    std::cout << omegax << " " << omegay << " " << omegaz << std::endl;
    std::cout << "---------------------------" << std::endl;

    for (int i = 0; i < number_of_modes; ++i)
    {

        std::cout << eigenvalues[i].real() << " " << eigenvalues[i].imag() << std::endl;
	if(calculate_eigenvectors)
	{

        	density_fluctuations[i].reinit(nx,ny,nz);
        	phase_fluctuations[i].reinit(nx,ny,nz);
        	for (int j = 0; j < nx*ny*nz; ++j)
        	{
            		density_fluctuations[i](j) = (u[i](j) + v[i](j)) * psi_real(j);
            		phase_fluctuations[i](j)   = (u[i](j) - v[i](j)) / psi_real(j);
        	}

        	output_fluctuations.set_output_name("density_fluctuations_mode_" + std::to_string(i));
        	output_fluctuations.write_vtk(x,y,z,density_fluctuations[i], "density_fluctuations","BINARY");

        	output_fluctuations.set_output_name("phase_fluctuations_mode_" + std::to_string(i));
        	output_fluctuations.write_vtk(x,y,z,phase_fluctuations[i], "phase_fluctuations","BINARY");

	}

    }

    return 0;

}
