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
 * @file example-2.cpp. Full documentation at \subpage example-2.
 * @brief Elementary excitations via Bogolyubov equations in a three-dimensional, harmonically trapped Bose gas.
 * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
 *
 */

#include "UltraCold.hpp"

using namespace UltraCold;

int main() {

    Tools::InputParser ip("example-2.prm");

    ip.read_input_file();

    double xmax = ip.retrieve_double("xmax");
    double ymax = ip.retrieve_double("ymax");

    const int nx = ip.retrieve_int("nx");
    const int ny = ip.retrieve_int("ny");

    double       scattering_length   = ip.retrieve_double("scattering length");
    const int    number_of_particles = ip.retrieve_int("number of particles");
    const double atomic_mass         = ip.retrieve_double("atomic mass");
    double omegax                    = ip.retrieve_double("omegax");
    double omegay                    = ip.retrieve_double("omegay");

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

    const double omega_ho = std::sqrt(omegax*omegay);

    omegax = omegax/omega_ho;
    omegay = omegay/omega_ho;

    const double a_ho = std::sqrt(hbar/(atomic_mass*omega_ho));

    scattering_length   *= bohr_radius/a_ho;

    xmax = xmax/a_ho;
    ymax = ymax/a_ho;

    Vector<double> x(nx);
    Vector<double> y(ny);

    double dx = 2.*xmax/nx;
    double dy = 2.*ymax/ny;

    for (size_t i = 0; i < nx; ++i) x(i) = -xmax + i*dx;
    for (size_t i = 0; i < ny; ++i) y(i) = -ymax + i*dy;

    double dv = dx*dy;

    Vector<std::complex<double>> psi(nx,ny);
    Vector<double> Vext(nx,ny);

    for (size_t i = 0; i < nx; ++i)
        for (size_t j = 0; j < ny; ++j)
            {

                psi(i,j)  = exp(- (pow(x(i),2) + pow(y(j),2)) );

                Vext(i,j) = 0.5*( std::pow(omegax,2)*pow(x(i),2) +
                                  std::pow(omegay,2)*pow(y(j),2));

            }

    double norm = 0.0;
    for (size_t i = 0; i < psi.size(); ++i) norm += std::norm(psi[i]);
    norm *= dv;
    for (size_t i = 0; i < psi.size(); ++i) psi[i] *= std::sqrt(number_of_particles/norm);

    GPSolvers::GPSolver gp_solver(x,
                                  y,
                                  psi,
                                  Vext,
                                  scattering_length);


    double chemical_potential;
    std::tie(psi,chemical_potential) = gp_solver.run_gradient_descent(number_of_gradient_descent_steps,
                                                                      residual,
                                                                      alpha,
                                                                      beta,
                                                                      std::cout,
                                                                      10);

    GraphicOutput::DataWriter output_wave_function;
    output_wave_function.set_output_name("ground_state_wave_function");
    output_wave_function.write_vtk(x,y,psi,"ground_state_wave_function","ASCII");

    Vector<double> psi0_real(nx,ny);
    for (int i = 0; i < psi.size() ; ++i) psi0_real[i] = psi[i].real();

    std::vector<std::complex<double>> eigenvalues(number_of_modes);
    std::vector< Vector<std::complex<double>> > u(number_of_modes),v(number_of_modes);

    BogolyubovSolvers::TrappedBogolyubovSolver bdg_solver(x,
                                                          y,
                                                          psi0_real,
                                                          Vext,
                                                          scattering_length,
                                                          chemical_potential,
                                                          number_of_modes,
                                                          tolerance,
                                                          maximum_number_arnoldi_iterations,
                                                          calculate_eigenvectors);

    std::tie(eigenvalues,u,v) =  bdg_solver.run();

    GraphicOutput::DataWriter                    output_fluctuations;
    std::vector<Vector< std::complex<double> >> density_fluctuations(number_of_modes);
    std::vector<Vector< std::complex<double> >> phase_fluctuations(number_of_modes);

    for (int i = 0; i < number_of_modes; ++i)
    {

        std::cout << eigenvalues[i].real() << " " << eigenvalues[i].imag() << std::endl;
        density_fluctuations[i].reinit(nx, ny);
        phase_fluctuations[i].reinit(nx, ny);
        for (int j = 0; j < nx * ny; ++j)
        {
            density_fluctuations[i](j) = (u[i](j) + v[i](j)) * psi0_real(j);
            phase_fluctuations[i](j)   = (u[i](j) - v[i](j)) / psi0_real(j);
        }

        output_fluctuations.set_output_name("density_fluctuations_mode_" + std::to_string(i));
        output_fluctuations.write_vtk(x, y, density_fluctuations[i], "density_fluctuations","ASCII");

        output_fluctuations.set_output_name("phase_fluctuations_mode_" + std::to_string(i));
        output_fluctuations.write_vtk(x, y, phase_fluctuations[i], "phase_fluctuations","ASCII");

    }

    return 0;

}