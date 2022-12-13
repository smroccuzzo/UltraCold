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
 * @file example-5.cpp. Full documentation at \subpage example-5.
 * @brief Simple vortex dynamics in a two-dimensional dipolar Bose gas.
 * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
 *
 */

#include "UltraCold.hpp"

using namespace UltraCold;

class Dipoles2d : public cudaSolvers::DipolarGPSolver
{
public:

    using DipolarGPSolver::DipolarGPSolver;

    void write_operator_splitting_output(size_t iteration_number,
                                         std::ostream& output_stream) override;

};

void Dipoles2d::write_operator_splitting_output(size_t iteration_number,
                                                std::ostream &output_stream)
{

    if(iteration_number % write_output_every == 0)
    {
        copy_out_wave_function();
        GraphicOutput::DataWriter data_out;
        data_out.set_output_name("psi"+std::to_string(iteration_number/write_output_every));
        data_out.write_vtk(x_axis,y_axis,wave_function_output,"psi","ASCII");
    }
}

int main() {

    Tools::InputParser ip("example-5.prm");

    ip.read_input_file();

    double xmax = ip.retrieve_double("xmax");
    double ymax = ip.retrieve_double("ymax");

    const int nx = ip.retrieve_int("nx");
    const int ny = ip.retrieve_int("ny");

    double scattering_length = ip.retrieve_double("scattering length");
    double dipolar_length    = ip.retrieve_double("dipolar length");
    const int    number_of_particles = ip.retrieve_int("number of particles");
    const double atomic_mass         = ip.retrieve_double("atomic mass");
    double omegaz                    = ip.retrieve_double("omegaz");
    double theta = ip.retrieve_double("theta");

    const int    number_of_gradient_descent_steps = ip.retrieve_int("number of gradient descent steps");
    const double alpha                            = ip.retrieve_double("alpha");
    const double beta                             = ip.retrieve_double("beta");

    const int write_output_every=ip.retrieve_int("write output every");

    const int    number_of_real_time_steps = ip.retrieve_int("number of real time steps");
    double time_step                       = ip.retrieve_double("time step");

    // These two constants are for fixing the units
    const double hbar        = 0.6347*1.E5; // hbar in amu*mum^2/s
    const double bohr_radius = 5.292E-5;    // bohr radius in mum

    // Lengths are measured in units of the harmonic oscillator length along the z-axis, times as 1/(2 PI omega_z)
    omegaz *= TWOPI;
    time_step     = time_step*omegaz/1000.0;

    const double a_ho = std::sqrt(hbar/(atomic_mass*(omegaz)));

    scattering_length *= bohr_radius/a_ho;
    dipolar_length    *= bohr_radius/a_ho;

    xmax = xmax/a_ho;
    ymax = ymax/a_ho;
    Vector<double> x(nx);
    Vector<double> y(ny);
    double dx = 2.*xmax/nx;
    double dy = 2.*ymax/ny;
    for (size_t i = 0; i < nx; ++i) x(i) = -xmax + i*dx;
    for (size_t i = 0; i < ny; ++i) y(i) = -ymax + i*dy;
    double dv = dx*dy;

    double density = number_of_particles/(4*xmax*ymax);
    double sqrt_density=std::sqrt(density);

    Vector<std::complex<double>> psi(nx,ny),psitilde(nx,ny);
    Vector<double> Vext(nx,ny);

    double x_vortex1=-1/a_ho+xmax/2;
    double y_vortex1=0./a_ho+ymax/2;
    double x_vortex2=1./a_ho+xmax/2;
    double y_vortex2=0./a_ho+ymax/2;
    double phase;

    std::complex<double> ci={0.0,1.0};
    for (size_t i = 0; i < nx; ++i)
        for (size_t j = 0; j < ny; ++j)
        {

            psi(i,j) = sqrt(density);

            // Imprint a vortex at x1,y1
            phase=atan2(y(j)-y_vortex1,x(i)-x_vortex1);
            psi(i,j) *= sqrt_density*exp(ci*phase);

            // Imprint the corresponding anti-vortices at -x1,y1 and x1,-y1, and a vortex at -x1, -y1

            // -x1,y1
            phase=atan2(y(j)-y_vortex1,x(i)+x_vortex1);
            psi(i,j) *= sqrt_density*exp(-ci*phase);

            // x1,-y1
            phase=atan2(y(j)+y_vortex1,x(i)-x_vortex1);
            psi(i,j) *= sqrt_density*exp(-ci*phase);

            // -x1,-y1
            phase=atan2(y(j)+y_vortex1,x(i)+x_vortex1);
            psi(i,j) *= sqrt_density*exp(ci*phase);

            // Now the same at x2,y2

            phase=atan2(y(j)-y_vortex2,x(i)-x_vortex2);
            psi(i,j) *= sqrt_density*exp(ci*phase);

            // -x2,y2
            phase=atan2(y(j)-y_vortex2,x(i)+x_vortex2);
            psi(i,j) *= sqrt_density*exp(-ci*phase);

            // x2,-y2
            phase=atan2(y(j)+y_vortex2,x(i)-x_vortex2);
            psi(i,j) *= sqrt_density*exp(-ci*phase);

            // -x2,-y2
            phase=atan2(y(j)+y_vortex2,x(i)+x_vortex2);
            psi(i,j) *= sqrt_density*exp(ci*phase);

            // Uncomment if you want to initialize a hard-wall box at the boundaries,
            // or to set the external potential that you want
//            if(i == 0 || j == 0 || i == nx-1 || j == ny-1 )
//                Vext(i,j) = 100;
        }

    double norm = 0.0;
    for (size_t i = 0; i < psi.size(); ++i) norm += std::norm(psi[i]);
    norm *= dv;
    for (size_t i = 0; i < psi.size(); ++i) psi[i] *= std::sqrt(number_of_particles/norm);

    Dipoles2d gp_solver(x,y,psi,Vext,scattering_length,dipolar_length,theta);

    std::fstream gradient_descent_output_stream;
    gradient_descent_output_stream.open("gradient_descent_output.csv",std::ios::out);
    double chemical_potential;
    std::tie(psi,chemical_potential) = gp_solver.run_gradient_descent(number_of_gradient_descent_steps,
                                                                      alpha,
                                                                      beta,
                                                                      std::cout,
                                                                      write_output_every);

    gradient_descent_output_stream.close();

    gp_solver.reinit(Vext,psi);
    gp_solver.run_operator_splitting(number_of_real_time_steps,time_step,std::cout,write_output_every);

    return 0;

}
