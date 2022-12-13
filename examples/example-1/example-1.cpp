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
 * @file example-1.cpp. Full documentation at \subpage example-1. 
 * @brief Ground state and simple dynamics of a three-dimensional, harmonically trapped Bose gas.
 * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
 *
 */

#include "UltraCold.hpp"

using namespace UltraCold;

class myGPSolver : public GPSolvers::GPSolver
{
    public:

        using GPSolver::GPSolver;
        void run_operator_splitting(int number_of_time_steps,
                                    double time_step,
                                    double ramp_duration,
                                    double initial_scattering_length,
                                    double final_scattering_length,
                                    std::ostream& output_stream) ;

        void write_operator_splitting_output(size_t iteration_number,
                                             double current_scattering_length,
                                             std::ostream& output_stream) override;

    protected:

        void solve_step_1_operator_splitting(double current_scattering_length) override;

};

void myGPSolver::run_operator_splitting(int number_of_time_steps,
                                        double time_step,
                                        double ramp_duration,
                                        double initial_scattering_length,
                                        double final_scattering_length,
                                        std::ostream &output_stream)
{

    // Initialize the member variable time_step

    this->time_step = time_step;

    // Since the G.P. equation is solved on a cartesian mesh with periodic boundary conditions, a
    // DFtCalculator is needed to calculate the laplacian of psi

    MKLWrappers::DFtCalculator dft_calculator_step_2(psi,psitilde);

    //----------------------------------------------------//
    //    Here the operator-splitting iterations start    //
    //----------------------------------------------------//

    double current_scattering_length=initial_scattering_length;
    double current_time=0;

    for (size_t iteration_number = 0; iteration_number < number_of_time_steps; ++iteration_number)
    {

        // Write outputs starting from the first time step

        write_operator_splitting_output(iteration_number,
                                        current_scattering_length,
                                        output_stream);

        // Update the current value of the scattering length

        current_time = iteration_number*time_step;
        if(current_time <= ramp_duration)
        {

            current_scattering_length = initial_scattering_length
                    + (final_scattering_length-initial_scattering_length) * current_time/ramp_duration;
        }

        // Solve step 1 of operator splitting

        solve_step_1_operator_splitting(current_scattering_length);

        // Solve step 2 of operator splitting

        solve_step_2_operator_splitting(dft_calculator_step_2);

    }

}

void myGPSolver::solve_step_1_operator_splitting(double current_scattering_length)
{

        for (size_t i = 0; i < psi.size(); ++i)
            psi(i) *= std::exp(-ci*time_step*(Vext(i)+ 4*PI*current_scattering_length*std::norm(psi(i))));

}

void myGPSolver::write_operator_splitting_output(size_t iteration_number,
                                                 double current_scattering_length,
                                                 std::ostream &output_stream)
{
    if(iteration_number % 100 == 0)
    {
        double r2m = 0.0;
        double norm = 0.0;

        for (size_t i = 0; i < psi.extent(0); ++i)
            for(size_t j = 0; j < psi.extent(1); ++j)
                for (size_t k = 0; k < psi.extent(2); ++k)
                {
                    r2m += (std::pow(x[i],2)+std::pow(y[j],2)+std::pow(z[k],2))*std::norm(psi(i,j,k));   
                    norm += std::norm(psi(i,j,k));         
                }
        r2m = std::sqrt(r2m/norm);
        output_stream << iteration_number*time_step << " " << current_scattering_length << " " << r2m << std::endl;
    }
}

int main() {

    Tools::InputParser ip("example-1.prm");

    ip.read_input_file();

    double xmax = ip.retrieve_double("xmax");
    double ymax = ip.retrieve_double("ymax");
    double zmax = ip.retrieve_double("zmax");

    const int nx = ip.retrieve_int("nx");
    const int ny = ip.retrieve_int("ny");
    const int nz = ip.retrieve_int("nz");

    double initial_scattering_length = ip.retrieve_double("initial scattering length");
    const int    number_of_particles = ip.retrieve_int("number of particles");
    const double atomic_mass         = ip.retrieve_double("atomic mass");
    double omegax                    = ip.retrieve_double("omegax");
    double omegay                    = ip.retrieve_double("omegay");
    double omegaz                    = ip.retrieve_double("omegaz");

    const int    number_of_gradient_descent_steps = ip.retrieve_int("number of gradient descent steps");
    const double residual                         = ip.retrieve_double("residual");
    const double alpha                            = ip.retrieve_double("alpha");
    const double beta                             = ip.retrieve_double("beta");

    const int    number_of_real_time_steps = ip.retrieve_int("number of real time steps");
    double time_step                       = ip.retrieve_double("time step");
    double final_scattering_length         = ip.retrieve_double("final scattering length");
    double ramp_duration                   = ip.retrieve_double("ramp duration");

    const double hbar        = 0.6347*1.E5;
    const double bohr_radius = 5.292E-5;

    omegax *= TWOPI;
    omegay *= TWOPI;
    omegaz *= TWOPI;

    const double omega_ho = std::cbrt(omegax*omegay*omegaz);

    time_step     = time_step*omega_ho/1000.0;
    ramp_duration = ramp_duration*omega_ho/1000.0;

    omegax = omegax/omega_ho;
    omegay = omegay/omega_ho;
    omegaz = omegaz/omega_ho;

    const double a_ho = std::sqrt(hbar/(atomic_mass*omega_ho));

    initial_scattering_length *= bohr_radius/a_ho;
    final_scattering_length   *= bohr_radius/a_ho;

    xmax = xmax/a_ho;
    ymax = ymax/a_ho;
    zmax = zmax/a_ho;

    Vector<double> x(nx);
    Vector<double> y(ny);
    Vector<double> z(nz);

    double dx = 2.*xmax/nx;
    double dy = 2.*ymax/ny;
    double dz = 2.*zmax/nz;

    for (size_t i = 0; i < nx; ++i) x(i) = -xmax + i*dx;
    for (size_t i = 0; i < ny; ++i) y(i) = -ymax + i*dy;
    for (size_t i = 0; i < nz; ++i) z(i) = -zmax + i*dz;

    double dv = dx*dy*dz;

    Vector<std::complex<double>> psi(nx,ny,nz);
    Vector<double> Vext(nx,ny,nz);

    for (size_t i = 0; i < nx; ++i)
        for (size_t j = 0; j < ny; ++j)
            for (size_t k = 0; k < nz; ++k)
            {
                psi(i,j,k)  = exp(- (pow(x(i),2) +
                                     pow(y(j),2) +
                                     pow(z(k),2)) );

                Vext(i,j,k) = 0.5*( std::pow(omegax,2)*pow(x(i),2) +
                                    std::pow(omegay,2)*pow(y(j),2) +
                                    std::pow(omegaz,2)*pow(z(k),2) );
            }

    double norm = 0.0;
    for (size_t i = 0; i < psi.size(); ++i) norm += std::norm(psi[i]);
    norm *= dv;
    for (size_t i = 0; i < psi.size(); ++i) psi[i] *= std::sqrt(number_of_particles/norm);

    myGPSolver gp_solver(x,y,z,psi,Vext,initial_scattering_length);

    std::fstream gradient_descent_output_stream;
    gradient_descent_output_stream.open("gradient_descent_output.csv",std::ios::out);

    double chemical_potential;
    std::tie(psi,chemical_potential) = gp_solver.run_gradient_descent(number_of_gradient_descent_steps,
                                                                      residual,
                                                                      alpha,
                                                                      beta,
                                                                      gradient_descent_output_stream,
                                                                      100);

    gradient_descent_output_stream.close();

    GraphicOutput::DataWriter data_out;
    data_out.set_output_name("ground_state_wave_function");
    data_out.write_vtk(x,y,z,psi,"psi","ASCII");

    gp_solver.reinit(Vext,psi);

    std::fstream output_file_stream;
    output_file_stream.open("real_time_output.csv",std::ios::out);

    gp_solver.run_operator_splitting(number_of_real_time_steps,
                                     time_step,
                                     ramp_duration,
                                     initial_scattering_length,
                                     final_scattering_length,
                                     output_file_stream);

    output_file_stream.close();

    return 0;

}