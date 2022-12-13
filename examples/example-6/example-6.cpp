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
 * @file example-6.cpp. Full documentation at \subpage example-6.
 * @brief Far from equilibrium dynamics of a dipolar BEC in 2D, using GPU acceleration.
 * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
 *
 */

#include "UltraCold.hpp"
#include <iomanip>
#include <cmath>

using namespace UltraCold;

/////////////////////////////////////////////////////////
// Our solver class, with customizable real-time output
/////////////////////////////////////////////////////////

class Dipoles2d : public cudaSolvers::DipolarGPSolver
{
public:

    using DipolarGPSolver::DipolarGPSolver;

    void write_operator_splitting_output(size_t iteration_number,
                                         std::ostream& output_stream) override;

    MKLWrappers::DFtCalculator ft_output_calculator;
    Vector<std::complex<double>> ft_wave_function;
    Vector<double> nbar_x;
    Vector<double> nbar_y;

};

void Dipoles2d::write_operator_splitting_output(size_t iteration_number,
                                                std::ostream &output_stream)
{

    if(iteration_number == 0)
    {
        ft_wave_function.reinit(nx,ny);
        nbar_x.reinit(nx);
        nbar_y.reinit(ny);
        ft_output_calculator.reinit(wave_function_output,ft_wave_function);
    }
    if(iteration_number % write_output_every == 0)
    {
        copy_out_wave_function();
        GraphicOutput::DataWriter psi_out;
        psi_out.set_output_name("psi"+std::to_string(iteration_number/write_output_every));
        psi_out.write_vtk(x_axis,y_axis,wave_function_output,"psi","BINARY");

        // Calculate momentum distribution
        ft_output_calculator.compute_forward();

        // Reset nx and ny
        for(int i = 0; i < nx; ++i) nbar_x(i) = 0.0;
        for(int i = 0; i < ny; ++i) nbar_y(i) = 0.0;

        // Calculate averages along x and y
        for(int i = 0; i < nx; ++i)
            for(int j = 0; j < ny; ++j)
            {
                nbar_x(i) += std::pow(std::abs(ft_wave_function(i,j)),2);
                nbar_y(j) += std::pow(std::abs(ft_wave_function(i,j)),2);
            }

        std::ofstream momentum_cuts;
        momentum_cuts.open("nkx"+std::to_string(iteration_number/write_output_every)+".txt",std::ios::out);
        for(size_t ix = 0; ix < int(nx/2+1); ++ix)
            momentum_cuts << kx_axis(ix) << " "
                          << nbar_x(ix)
                          << std::endl;
        momentum_cuts.close();

        momentum_cuts.open("nky"+std::to_string(iteration_number/write_output_every)+".txt",std::ios::out);
        for(size_t iy = 0; iy < int(ny/2+1); ++iy)
            momentum_cuts << ky_axis(iy) << " "
                          << nbar_y(iy)
                          << std::endl;

    }
}

int main() {

    // Read input parameters from file "dipolars2d.prm". This must be placed in the same directory as
    // the executable

    Tools::InputParser ip("example-6.prm");

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

    const int    number_of_real_time_steps = ip.retrieve_int("number of real time steps");
    double time_step                       = ip.retrieve_double("time step");

    const int write_output_every=ip.retrieve_int("write output every");

    const bool box_initial_conditions = ip.retrieve_bool("box initial conditions");
    const bool vortex_lattice_initial_conditions = ip.retrieve_bool("vortex lattice initial conditions");
    const bool vortex_random_distribution_initial_conditions =
            ip.retrieve_bool("vortex random distribution initial conditions");
    const int winding_number = ip.retrieve_int("winding number");
    const int number_of_imprinted_defects = ip.retrieve_int("number of imprinted defects");

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

    //////////////////////
    // Create the mesh
    /////////////////////

    Vector<double> x(nx);
    Vector<double> y(ny);

    double dx = 2.*xmax/nx;
    double dy = 2.*ymax/ny;

    for (size_t i = 0; i < nx; ++i) x(i) = -xmax + i*dx;
    for (size_t i = 0; i < ny; ++i) y(i) = -ymax + i*dy;

    double dv = dx*dy;

    //////////////////////////////
    // Create the momentum mesh
    //////////////////////////////

    Vector<double> kx(nx),ky(ny);
    create_mesh_in_Fourier_space(x,y,kx,ky);

    ///////////////////////////////////////////////////////////////////
    // Initialize wave function and external potential
    // Wave function is initialized in momentum space populating
    // low momentum modes each with a random phase
    ///////////////////////////////////////////////////////////////////

    double density = number_of_particles/(4*xmax*ymax);

    Vector<std::complex<double>> psi(nx,ny),psitilde(nx,ny);
    Vector<double> Vext(nx,ny);

    // Print out some useful parameters
    double epsilon_dd=dipolar_length/scattering_length;
    double expected_chemical_potential = sqrt(8*PI)*(scattering_length)*density*(1+epsilon_dd*(3*pow(cos(theta),2)-1));
    std::cout << "density = " << number_of_particles/(4*xmax*ymax*pow(a_ho,2)) << " mum^{-2}" << std::endl;
    std::cout << "expected chemical potential   = " << expected_chemical_potential << " hbar omegaz" << std::endl;
    double healing_length = 1./sqrt(std::abs(expected_chemical_potential));
    std::cout << "az= " << a_ho << " mum" << std::endl;
    std::cout << "healing_length    = " << healing_length*a_ho << " mum" << std::endl;
    std::cout << "mesh size along x = " << 2*xmax*a_ho << " mum" << std::endl;
    std::cout << "mesh size along y = " << 2*ymax*a_ho << " mum" << std::endl;
    std::cout << "step size along x = " << dx*a_ho << " mum" << std::endl;
    std::cout << "step size along y = " << dy*a_ho << " mum" << std::endl;
    std::cout << "natural time step = " << 1./(expected_chemical_potential * omegaz/1000) << " ms" << std::endl;
    std::cout << "kxi = " << 1./healing_length << std::endl;
    std::cout << "dipolar factor = " << (1+epsilon_dd*(3*pow(cos(theta),2)-1)) << std::endl;
    std::cout << "kxi nodipo = " << sqrt((1+epsilon_dd*(3*pow(cos(theta),2)-1))) * 1./healing_length << std::endl;
    std::cout << "kmax= " << 1./(2*dx) << std::endl;

    std::default_random_engine generator;

    std::complex<double> ci={0.0,1.0};

    if(box_initial_conditions)
    {

        typedef std::chrono::high_resolution_clock clock;
        clock::time_point beginning = clock::now();
        clock::duration d = clock::now() - beginning;
        generator.seed(d.count());

        std::uniform_real_distribution<double> phase_distribution(0,TWOPI);

        for (size_t i = 0; i < nx; ++i)
            for (size_t j = 0; j < ny; ++j)
            {
                double random_phase = phase_distribution(generator);
                if (std::abs(kx(i)) <= 0.1 / (TWOPI*healing_length) && std::abs(ky(j)) <= 0.1 /(TWOPI*healing_length))
                    psitilde(i, j) = exp(-ci * random_phase);
            }

        MKLWrappers::DFtCalculator dft(psi,psitilde);
        dft.compute_backward();

        double norm = 0.0;
        for (size_t i = 0; i < psi.size(); ++i) norm += std::norm(psi[i]);
        norm *= dv;
        for (size_t i = 0; i < psi.size(); ++i) psi[i] *= std::sqrt(number_of_particles/norm);

    }

    else if(vortex_lattice_initial_conditions)
    {
        double sqrt_density = std::sqrt(number_of_particles/(4*xmax*ymax));

        Vector<double> x_defect_positions(number_of_imprinted_defects);
        Vector<double> y_defect_positions(number_of_imprinted_defects);

        std::uniform_real_distribution<double> x_displacement(-2*dx,2*dx);
        std::uniform_real_distribution<double> y_displacement(-2*dy,2*dy);

        typedef std::chrono::high_resolution_clock clock;
        clock::time_point beginning = clock::now();
        clock::duration d = clock::now() - beginning;
        generator.seed(d.count());
        for(int i = 0; i < number_of_imprinted_defects; ++i)
            x_defect_positions(i) = -xmax+dx+i*nx/number_of_imprinted_defects*dx;

        beginning = clock::now();
        d = clock::now() - beginning;
        generator.seed(d.count());
        for(int i = 0; i < number_of_imprinted_defects; ++i)
            y_defect_positions(i) = -ymax+dy+i*ny/number_of_imprinted_defects*dy;

        for (size_t i = 0; i < nx; ++i)
            for (size_t j = 0; j < ny; ++j)
            {
                psi(i,j) = sqrt_density;

                for(int i_defect=0; i_defect < number_of_imprinted_defects; ++i_defect)
                    for(int j_defect=0; j_defect < number_of_imprinted_defects; ++j_defect)
                    {
                        double phase = atan2(y(j)-y_defect_positions(j_defect)+y_displacement(generator),
                                             x(i)-x_defect_positions(i_defect)+x_displacement(generator));
                        phase *= (i_defect%2 == 0 && j_defect%2 == 0 ||
                                  i_defect%2 != 0 && j_defect%2 != 0) ? winding_number : -winding_number;
                        psi(i,j) *= exp(ci*phase);
                    }
            }

        double norm = 0.0;
        for (size_t i = 0; i < psi.size(); ++i) norm += std::norm(psi[i]);
        norm *= dv;
        for (size_t i = 0; i < psi.size(); ++i) psi[i] *= std::sqrt(number_of_particles/norm);
    }

    else if(vortex_random_distribution_initial_conditions)
    {

        double sqrt_density = std::sqrt(number_of_particles/(4*xmax*ymax));

        Vector<double> x_defect_positions(number_of_imprinted_defects);
        Vector<double> y_defect_positions(number_of_imprinted_defects);

        std::uniform_real_distribution<double> x_defect_distribution(-xmax,xmax);
        std::uniform_real_distribution<double> y_defect_distribution(-ymax,ymax);

        typedef std::chrono::high_resolution_clock clock;
        clock::time_point beginning = clock::now();
        clock::duration d = clock::now() - beginning;
        generator.seed(d.count());
        for(int i = 0; i < number_of_imprinted_defects; ++i) x_defect_positions(i) = x_defect_distribution(generator);

        typedef std::chrono::high_resolution_clock clock;
        beginning = clock::now();
        d = clock::now() - beginning;
        generator.seed(d.count());
        for(int i = 0; i < number_of_imprinted_defects; ++i) y_defect_positions(i) = y_defect_distribution(generator);

        for (size_t i = 0; i < nx; ++i)
            for (size_t j = 0; j < ny; ++j)
            {
                psi(i,j) = sqrt_density;
                for(int i_defect=0; i_defect < number_of_imprinted_defects; ++i_defect)
                {
                    double phase = atan2(y(j)-y_defect_positions(i_defect),x(i)-x_defect_positions(i_defect));
                    phase *= (i_defect%2==0) ? winding_number : -winding_number;
                    psi(i,j) *= exp(ci*phase);
                }

            }

        double norm = 0.0;
        for (size_t i = 0; i < psi.size(); ++i) norm += std::norm(psi[i]);
        norm *= dv;
        for (size_t i = 0; i < psi.size(); ++i) psi[i] *= std::sqrt(number_of_particles/norm);

    }
    else
    {

        double sqrt_density = std::sqrt(number_of_particles/(4*xmax*ymax));

        std::uniform_real_distribution<double> random_distribution(0,1);

        typedef std::chrono::high_resolution_clock clock;
        clock::time_point beginning = clock::now();
        clock::duration d = clock::now() - beginning;
        generator.seed(d.count());

        for (size_t i = 0; i < nx; ++i)
            for (size_t j = 0; j < ny; ++j)
                psi(i,j) = (1.0+random_distribution(generator))*sqrt_density;

        double norm = 0.0;
        for (size_t i = 0; i < psi.size(); ++i) norm += std::norm(psi[i]);
        norm *= dv;
        for (size_t i = 0; i < psi.size(); ++i) psi[i] *= std::sqrt(number_of_particles/norm);

    }

    //////////////////////////////////////////////////
    // Write the initial wave function in a .vtk file
    //////////////////////////////////////////////////

    GraphicOutput::DataWriter data_out;
    data_out.set_output_name("initial_wave_function_gradient_descent");
    data_out.write_vtk(x,y,psi,"psi","ASCII");

    /////////////////////////////
    // Initialize the solver
    /////////////////////////////

    Dipoles2d gp_solver(x,y,psi,Vext,scattering_length,dipolar_length,theta);

    /////////////////////////////////////////////////////////////////////////////////
    // The output from the short imaginary time propagation will be written in the
    // file "gradient_descent_output.csv"
    /////////////////////////////////////////////////////////////////////////////////

    if(vortex_lattice_initial_conditions || vortex_random_distribution_initial_conditions)
    {
        std::fstream gradient_descent_output_stream;
        gradient_descent_output_stream.open("gradient_descent_output.csv",std::ios::out);
        double chemical_potential;
        std::tie(psi,chemical_potential) = gp_solver.run_gradient_descent(number_of_gradient_descent_steps,
                                                                          alpha,
                                                                          beta,
                                                                          std::cout,
                                                                          1);

        gradient_descent_output_stream.close();
        gp_solver.reinit(Vext,psi);

    }

    ////////////////////////////////////////////////
    // Run the real time simulation!
    // The output will be written to the shell
    ////////////////////////////////////////////////

    gp_solver.set_tw_initial_conditions(false);
    gp_solver.run_operator_splitting(number_of_real_time_steps,time_step,std::cout,write_output_every);

    return 0;

}