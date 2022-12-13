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
 * @file example-7.cpp. Full documentation at \subpage example-7.
 * @brief Far from equilibrium dynamics of a dipolar BEC in 3D, using GPU acceleration.
 * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
 *
 */

#include "UltraCold.hpp"
#include <iomanip>
#include <cmath>

using namespace UltraCold;

////////////////////////////////////////////////////////////////
// Our solver class, with customizable real-time output
////////////////////////////////////////////////////////////////

class Dipoles3d : public cudaSolvers::DipolarGPSolver
{
public:

    using DipolarGPSolver::DipolarGPSolver;

    void write_operator_splitting_output(size_t iteration_number,
                                         std::ostream& output_stream) override;

    void set_healing_length(double healing_length){ internal_healing_length = healing_length; };

private:

    double kc;
    double internal_healing_length;

};

void Dipoles3d::write_operator_splitting_output(size_t iteration_number,
                                                std::ostream &output_stream)
{

    evaluate_integrated_occupation_number();
    std::ofstream momentum_cuts;
    momentum_cuts.open("nkx"+std::to_string(iteration_number/write_output_every)+".txt",std::ios::out);
    for(size_t ix = 0; ix < int(nx/2); ++ix)
        momentum_cuts << kx_axis(ix) << " "
                      << integrated_occupation_number[0](ix)/(npoints)
                      << std::endl;
    momentum_cuts.close();

    momentum_cuts.open("nky"+std::to_string(iteration_number/write_output_every)+".txt",std::ios::out);
    for(size_t iy = 0; iy < int(ny/2); ++iy)
        momentum_cuts << ky_axis(iy) << " "
                      << integrated_occupation_number[1](iy)/npoints
                      << std::endl;
    momentum_cuts.close();

    momentum_cuts.open("nkz"+std::to_string(iteration_number/write_output_every)+".txt",std::ios::out);
    for(size_t iz = 0; iz < int(nz/2); ++iz)
        momentum_cuts << kz_axis(iz) << " "
                      << integrated_occupation_number[2](iz)/npoints
                      << std::endl;
    momentum_cuts.close();

    apply_momentum_cutoff(0.1);

    double density_threshold = 0.01*(initial_norm_d[0]/(8*x_axis[nx-1]*y_axis[ny-1]*z_axis[nz-1]));
    calculate_vortex_tangle_length(density_threshold);
    GraphicOutput::DataWriter quasi_condensate_out;
    quasi_condensate_out.set_output_name("vortex_tangle"+std::to_string(iteration_number/write_output_every));
    quasi_condensate_out.write_vtk(x_axis,y_axis,z_axis,vortex_tangle_density,"psi","BINARY");

}

int main() {

    // Read input parameters from file "dipolars3d.prm". This must be placed in the same directory as
    // the executable

    Tools::InputParser ip("example-7.prm");

    ip.read_input_file();

    double xmax = ip.retrieve_double("xmax");
    double ymax = ip.retrieve_double("ymax");
    double zmax = ip.retrieve_double("zmax");

    const int nx = ip.retrieve_int("nx");
    const int ny = ip.retrieve_int("ny");
    const int nz = ip.retrieve_int("ny");

    double scattering_length = ip.retrieve_double("scattering length");
    double dipolar_length    = ip.retrieve_double("dipolar length");
    const double    number_of_particles = ip.retrieve_double("number of particles");
    const double atomic_mass         = ip.retrieve_double("atomic mass");
    double theta = ip.retrieve_double("theta");

    const int    number_of_gradient_descent_steps = ip.retrieve_int("number of gradient descent steps");
    const double alpha                            = ip.retrieve_double("alpha");
    const double beta                             = ip.retrieve_double("beta");

    const int    number_of_real_time_steps = ip.retrieve_int("number of real time steps");
    double time_step                       = ip.retrieve_double("time step");

    const int write_output_every=ip.retrieve_int("write output every");

    const double kcut = ip.retrieve_double("kcut");

    // These two constants are for fixing the units
    const double hbar        = 0.6347*1.E5; // hbar in amu*mum^2/s
    const double bohr_radius = 5.292E-5;    // bohr radius in mum

    // Some constants

    double density = number_of_particles/(8*xmax*ymax*zmax);
    double epsilon_dd = dipolar_length/scattering_length;

    double chemical_potential = 4*PI*scattering_length*bohr_radius*density*(1-epsilon_dd);
    double healing_length = 1.0/std::sqrt(chemical_potential);
    std::cout << "Atomic density = " << density/100.0 << " 10^14 cm^{-3}" << std::endl;
    std::cout << "Chemical potential = " << chemical_potential << std::endl;
    std::cout << "Healing length = " << healing_length << " mum" << std::endl;
    std::cout << "dx = " << 2.*xmax/nx << " mum" << std::endl;
    std::cout << "kmax = " << nx/(4*xmax) << std::endl;
    std::cout << "dk = " <<  1./xmax << std::endl;
    std::cout << "kxi dipo = " << 1./healing_length << std::endl;
    std::cout << "kxi no dipo = " << 1./sqrt(1.-epsilon_dd)*1./healing_length << std::endl;

    scattering_length *= bohr_radius/healing_length;
    dipolar_length    *= bohr_radius/healing_length;

    //////////////////////
    // Create the mesh
    /////////////////////

    Vector<double> x(nx);
    Vector<double> y(ny);
    Vector<double> z(nz);

    xmax = xmax/healing_length;
    ymax = ymax/healing_length;
    zmax = zmax/healing_length;

    double dx = 2.*xmax/nx;
    double dy = 2.*ymax/ny;
    double dz = 2.*zmax/nz;

    for (size_t i = 0; i < nx; ++i) x(i) = -xmax + i*dx;
    for (size_t i = 0; i < ny; ++i) y(i) = -ymax + i*dy;
    for (size_t i = 0; i < nz; ++i) z(i) = -zmax + i*dz;

    double dv = dx*dy*dz;

    //////////////////////////////
    // Create the momentum mesh
    //////////////////////////////

    Vector<double> kx(nx),ky(ny),kz(nz);
    create_mesh_in_Fourier_space(x,y,z,kx,ky,kz);

    ///////////////////////////////////////////////////////////////////
    // Initialize wave function and external potential
    // Wave function is initialized in momentum space populating
    // low momentum modes each with a random phase
    ///////////////////////////////////////////////////////////////////

    Vector<std::complex<double>> psi(nx,ny,nz),psitilde(nx,ny,nz);
    Vector<double> Vext(nx,ny,nz);

    std::default_random_engine generator;

    std::complex<double> ci={0.0,1.0};

    typedef std::chrono::high_resolution_clock clock;
    clock::time_point beginning = clock::now();
    clock::duration d = clock::now() - beginning;
    generator.seed(d.count());
    std::uniform_real_distribution<double> phase_distribution(0,TWOPI);
    std::uniform_real_distribution<double> density_distribution(0,1);

    for (size_t i = 0; i < nx; ++i)
        for (size_t j = 0; j < ny; ++j)
            for (size_t k = 0; k < nz; ++k)
            {
                double random_phase = phase_distribution(generator);
                if (std::abs(kx(i)) <= kcut/TWOPI &&
                    std::abs(ky(j)) <= kcut/TWOPI &&
                    std::abs(kz(k)) <= kcut/TWOPI )
                    psitilde(i,j,k) = sqrt(density)*exp(-ci * random_phase);
            }
    MKLWrappers::DFtCalculator dft(psi,psitilde);
    dft.compute_backward();

    double norm = 0.0;
    for (size_t i = 0; i < psi.size(); ++i) norm += std::norm(psi[i]);
    norm *= dv;
    for (size_t i = 0; i < psi.size(); ++i) psi[i] *= std::sqrt(number_of_particles/norm);

    /////////////////////////////
    // Initialize the solver
    /////////////////////////////

    Dipoles3d gp_solver(x,y,z,psi,Vext,scattering_length,dipolar_length,theta,0,false);

    gp_solver.set_healing_length(healing_length);

    ////////////////////////////////////////////////
    // Run the real time simulation!
    // The output will be written to the shell
    ////////////////////////////////////////////////

    gp_solver.set_tw_initial_conditions(false);
    gp_solver.run_operator_splitting(number_of_real_time_steps,time_step,std::cout,write_output_every);

    return 0;

}