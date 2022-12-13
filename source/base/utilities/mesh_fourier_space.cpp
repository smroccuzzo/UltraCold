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

#include "mesh_fourier_space.hpp"

namespace UltraCold
{

    /**
     * @brief Creation of a mesh in Fourier space for a one-dimensional problem
     *
     * @param x *Vector<double>* representing the x-axis of a cartesian reference frame
     * @param kx *Vector<double>* representing the kx-axis of the corresponding Fourier space
     *
     */

    void create_mesh_in_Fourier_space(Vector<double> &x,
                                      Vector<double> &kx)
    {

        // Check that the dimensions are equal

        int nx = x.size();
        int nkx = kx.size();

        assert(x.order()==1);
        assert(kx.order()==1);
        assert(nx==nkx);

        double dx = abs(x(1) - x(0));
        double xmax = x(nx - 1);

        double kxmax = 1./(2.*dx);
        double dkx = 1./(nx*dx);

        for (size_t i = 0; i < nx/2; ++i)  kx(i) = i * dkx;
        for (size_t i = nx/2; i < nx; ++i) kx(i) = -kxmax + (i-(nx/2)) * dkx;

    }

    /**
     * @brief Creation of a mesh in Fourier space for a two-dimensional problem
     *
     * @param x *Vector<double>* representing the x-axis of a cartesian reference frame
     * @param y *Vector<double>* representing the y-axis of a cartesian reference frame
     * @param kx *Vector<double>* representing the kx-axis of the corresponding Fourier space
     * @param ky *Vector<double>* representing the ky-axis of the corresponding Fourier space
     *
     */

    void create_mesh_in_Fourier_space(Vector<double> &x,
                                      Vector<double> &y,
                                      Vector<double> &kx,
                                      Vector<double> &ky)
    {

        // Check that the dimensions are consistent

        int nx = x.size();
        int nkx = kx.size();

        int ny = y.size();
        int nky = ky.size();

        assert(x.order()==1);
        assert(y.order()==1);
        assert(kx.order()==1);
        assert(ky.order()==1);
        assert(nx==nkx);
        assert(ny==nky);

        // Dimensions are fine, generate the kx and ky axis

        double dx = abs(x(1) - x(0));
        double xmax = x(nx - 1);

        double kxmax = 1. / (2. * dx);
        double dkx = 1. / (nx * dx);

        for (size_t i = 0; i < nx/2; ++i) kx(i) = i * dkx;
        for (size_t i = nx/2; i < nx; ++i) kx(i) = -kxmax + (i-(nx/2)) * dkx;

        double dy = abs(y(1) - y(0));
        double ymax = y(ny - 1);

        double kymax = 1. / (2. * dy);
        double dky = 1. / (ny * dy);

        for (size_t i = 0; i < ny / 2; ++i) ky(i) = i * dky;
        for (size_t i = ny / 2; i < ny; ++i) ky(i) = -kymax + (i-(ny/2)) * dky;


    }

    /**
     * @brief Creation of a mesh in Fourier space for a three-dimensional problem
     *
     * @param x *Vector<double>* representing the x-axis of a cartesian reference frame
     * @param y *Vector<double>* representing the y-axis of a cartesian reference frame
     * @param z *Vector<double>* representing the z-axis of a cartesian reference frame
     * @param kx *Vector<double>* representing the kx-axis of the corresponding Fourier space
     * @param ky *Vector<double>* representing the ky-axis of the corresponding Fourier space
     * @param kz *Vector<double>* representing the kz-axis of the corresponding Fourier space
     *
     */

    void create_mesh_in_Fourier_space(Vector<double> &x, Vector<double> &y, Vector<double> &z,
                                      Vector<double> &kx, Vector<double> &ky, Vector<double> &kz)
    {

        // Check that the dimensions are consistent

        int nx = x.size();
        int nkx = kx.size();

        int ny = y.size();
        int nky = ky.size();

        int nz = z.size();
        int nkz = kz.size();

        assert(x.order()==1);
        assert(y.order()==1);
        assert(z.order()==1);
        assert(kx.order()==1);
        assert(ky.order()==1);
        assert(kz.order()==1);
        assert(nx==nkx);
        assert(ny==nky);
        assert(nz==nkz);

        double dx = abs(x(1) - x(0));
        double xmax = x(nx - 1);

        double kxmax = 1. / (2. * dx);
        double dkx = 1. / (nx * dx);

        for (size_t i = 0; i < nx / 2; ++i) kx(i) = i * dkx;
        for (size_t i = nx / 2; i < nx; ++i) kx(i) = -kxmax + (i - (nx / 2)) * dkx;

        double dy = abs(y(1) - y(0));
        double ymax = y(ny - 1);

        double kymax = 1. / (2. * dy);
        double dky = 1. / (ny * dy);

        for (size_t i = 0; i < ny / 2; ++i) ky(i) = i * dky;
        for (size_t i = ny / 2; i < ny; ++i) ky(i) = -kymax + (i - (ny / 2)) * dky;

        double dz = abs(z(1) - z(0));
        double zmax = z(nz - 1);

        double kzmax = 1. / (2. * dz);
        double dkz = 1. / (nz * dz);

        for (size_t i = 0; i < nz / 2; ++i) kz(i) = i * dkz;
        for (size_t i = nz / 2; i < nz; ++i) kz(i) = -kzmax + (i - (nz / 2)) * dkz;

    }

} // UltraCold