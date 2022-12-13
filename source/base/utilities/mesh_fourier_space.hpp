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

#ifndef ULTRACOLD_MESH_FOURIER_SPACE
#define ULTRACOLD_MESH_FOURIER_SPACE

#include <iostream>
#include <assert.h>

#include "Vector.hpp"

namespace UltraCold
{

    /**
     * @brief Function to create a mesh in Fourier space
     * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
     *
     * This function allows to generate a mesh in Fourier space, starting from a corresponding mesh in real,
     * cartesian space. It is also useful to plot the Fourier transform of a Vector. \n
     * The function is used as follows
     *
     * \code {.cpp}
     *
     * create_mesh_in_Fourier_space(x,kx);          // for a 1D problem
     * create_mesh_in_Fourier_space(x,y,kx,ky);     // for a 2D problem
     * create_mesh_in_Fourier_space(x,y,z,kx,ky,kz) // for a 3D problem
     *
     * \endcode
     *
     * This will take the three Vectors x,y and z, representing a mesh in real space, and use them to generate the
     * corresponding mesh in Fourier space. \n
     *
     */

    void create_mesh_in_Fourier_space(Vector<double>& x,
                                      Vector<double>& kx);

    void create_mesh_in_Fourier_space(Vector<double>& x,
                                      Vector<double>& y,
                                      Vector<double>& kx,
                                      Vector<double>& ky);

    void create_mesh_in_Fourier_space(Vector<double>& x,
                                      Vector<double>& y,
                                      Vector<double>& z,
                                      Vector<double>& kx,
                                      Vector<double>& ky,
                                      Vector<double>& kz);

} // UltraCold

#endif //ULTRACOLD_MESH_FOURIER_SPACE
