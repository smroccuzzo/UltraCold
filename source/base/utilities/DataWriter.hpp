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

#ifndef ULTRACOLD_GRAPHIC_OUT
#define ULTRACOLD_GRAPHIC_OUT

#include <iostream>
#include <fstream>
#include <cstring>
#include <complex>
#include <assert.h>

#include "Vector.hpp"
#include "InputParser.hpp"

namespace UltraCold
{

    /**
     * @brief Classes and functions to output a data Vector in real space
     * 
     */

    namespace GraphicOutput
    {

        /**
         * @brief A class to output a data Vector in real space
         * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
         *
         * This class allows to write output datafiles in different formats, ready for plots
         * with common visualization tools such as gnuplot, matplotlib, Paraview or Visit. 
         * The currently supported data output formats are:
         *
         *  - .csv
         *  - .vtk
         * 
         * The class is used as follows:
         *
         * \code {.cpp}
         * 
         * DataWriter data_output_writer();
         * data_output_writer.set_output_name("MyOutputFile");
         * data_output_writer.write_csv(axis, Vector); // Or other appropriate member function
         * 
         * \endcode
         * 
         * The first line creates an object of type <code>DataWriter</code> named <code>data_output_writer</code> . \n
         * The second line sets the output file name to "MyOutputFile". \n 
         * The third line provides an example for writing a Vector in a simple .csv file, ready to be plotted for
         * example with gnuplot. \n
         * See the member function documentation for details on the different ways to output your data.
         * 
         */

        class DataWriter
        {

            public:

                DataWriter() = default;
                ~DataWriter() = default;

                void set_output_name(const std::string& output_file_name);
                void set_output_name(const char* output_file_name);

                ////////////////////
                // .csv output    //
                ////////////////////


                // Output a 1D real Vector in a .csv file
                void write_csv(Vector<double>& x_axis,
                               Vector<double>& real_output_vector);

                // Output a 1D complex Vector in a .csv file
                void write_csv(Vector<double>& x_axis,
                               Vector<std::complex<double>>& complex_output_vector);

                // Output a 2D real Vector in a .csv file
                void write_csv(Vector<double>& x_axis,
                               Vector<double>& y_axis,
                               Vector<double>& real_output_vector);

                // Output a 2D complex Vector in a .csv file
                void write_csv(Vector<double>& x_axis,
                               Vector<double>& y_axis,
                               Vector<std::complex<double>>& complex_output_vector);

                // Output a 1D slice of a 2D or 3D real Vector in a .csv file
                void write_slice1d_csv(Vector<double>& ax,
                                       Vector<double>& real_output_vector,
                                       const char* axis_name);

                // Output a 1D slice of a 2D or 3D complex Vector in a .csv file
                void write_slice1d_csv(Vector<double>& axis,
                                       Vector<std::complex<double>>& complex_output_vector,
                                       const char* axis_name);

                // Output a 2D slice of a 3D real Vector in a .csv file
                void write_slice2d_csv(Vector<double>& ax1,
                                       Vector<double>& ax2,
                                       Vector<double>& real_output_vector,
                                       const char* plane_name);

                // Output a 2D slice of a 3D complex Vector in a .csv file
                void write_slice2d_csv(Vector<double>& x_axis,
                                       Vector<double>& y_axis,
                                       Vector<std::complex<double>>& complex_output_vector,
                                       const char* plane_name);

                // Stack 1-dimensional time data and output for a 2-dimensional plot
                void stack1d_csv(double time,
                                 Vector<double>& x_axis,
                                 Vector<double>& real_output_vector);
                void stack1d_csv(double time,
                                 Vector<double>& x_axis,
                                 Vector<std::complex<double>>& complex_output_vector);

                ////////////////////
                // .vtk output    //
                ////////////////////

                // Output a 2D real Vector in a .vtk file
                void write_vtk(Vector<double>& x_axis,
                               Vector<double>& y_axis,
                               Vector<double>& real_output_vector,
                               const char* vector_name,
                               const char* format);

                // Output a 2D complex Vector in a .vtk file
                void write_vtk(Vector<double>& x_axis,
                               Vector<double>& y_axis,
                               Vector<std::complex<double>>& complex_output_vector,
                               const char* vector_name,
                               const char* format);

                // Output a 3D real Vector in a .vtk file
                void write_vtk(Vector<double>& x_axis,
                               Vector<double>& y_axis,
                               Vector<double>& z_axis,
                               Vector<double>& real_output_vector,
                               const char* vector_name,
                               const char* format);

                // Output a 3D complex Vector in a .vtk file
                void write_vtk(Vector<double>& x_axis,
                               Vector<double>& y_axis,
                               Vector<double>& z_axis,
                               Vector<std::complex<double>>& complex_output_vector,
                               const char* vector_name,
                               const char* format);

                // Output a 2D slice of a 3D real Vector in a .vtk file
                void write_slice2d_vtk(Vector<double>& ax1,
                                       Vector<double>& ax2,
                                       Vector<double>& real_output_vector,
                                       const char* vector_name,
                                       const char* plane,
                                       const char* format);

                // Output a 2D slice of a 3D complex Vector in a .vtk file
                void write_slice2d_vtk(Vector<double>& x_axis,
                                       Vector<double>& y_axis,
                                       Vector<std::complex<double>>& complex_output_name,
                                       const char* vector_name,
                                       const char* plane,
                                       const char* format);

            private:

                std::string  output_file_name;
                std::fstream output_stream;

        };

    } // namespace GraphicOutput
} // namespace UltraCold

#endif // ULTRACOLD_GRAPHIC_OUT
