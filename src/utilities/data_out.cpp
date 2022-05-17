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

#include "data_out.hpp"

namespace UltraCold{

    namespace RealSpaceOutput
    {

        /**
         * @brief Set the name for the output data file, input as an std::string
         *
         */

        void DataOut::set_output_name(const std::string& file_name)
        {
            output_file_name = file_name;
        }

        /**
         * @brief Set the name for the output data file, input as simple text
         *
         */

        void DataOut::set_output_name(const char* file_name)
        {
            output_file_name = file_name;
        }

        //////////////////////////////////////////////////
        //                  .csv                        //
        //////////////////////////////////////////////////

        /**
         * @brief Write an output data file in .csv format, for real 1D output
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param v *Vector<double>* the second Vector for the output.
         * This is considered as the data to be plotted on the y-axis of the plot.
         *
         * The output file will have the format
         *
         \verbatim
            x(0),v(0)
            x(1),v(1)
            ...  ...
            x(n),v(n)
        \endverbatim
        *
        * where <code>n</code> is the leading dimension of the two Vectors (which of course must
        * have the same size).  \n
        * To visualize the data written in the file using, for example, <code>gnuplot</code>, open
        * <code>gnuplot</code> in a terminal and type the commands
        * \code
        *  set datafile separator ",";
        *  plot "filename.csv"
        * \endcode
        *
        */

        void DataOut::write_csv(Vector<double>& x,
                                Vector<double>& v)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzv = v.extent(2);


            if(nxx != nxv || nyx != 0 || nzx != 0)
            {

                // The two Vectors do not have the same leading dimension

                std::cout
                    << "\n\n"
                    << "*****************************************************************************\n"
                    << "Error found while trying to write a .csv file for a 1D plot of a real Vector.\n"
                    << "The dimensions of the Vectors provided are not consistent.\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
            else if( nyx != 0 || nyv != 0 || nzx != 0 || nzv != 0)
            {

                // Trying to write a one-dimensional plot for non-one-dimensional Vectors

                std::cout
                    << "\n\n"
                    << "*****************************************************************************\n"
                    << "Error found while trying to write a .csv file for a 1D plot of a real Vector.\n"
                    << "The Vector(s) provided are not one-dimensional.\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }
            else
            {

                // Everything is fine, write the output

                output_stream.open(output_file_name+".csv",std::ios::out);

                    for (size_t i = 0; i < nxx; ++i)
                    {
                        output_stream << x(i) << "," << v(i) << std::endl;
                    }

                output_stream.close();
            }

        }

        /**
         * @brief Write an output data file in .csv format, for complex 1D output
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param v *Vector<std::complex<double>>* the second Vector for the output.
         * This is considered as the data to be plotted on the y-axis of the plot.
         *
         * The output file will have the format
         *
         \verbatim
            x(0),v'(0),v''(0)
            x(1),v'(1),v''(1)
            ...  ...
            x(n),v'(n),v''(n)
        \endverbatim
        *
        * where <code>n</code> is the leading dimension of the two Vectors (which of course must
        * have the same size), v' is the real part of Vector v, and v'' is its imaginary part. \n
        * To visualize the data written in the file using, for example, <code>gnuplot</code>, open
        * <code>gnuplot</code> in a terminal and type the commands
        * \code
        *  set datafile separator ",";
        *  plot "filename.csv" u 1:2 # for the real part, u 1:3 for the imaginary part
        * \endcode
        *
        */

        void DataOut::write_csv(Vector<double>& x,
                                Vector<std::complex<double>>& v)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzv = v.extent(2);

            if(nxx != nxv || nyx != 0 || nzx != 0)
            {

                // The two Vectors do not have the same leading dimension

                std::cout
                    << "\n\n"
                    << "********************************************************************************\n"
                    << "Error found while trying to write a .csv file for a 1D plot of a complex Vector.\n"
                    << "The dimensions of the Vectors provided are not equal.\n"
                    << "Terminating program now...\n"
                    << "********************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
            else if( nyx != 0 || nyv != 0 || nzx != 0 || nzv != 0)
            {

                // Trying to write a one-dimensional plot for non-one-dimensional Vectors

                std::cout
                    << "\n\n"
                    << "********************************************************************************\n"
                    << "Error found while trying to write a .csv file for a 1D plot of a complex Vector.\n"
                    << "The Vector(s) provided are not one-dimensional.\n"
                    << "Terminating program now...\n"
                    << "********************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }
            else
            {

                // Everything is fine, write the output

                output_stream.open(output_file_name+".csv",std::ios::out);

                    for (size_t i = 0; i < nxx; ++i)
                    {
                        output_stream << x(i)        << ","
                                      << v(i).real() << ","
                                      << v(i).imag() << std::endl;
                    }

                output_stream.close();
            }

        }

        /**
         * @brief Write an output data file in .csv format, for real 2D output
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param y *Vector<double>* the second Vector for the output.
         * This is considered as the y-axis for the plot.
         * @param v *Vector<double>* the third Vector for the output.
         * This is considered as the data to be plotted on the z-axis of the plot.
         *
         * The output file will have the format
         *
         \verbatim
            x(0),y(0),v(0,0)
            x(0),y(1),v(0,1)
            ...  ...
            x(0),y(ny-1),v(0,ny-1)

            x(1),y(0),v(1,0)
            x(1),y(1),v(1,1)
            ...  ...  ...
            ...  ...  ...

            x(nx-1),y(ny-1),v(nx-1,ny-1)
        \endverbatim
        *
        * where <code>nx</code> and <code>ny</code> are the dimensions of the Vectors (which of course must
        * be consistent). \n
        * To visualize the data written in the file using, for example, <code>gnuplot</code>, open
        * <code>gnuplot</code> in a terminal and type the commands
        *
        * \code
        *  set datafile separator ",";
        *  set pm3d;
        *  set palette;
        *  unset surface;
        *  set view map;
        *  splot "filename.csv" u 1:2:3
        * \endcode
        *
        */

        void DataOut::write_csv(Vector<double>& x,
                                Vector<double>& y,
                                Vector<double>& v)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxy = y.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyy = y.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzy = y.extent(2);
            int nzv = v.extent(2);

            if(nxx != nxv || nxy != nyv || nyx != 0 || nzx != 0 || nyy != 0 || nzy != 0)
            {

                // The Vectors do not have consistent dimensions

                std::cout
                    << "\n\n"
                    << "*****************************************************************************\n"
                    << "Error found while trying to write a .csv file for a 2D plot of a real Vector.\n"
                    << "The dimensions of the Vectors provided are not consistent.\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
            else if( nyx != 0 || nyy != 0 || nzx != 0 || nzy != 0 || nzv != 0)
            {

                // Trying to write a two-dimensional plot for non-two-dimensional Vector

                std::cout
                    << "\n\n"
                    << "*********************************************************************************\n"
                    << "Error found while trying to write a .csv file for a 2D plot of a real Vector.    \n"
                    << "The Vector provided is not two-dimensional, or the axis are not one-dimensional. \n"
                    << "Terminating program now...\n"
                    << "*********************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }
            else
            {

                // Everything is fine, write the output

                output_stream.open(output_file_name+".csv",std::ios::out);

                    for (size_t i = 0; i < nxx; ++i)
                    {
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            output_stream << x(i) << "," << y(j)
                                                  << "," << v(i,j) << std::endl;
                        }
                        output_stream << std::endl;
                    }

                output_stream.close();
            }

        }

        /**
         * @brief Write an output data file in .csv format, for complex 2D output
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param y *Vector<double>* the second Vector for the output.
         * This is considered as the y-axis for the plot.
         * @param v *Vector<std::complex<double>>* the third Vector for the output.
         * This is considered as the data to be plotted on the z-axis of the plot.
         *
         * The output file will have the format
         *
         \verbatim
            x(0),y(0),v'(0,0),v''(0,0)
            x(0),y(1),v'(0,1),v''(0,1)
            ...  ...
            x(0),y(ny-1),v'(0,ny-1),v''(0,ny-1)

            x(1),y(0),v'(1,0),v''(1,0)
            x(1),y(1),v'(1,1),v''(1,1)
            ...  ...  ...
            ...  ...  ...

            x(nx-1),y(ny-1),v'(nx-1,ny-1),v''(nx-1,ny-1)
        \endverbatim
        *
        * where <code>nx</code> and <code>ny</code> are the dimensions of the Vectors (which of course must
        * be consistent), v' is the real part of the Vector v, and v'' is its imaginary part. \n
        * To visualize the data written in the file using, for example, <code>gnuplot</code>, open
        * <code>gnuplot</code> in a terminal and type the commands
        *
        * \code
        *  set datafile separator ",";
        *  set pm3d;
        *  set palette;
        *  unset surface;
        *  set view map;
        *  splot "filename.csv" u 1:2:3 # for real() part, for the imaginary u 1:2:4
        * \endcode
        *
        */

        void DataOut::write_csv(Vector<double>& x,
                                Vector<double>& y,
                                Vector<std::complex<double>>& v)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxy = y.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyy = y.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzy = y.extent(2);
            int nzv = v.extent(2);

            if(nxx != nxv || nxy != nyv || nyx != 0 || nzx != 0 || nyy != 0 || nzy != 0)
            {

                // The Vectors do not have consistent dimensions

                std::cout
                    << "\n\n"
                    << "********************************************************************************\n"
                    << "Error found while trying to write a .csv file for a 2D plot of a complex Vector.\n"
                    << "The dimensions of the Vectors provided are not consistent.\n"
                    << "Terminating program now...\n"
                    << "********************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
            else if( nyx != 0 || nyy != 0 || nzx != 0 || nzy != 0 || nzv != 0)
            {

                // Trying to write a two-dimensional plot for non-two-dimensional Vector

                std::cout
                    << "\n\n"
                    << "********************************************************************************\n"
                    << "Error found while trying to write a .csv file for a 2D plot of a complex Vector.\n"
                    << "The Vector provided is not two-dimensional, or the axis are not one-dimensional.\n"
                    << "Terminating program now...\n"
                    << "********************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }
            else
            {

                // Everything is fine, write the output

                output_stream.open(output_file_name+".csv",std::ios::out);

                    for (size_t i = 0; i < nxx; ++i)
                    {
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            output_stream << x(i) << "," << y(j)
                                                  << "," << v(i,j).real()
                                                  << "," << v(i,j).imag() << std::endl;
                        }
                        output_stream << std::endl;
                    }

                output_stream.close();
            }

        }

        /**
         * @brief Write an output data file in .csv format, for 1D slice of real 2D or 3D Vector
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param v *Vector<double>* the second Vector for the output.
         * One of its cuts is considered as the y-axis of the plot.
         * @param axis the axis along which the cut is taken. This can either be <code>"x"</code>,
         * <code>"y"</code>, or <code>"z"</code>. For a 2D Vector <code>v</code>, if <code>axis="x"</code> a cut of
         * the 2D Vector <code>v</code> will be taken along the y=0 axis, in the second along the x=0 axis.
         * For a 3D Vector <code>v</code>, if <code>axis="x"</code>, a cut along the intersection of the two planes
         * z=0 and y=0 will be taken, ans similarly if <code>axis="y"</code> or <code>axis="z"</code>.
         *
         * The format of the output data file resembles the one of 1D plots.
         *
         */

        void DataOut::write_slice1d_csv(Vector<double>& x,
                                        Vector<double>& v,
                                        const char* axis)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzv = v.extent(2);

            // Cut of v along the x-axis

            if(strcmp(axis,"x")==0)
            {
                if(nxx != nxv || nyx != 0 || nzx != 0)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for the plot of a 1D slice of a real  "
                           "Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    output_stream.open(output_file_name+".csv",std::ios::out);

                        if(nzv == 0)
                        {

                            // if v is 2D

                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << x(i) << "," << v(i,(int) nyv/2) << std::endl;
                            }

                        }
                        else
                        {
                            // if v is 3D

                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << x(i) << "," << v(i,(int) nyv/2, (int) nzv/2) << std::endl;
                            }
                        }

                    output_stream.close();
                }
            }

            // Cut of v along the y-axis

            else if(strcmp(axis,"y")==0)
            {
                if(nxx != nyv)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for the plot of a 1D slice of a real  "
                           "Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    output_stream.open(output_file_name,std::ios::out);

                        if(nzv == 0)
                        {

                            // if v is 2D

                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << x(i) << "," << v((int) nxv/2,i) << std::endl;
                            }

                        }
                        else
                        {
                            // if v is 3D

                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << x(i) << "," << v((int) nxv/2,i,(int) nzv/2) << std::endl;
                            }
                        }

                    output_stream.close();
                }
            }

            // Cut of v along the z-axis

            else if(strcmp(axis,"z")==0)
            {

                if(nxx != nzv)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for the plot of a 1D slice of a real  "
                           "Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);

                }
                else
                {

                    // Everything is fine, write the output

                    output_stream.open(output_file_name,std::ios::out);

                        if(nzv == 0)
                        {

                            // if v is 2D

                            std::cout
                                << "\n\n"
                                << "*****************************************************************************\n"
                                << "Error found while trying to write a .csv file for the plot of a 1D slice of a"
                                   "  real Vector.\n"
                                << "A 2-dimensional Vector cannot be cut along the z-axis.\n"
                                << "Terminating program now...\n"
                                << "*****************************************************************************\n"
                                << "\n\n"
                                <<
                            std::endl;
                            exit(1);

                        }
                        else
                        {
                            // if v is 3D

                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << x(i) << "," << v((int) nxv/2,(int) nyv/2,i) << std::endl;
                            }
                        }

                    output_stream.close();
                }

            }
            else
            {

                // Writing a wrong axis name

                std::cout
                    << "\n\n"
                    << "*****************************************************************************************\n"
                    << "Error found while trying to write a .csv file for the plot of a 1D slice of a real  "
                       "Vector.\n"
                    << "Name of axis invalid. Please use either \"x\", \"y\" or \"z\".\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }

        }

        /**
         * @brief Write an output data file in .csv format, for 1D slice of complex 2D or 3D Vector
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param v *Vector<double>* the second Vector for the output.
         * One of its cuts is considered as the y-axis of the plot.
         * @param axis the axis along which the cut is taken. This can either be <code>"x"</code>,
         * <code>"y"</code>, or <code>"z"</code>. For a 2D Vector <code>v</code>, if <code>axis="x"</code> a cut of
         * the 2D Vector <code>v</code> will be taken along the y=0 axis, in the second along the x=0 axis.
         * For a 3D Vector <code>v</code>, if <code>axis="x"</code>, a cut along the intersection of the two planes
         * z=0 and y=0 will be taken, ans similarly if <code>axis="y"</code> or <code>axis="z"</code>.
         *
         * The format of the output data file resembles the one of 1D plots.
         *
         */

        void DataOut::write_slice1d_csv(Vector<double>& x,
                                        Vector<std::complex<double>>& v,
                                        const char* axis)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzv = v.extent(2);

            // Cut of v along the x-axis

            if(strcmp(axis,"x")==0)
            {
                if(nxx != nxv || nyx != 0 || nzx != 0)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for the plot of a 1D slice of a complex"
                           "  Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    output_stream.open(output_file_name+".csv",std::ios::out);

                        if(nzv == 0)
                        {

                            // if v is 2D

                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << x(i)                    << ","
                                              << v(i,(int) nyv/2).real() << ","
                                              << v(i,(int) nyv/2).imag() << std::endl;
                            }

                        }
                        else
                        {
                            // if v is 3D

                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << x(i)                                 << ","
                                              << v(i,(int) nyv/2, (int) nzv/2).real() << ","
                                              << v(i,(int) nyv/2, (int) nzv/2).imag() << std::endl;
                            }
                        }

                    output_stream.close();
                }
            }

            // Cut of v along the y-axis

            else if(strcmp(axis,"y")==0)
            {
                if(nxx != nyv)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for the plot of a 1D slice of a complex"
                           "  Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    output_stream.open(output_file_name+".csv",std::ios::out);

                        if(nzv == 0)
                        {

                            // if v is 2D

                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << x(i) << "," << v((int) nxv/2,i).real()
                                                      << "," << v((int) nxv/2,i).imag() << std::endl;
                            }

                        }
                        else
                        {
                            // if v is 3D

                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << x(i) << "," << v((int) nxv/2,i,(int) nzv/2).real()
                                                      << "," << v((int) nxv/2,i,(int) nzv/2).imag() << std::endl;
                            }
                        }

                    output_stream.close();
                }
            }

            // Cut of v along the z-axis

            else if(strcmp(axis,"z")==0)
            {
                if(nxx != nzv)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for the plot of a 1D slice of a complex"
                           "  Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);

                }
                else
                {

                    // Everything is fine, write the output

                    output_stream.open(output_file_name+".csv",std::ios::out);

                        if(nzv == 0)
                        {

                            // if v is 2D

                            std::cout
                                << "\n\n"
                                << "*****************************************************************************\n"
                                << "Error found while trying to write a .csv file for the plot of a 1D slice of a"
                                   "  complex Vector.\n"
                                << "A 2-dimensional Vector cannot be cut along the z-axis.\n"
                                << "Terminating program now...\n"
                                << "*****************************************************************************\n"
                                << "\n\n"
                                <<
                            std::endl;
                            exit(1);

                        }
                        else
                        {
                            // if v is 3D

                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << x(i) << "," << v((int) nxv/2,(int) nyv/2,i).real()
                                                      << "," << v((int) nxv/2,(int) nyv/2,i).imag() << std::endl;
                            }
                        }

                    output_stream.close();
                }

            }
            else
            {

                // Writing a wrong axis name

                std::cout
                    << "\n\n"
                    << "*****************************************************************************************\n"
                    << "Error found while trying to write a .csv file for the plot of a 1D slice of a complex  "
                       "Vector.\n"
                    << "Name of axis invalid. Please use either \"x\", \"y\" or \"z\".\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }

        }

        /**
         * @brief Write an output data file in .csv format, for 2D slice of real 3D Vector
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param y *Vector<double>* the second Vector for the output.
         * This is considered as the y-axis for the plot.
         * @param v *Vector<double>* the third Vector for the output.
         * One of its cuts is considered as the z-axis of the plot.
         * @param plane the plane along which the cut is taken. This can either be <code>"xy"</code>,
         * <code>"yz"</code>, or <code>"xz"</code>. If <code>plane="xy"</code>, a slice along
         * the z=0 plane will be taken, and similarly if <code>axis="yz"</code> or <code>axis="xz"</code>.
         *
         * The format of the output data file resembles the one of 2D plots.
         *
         */

        void DataOut::write_slice2d_csv(Vector<double>& x,
                                        Vector<double>& y,
                                        Vector<double>& v,
                                        const char* plane)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxy = y.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyy = y.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzy = y.extent(2);
            int nzv = v.extent(2);

            // slice of v along the z=0 plane

            if(strcmp(plane,"xy")==0)
            {
                if(nxx != nxv || nxy != nyv || nyx != 0 || nzx != 0 || nyy != 0 || nzy != 0)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for the plot of a 2D slice of a real  "
                           "Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    output_stream.open(output_file_name+".csv",std::ios::out);

                        for (size_t i = 0; i < nxx; ++i)
                        {
                            for(size_t j = 0; j < nxy; ++j)
                            {
                                output_stream << x(i) << ","
                                              << y(j) << ","
                                              << v(i, j, (int) nzv/2) << std::endl;
                            }
                            output_stream << std::endl;
                        }

                    output_stream.close();
                }
            }

            // Slice of v along the x=0 plane

            else if(strcmp(plane,"yz")==0)
            {
                if(nxx != nyv || nxy != nzv)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for the plot of a 2D slice of a real  "
                           "Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    output_stream.open(output_file_name+".csv",std::ios::out);

                        for (size_t i = 0; i < nxx; ++i)
                        {
                            for(size_t j = 0; j < nxy; ++j)
                            {
                                output_stream << x(i) << ","
                                              << y(j) << ","
                                              << v((int) nxv/2, i, j) << std::endl;
                            }
                            output_stream << std::endl;
                        }

                    output_stream.close();
                }
            }

            // Slice of v along the y=0 plane

            else if(strcmp(plane,"xz")==0)
            {
                if(nxx != nxv || nxy != nzv)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "************************************************************************************\n"
                        << "Error found while trying to write a .csv file for the plot of a 2D slice of a real  "
                           "Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    output_stream.open(output_file_name+".csv",std::ios::out);

                        for (size_t i = 0; i < nxx; ++i)
                        {
                            for(size_t j = 0; j < nxy; ++j)
                            {
                                output_stream << x(i) << ","
                                              << y(j) << ","
                                              << v(i, (int) nyv/2, j) << std::endl;
                            }
                            output_stream << std::endl;
                        }

                    output_stream.close();
                }

            }
            else
            {

                // Writing a wrong axis name

                std::cout
                    << "\n\n"
                    << "*****************************************************************************************\n"
                    << "Error found while trying to write a .csv file for the plot of a 1D slice of a real Vector\n"
                    << "Name of plane invalid. Please use either \"xy\", \"yz\" or \"xz\".\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }

        }

        /**
         * @brief Write an output data file in .csv format, for 2D slice of complex 3D Vector
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param y *Vector<double>* the second Vector for the output.
         * This is considered as the y-axis for the plot.
         * @param v *Vector<double>* the third Vector for the output.
         * One of its cuts is considered as the z-axis of the plot.
         * @param plane the plane along which the cut is taken. This can either be <code>"xy"</code>,
         * <code>"yz"</code>, or <code>"xz"</code>. If <code>plane="xy"</code>, a slice along
         * the z=0 plane will be taken, and similarly if <code>axis="yz"</code> or <code>axis="xz"</code>.
         *
         * The format of the output data file resembles the one of 2D plots.
         *
         */

        void DataOut::write_slice2d_csv(Vector<double>& x,
                                        Vector<double>& y,
                                        Vector<std::complex<double>>& v,
                                        const char* plane)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxy = y.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyy = y.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzy = y.extent(2);
            int nzv = v.extent(2);

            // slice of v along the z=0 plane

            if(strcmp(plane,"xy")==0)
            {
                if(nxx != nxv || nxy != nyv || nyx != 0 || nzx != 0 || nyy != 0 || nzy != 0)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for the plot of a 2D slice of a complex"
                           "  Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    output_stream.open(output_file_name+".csv",std::ios::out);

                        for (size_t i = 0; i < nxx; ++i)
                        {
                            for(size_t j = 0; j < nxy; ++j)
                            {
                                output_stream << x(i) << ","
                                              << y(j) << ","
                                              << v(i, j, (int) nzv/2).real() << ","
                                              << v(i, j, (int) nzv/2).imag() << std::endl;
                            }
                            output_stream << std::endl;
                        }
                    output_stream.close();
                }
            }

            // Slice of v along the x=0 plane

            else if(strcmp(plane,"yz")==0)
            {
                if(nxx != nyv || nxy != nzv)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for the plot of a 2D slice of a complex"
                           "  Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    output_stream.open(output_file_name+".csv",std::ios::out);

                        for (size_t i = 0; i < nxx; ++i)
                        {
                            for(size_t j = 0; j < nxy; ++j)
                            {
                                output_stream << x(i) << ","
                                              << y(j) << ","
                                              << v((int) nxv/2, i, j).real() << ","
                                              << v((int) nxv/2, i, j).imag() << std::endl;
                            }
                            output_stream << std::endl;
                        }

                    output_stream.close();
                }
            }

            // Slice of v along the y=0 plane

            else if(strcmp(plane,"xz")==0)
            {
                if(nxx != nxv || nxy != nzv)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for the plot of a 2D slice of a complex"
                           "  Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    output_stream.open(output_file_name+".csv",std::ios::out);

                        for (size_t i = 0; i < nxx; ++i)
                        {
                            for(size_t j = 0; j < nxy; ++j)
                            {
                                output_stream << x(i) << ","
                                              << y(j) << ","
                                              << v(i, (int) nyv/2, j).real() << ","
                                              << v(i, (int) nyv/2, j).imag() << std::endl;
                            }
                            output_stream << std::endl;
                        }
                    output_stream.close();
                }

            }
            else
            {

                // Writing a wrong axis name

                std::cout
                    << "\n\n"
                    << "*****************************************************************************************\n"
                    << "Error found while trying to write a .csv file for the plot of a 1D slice of a complex  "
                       "Vector.\n"
                    << "Name of plane invalid. Please use either \"xy\", \"yz\" or \"xz\".\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }

        }

        /**
         * @brief Stack a one-dimensional, time-varying real vector for a two-dimensional plot.
         * @param time *double* the current instant of time
         * @param x *Vector<double>* space axis
         * @param v *Vector<double>* the vector representing the data for output
         *
         * The output file will be a .csv file with a structure resembling the one of the files generated
         * with, for example, write_csv() for two-dimensional plots of real Vectors. However, here, the first column
         * will contain time, the second space, and third the data to output. This function should be placed inside a
         * time-loop, since it appends new data at the bottom of the file every time it is called.
         *
         */

        void DataOut::stack1d_csv(double time,
                                  Vector<double> & x,
                                  Vector<double> & v)
        {
            // Check the extents are correct

            int nxx = x.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzv = v.extent(2);


            if(nxx != nxv || nyx != 0 || nzx != 0)
            {

                // The two Vectors do not have the same leading dimension

                std::cout
                        << "\n\n"
                        << "*******************************************************************************\n"
                        << "Error found while trying to append a .csv file for a 1D stack of a real Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*******************************************************************************\n"
                        << "\n\n"
                        <<
                        std::endl;
                exit(1);
            }
            else if( nyx != 0 || nyv != 0 || nzx != 0 || nzv != 0)
            {

                // Trying to write a one-dimensional plot for non-one-dimensional Vectors

                std::cout
                        << "\n\n"
                        << "*******************************************************************************\n"
                        << "Error found while trying to append a .csv file for a 1D stack of a real Vector.\n"
                        << "The Vector(s) provided are not one-dimensional.\n"
                        << "Terminating program now...\n"
                        << "*******************************************************************************\n"
                        << "\n\n"
                        <<
                        std::endl;
                exit(1);

            }
            else
            {

                // Everything is fine, write the output

                output_stream.open(output_file_name+".csv",std::ios_base::app);

                output_stream << std::endl;
                for (size_t i = 0; i < nxx; ++i)
                {
                    output_stream << time << "," << x(i) << "," << v(i) << std::endl;
                }

                output_stream.close();
            }

        }

        /**
         * @brief Stack a one-dimensional, time-varying complex vector for a two-dimensional plot.
         * @param time *double* the current instant of time
         * @param x *Vector<double>* space axis
         * @param v *Vector<std::complex<double>>* the vector representing the data for output
         *
         * The output file will be a .csv file with a structure resembling the one of the files generated
         * with, for example, write_csv() for two-dimensional plots of complex Vectors. However, here, the first column
         * will contain time, the second space, and third and fourth the data to output. This function should be placed
         * inside a time-loop, since it appends new data at the bottom of the file every time it is called.
         *
         */

        void DataOut::stack1d_csv(double time,
                                  Vector<double> & x,
                                  Vector<std::complex<double>> & v)
        {
            // Check the extents are correct

            int nxx = x.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzv = v.extent(2);


            if(nxx != nxv || nyx != 0 || nzx != 0)
            {

                // The two Vectors do not have the same leading dimension

                std::cout
                        << "\n\n"
                        << "*******************************************************************************\n"
                        << "Error found while trying to append a .csv file for a 1D stack of a real Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*******************************************************************************\n"
                        << "\n\n"
                        <<
                        std::endl;
                exit(1);
            }
            else if( nyx != 0 || nyv != 0 || nzx != 0 || nzv != 0)
            {

                // Trying to write a one-dimensional plot for non-one-dimensional Vectors

                std::cout
                        << "\n\n"
                        << "*******************************************************************************\n"
                        << "Error found while trying to append a .csv file for a 1D stack of a real Vector.\n"
                        << "The Vector(s) provided are not one-dimensional.\n"
                        << "Terminating program now...\n"
                        << "*******************************************************************************\n"
                        << "\n\n"
                        <<
                        std::endl;
                exit(1);

            }
            else
            {

                // Everything is fine, write the output

                output_stream.open(output_file_name+".csv",std::ios_base::app);

                output_stream << std::endl;
                for (size_t i = 0; i < nxx; ++i)
                {
                    output_stream << time << "," << x(i) << "," << v(i).real() << "," << v(i).imag() << std::endl;
                }

                output_stream.close();
            }

        }

        //////////////////////////////////////////////////
        //                  .vtk                        //
        //////////////////////////////////////////////////

        /**
         * @brief Write an output data file in .vtk format, for real 2D output
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param y *Vector<double>* the second Vector for the output.
         * This is considered as the y-axis for the plot.
         * @param v *Vector<double>* the third Vector for the output.
         * This is considered as the data to be plotted on the z-axis of the plot.
         * @param output_vector_name *char* the name of the output vector
         *
         * The output file will be in standard .vtk format for structured data points.
         * It can be readily visualized using programs like
         * <a href="https://www.paraview.org/"> Paraview </a> or
         * <a href="https://visit-dav.github.io/visit-website/index.html">Visit</a>.
         *
         */

        void DataOut::write_vtk(Vector<double>& x,
                                Vector<double>& y,
                                Vector<double>& v,
                                const char* output_vector_name)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxy = y.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyy = y.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzy = y.extent(2);
            int nzv = v.extent(2);

            if(nxx != nxv || nxy != nyv || nyx != 0 || nzx != 0 || nyy != 0 || nzy != 0)
            {

                // The Vectors do not have consistent dimensions

                std::cout
                    << "\n\n"
                    << "*****************************************************************************\n"
                    << "Error found while trying to write a .vtk file for a 2D plot of a real Vector.\n"
                    << "The dimensions of the Vectors provided are not consistent.\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
            else if( nyx != 0 || nyy != 0 || nzx != 0 || nzy != 0 || nzv != 0)
            {

                // Trying to write a two-dimensional plot for non-two-dimensional Vector

                std::cout
                    << "\n\n"
                    << "*********************************************************************************\n"
                    << "Error found while trying to write a .vtk file for a 2D plot of a real Vector.    \n"
                    << "The Vector provided is not two-dimensional, or the axis are not one-dimensional. \n"
                    << "Terminating program now...\n"
                    << "*********************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }
            else
            {

                // Everything is fine, write the output

                double xmin = x(0);
                double ymin = y(0);
                double dx = std::abs(x(1)-x(0));
                double dy = std::abs(y(1)-y(0));

                output_stream.open(output_file_name+".vtk",std::ios::out);

                    // Header of the .vtk file

                    output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                    output_stream << "# 2D plot of a 2D Real Vector"                    << std::endl;

                    // Data type is text

                    output_stream << "ASCII"                                            << std::endl;

                    // Dataset and mesh information

                    output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                    output_stream << "DIMENSIONS " << nxx  << " " << nxy  << " " << "1" << std::endl;
                    output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                    output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                    // Output the vector

                    output_stream << "POINT_DATA " << nxx*nxy                           << std::endl;
                    output_stream << "SCALARS "    << output_vector_name << " double 1" << std::endl;
                    output_stream << "LOOKUP_TABLE default"                             << std::endl;
                    for (size_t j = 0; j < nxy; ++j)
                    {
                        for (size_t i = 0; i < nxx; ++i)
                        {
                            output_stream << v(i,j) << " ";
                        }
                    }

                output_stream.close();

            }

        }

        /**
         * @brief Write an output data file in .vtk format, for complex 2D output
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param y *Vector<double>* the second Vector for the output.
         * This is considered as the y-axis for the plot.
         * @param v *Vector<std::complex<double>>* the third Vector for the output.
         * This is considered as the data to be plotted on the z-axis of the plot.
         * @param output_vector_name *char* the name of the output vector.
         *
         * The output file will be in standard .vtk format for structured data points.
         * It can be readily visualized using programs like
         * <a href="https://www.paraview.org/"> Paraview </a> or
         * <a href="https://visit-dav.github.io/visit-website/index.html">Visit</a>.
         *
         */

        void DataOut::write_vtk(Vector<double>& x,
                                Vector<double>& y,
                                Vector<std::complex<double>>& v,
                                const char* output_vector_name)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxy = y.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyy = y.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzy = y.extent(2);
            int nzv = v.extent(2);

            if(nxx != nxv || nxy != nyv || nyx != 0 || nzx != 0 || nyy != 0 || nzy != 0)
            {

                // The Vectors do not have consistent dimensions

                std::cout
                    << "\n\n"
                    << "********************************************************************************\n"
                    << "Error found while trying to write a .vtk file for a 2D plot of a complex Vector.\n"
                    << "The dimensions of the Vectors provided are not consistent.\n"
                    << "Terminating program now...\n"
                    << "********************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
            else if( nyx != 0 || nyy != 0 || nzx != 0 || nzy != 0 || nzv != 0)
            {

                // Trying to write a two-dimensional plot for non-two-dimensional Vector

                std::cout
                    << "\n\n"
                    << "********************************************************************************\n"
                    << "Error found while trying to write a .vtk file for a 2D plot of a complex Vector.\n"
                    << "The Vector provided is not two-dimensional, or the axis are not one-dimensional.\n"
                    << "Terminating program now...\n"
                    << "********************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }
            else
            {

                // Everything is fine, write the output

                double xmin = x(0);
                double ymin = y(0);
                double dx = std::abs(x(1)-x(0));
                double dy = std::abs(y(1)-y(0));

                output_stream.open(output_file_name+".vtk",std::ios::out);

                    // Header of the .vtk file

                    output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                    output_stream << "# 2D plot of a 2D Complex Vector"                 << std::endl;

                    // Data type is text

                    output_stream << "ASCII"                                            << std::endl;

                    // Dataset and mesh information

                    output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                    output_stream << "DIMENSIONS " << nxx  << " " << nxy  << " " << "1" << std::endl;
                    output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                    output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                    // Real part of the input vector

                    output_stream << "POINT_DATA "   << nxx*nxy                           << std::endl;
                    output_stream << "SCALARS real_" << output_vector_name << " double 1" << std::endl;
                    output_stream << "LOOKUP_TABLE default"                               << std::endl;
                    for (size_t j = 0; j < nxy; ++j)
                    {
                        for (size_t i = 0; i < nxx; ++i)
                        {
                            output_stream << v(i,j).real() << " ";
                        }
                    }

                    output_stream << std::endl;
                    output_stream << std::endl;

                    // Imaginary part of the input vector

                    output_stream << "SCALARS imag_" << output_vector_name << " double 1" << std::endl;
                    output_stream << "LOOKUP_TABLE default"                               << std::endl;

                    for (size_t j = 0; j < nxy; ++j)
                    {
                        for (size_t i = 0; i < nxx; ++i)
                        {
                            output_stream << v(i,j).imag() << " ";
                        }
                    }
                output_stream.close();

            }

        }

        /**
         * @brief Write an output data file in .vtk format, for real 3D output
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param y *Vector<double>* the second Vector for the output.
         * This is considered as the y-axis for the plot.
         * @param z *Vector<double>* the third Vector for the output.
         * This is considered as the z-axis for the plot.
         * @param v *Vector<double>* the fourth Vector for the output.
         * @param output_vector_name *char* the name of the output vector.
         *
         * The output file will be in standard .vtk format for structured data points.
         * It can be readily visualized using programs like
         * <a href="https://www.paraview.org/"> Paraview </a> or
         * <a href="https://visit-dav.github.io/visit-website/index.html">Visit</a>.
         *
         */

        void DataOut::write_vtk(Vector<double>& x,
                                Vector<double>& y,
                                Vector<double>& z,
                                Vector<double>& v,
                                const char* output_vector_name)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxy = y.extent(0);
            int nxz = z.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyy = y.extent(1);
            int nyz = z.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzy = y.extent(2);
            int nzz = z.extent(2);
            int nzv = v.extent(2);

            if(nxx != nxv || nxy != nyv || nxz != nzv || nyx != 0 ||
               nzx != 0 || nyy != 0 || nzy != 0 || nyz != 0 || nzz != 0)
            {

                // The Vectors do not have consistent dimensions

                std::cout
                    << "\n\n"
                    << "*****************************************************************************\n"
                    << "Error found while trying to write a .vtk file for a 3D plot of a real Vector.\n"
                    << "The dimensions of the Vectors provided are not consistent.\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
            else
            {

                // Everything is fine, write the output

                double xmin = x(0);
                double ymin = y(0);
                double zmin = z(0);
                double dx = std::abs(x(1)-x(0));
                double dy = std::abs(y(1)-y(0));
                double dz = std::abs(z(1)-z(0));

                output_stream.open(output_file_name+".vtk",std::ios::out);

                    // Header of the .vtk file

                    output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                    output_stream << "# 3D plot of a 3D Real Vector"                    << std::endl;

                    // Data type is text

                    output_stream << "ASCII"                                            << std::endl;

                    // Dataset and mesh information

                    output_stream << "DATASET STRUCTURED_POINTS"                         << std::endl;
                    output_stream << "DIMENSIONS " << nxx  << " " << nxy  << " " << nxz  << std::endl;
                    output_stream << "ORIGIN "     << xmin << " " << ymin << " " << zmin << std::endl;
                    output_stream << "SPACING "    << dx   << " " << dy   << " " << dz   << std::endl;

                    // Output the vector

                    output_stream << "POINT_DATA " << nxx*nxy*nxz                       << std::endl;
                    output_stream << "SCALARS "    << output_vector_name << " double 1" << std::endl;
                    output_stream << "LOOKUP_TABLE default"                             << std::endl;
                    for (size_t k = 0; k < nxz; ++k)
                    {
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << v(i,j,k) << " ";
                            }
                        }
                    }

                output_stream.close();

            }

        }

        /**
         * @brief Write an output data file in .vtk format, for complex 3D output
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param y *Vector<double>* the second Vector for the output.
         * This is considered as the y-axis for the plot.
         * @param z *Vector<double>* the third Vector for the output.
         * This is considered as the z-axis for the plot.
         * @param v *Vector<std::complex<double>>* the fourth Vector for the output.
         * @param output_vector_name *char* the name of the output vector.
         *
         * The output file will be in standard .vtk format for structured data points.
         * It can be readily visualized using programs like
         * <a href="https://www.paraview.org/"> Paraview </a> or
         * <a href="https://visit-dav.github.io/visit-website/index.html">Visit</a>.
         *
         */

        void DataOut::write_vtk(Vector<double>& x,
                                Vector<double>& y,
                                Vector<double>& z,
                                Vector<std::complex<double>>& v,
                                const char* output_vector_name)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxy = y.extent(0);
            int nxz = z.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyy = y.extent(1);
            int nyz = z.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzy = y.extent(2);
            int nzz = z.extent(2);
            int nzv = v.extent(2);

            if(nxx != nxv || nxy != nyv || nxz != nzv || nyx != 0 || nzx != 0 ||
               nyy != 0 || nzy != 0 || nyz != 0 || nzz != 0)
            {

                // The Vectors do not have consistent dimensions

                std::cout
                    << "\n\n"
                    << "********************************************************************************\n"
                    << "Error found while trying to write a .vtk file for a 3D plot of a complex Vector.\n"
                    << "The dimensions of the Vectors provided are not consistent.\n"
                    << "Terminating program now...\n"
                    << "********************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
            else
            {

                // Everything is fine, write the output

                double xmin = x(0);
                double ymin = y(0);
                double zmin = z(0);
                double dx = std::abs(x(1)-x(0));
                double dy = std::abs(y(1)-y(0));
                double dz = std::abs(z(1)-z(0));

                output_stream.open(output_file_name+".vtk",std::ios::out);

                    // Header of the .vtk file

                    output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                    output_stream << "# 3D plot of a 3D Real Vector"                    << std::endl;

                    // Data type is text

                    output_stream << "ASCII"                                            << std::endl;

                    // Dataset and mesh information

                    output_stream << "DATASET STRUCTURED_POINTS"                         << std::endl;
                    output_stream << "DIMENSIONS " << nxx  << " " << nxy  << " " << nxz  << std::endl;
                    output_stream << "ORIGIN "     << xmin << " " << ymin << " " << zmin << std::endl;
                    output_stream << "SPACING "    << dx   << " " << dy   << " " << dz   << std::endl;

                    // Output the real() part of the vector

                    output_stream << "POINT_DATA "   << nxx*nxy*nxz                         << std::endl;
                    output_stream << "SCALARS real_" << output_vector_name << " double 1"   << std::endl;
                    output_stream << "LOOKUP_TABLE default"                                 << std::endl;
                    for (size_t k = 0; k < nxz; ++k)
                    {
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << v(i,j,k).real() << " ";
                            }
                        }
                    }

                    // Output the imaginary part of the vector

                    output_stream << "SCALARS imag_" << output_vector_name << " double 1"   << std::endl;
                    output_stream << "LOOKUP_TABLE default"                                 << std::endl;
                    for (size_t k = 0; k < nxz; ++k)
                    {
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << v(i,j,k).imag() << " ";
                            }
                        }
                    }

                output_stream.close();

            }

        }

        /**
         * @brief Write an output data file in .vtk format, for 2D slice of real 3D Vector
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param y *Vector<double>* the second Vector for the output.
         * This is considered as the y-axis for the plot.
         * @param v *Vector<double>* the third Vector for the output.
         * One of its cuts is considered as the z-axis of the plot.
         * @param plane the plane along which the cut is taken. This can either be <code>"xy"</code>,
         * <code>"yz"</code>, or <code>"xz"</code>. If <code>plane="xy"</code>, a slice along
         * the z=0 plane will be taken, and similarly if <code>axis="yz"</code> or <code>axis="xz"</code>.
         * @param output_vector_name *char* the name of the output
         *
         * The output file will be in standard .vtk format for structured data points.
         * It can be readily visualized using programs like
         * <a href="https://www.paraview.org/"> Paraview </a> or
         * <a href="https://visit-dav.github.io/visit-website/index.html">Visit</a>.
         *
         */

        void DataOut::write_slice2d_vtk(Vector<double>& x,
                                        Vector<double>& y,
                                        Vector<double>& v,
                                        const char* plane,
                                        const char* output_vector_name)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxy = y.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyy = y.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzy = y.extent(2);
            int nzv = v.extent(2);

            // slice of v along the z=0 plane

            if(strcmp(plane,"xy")==0)
            {
                if(nxx != nxv || nxy != nyv || nyx != 0 || nzx != 0 || nyy != 0 || nzy != 0)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .vtk file for the plot of a 2D slice of a real "
                           "Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {
                    // Everything is fine, write the output

                    double xmin = x(0);
                    double ymin = y(0);
                    double dx = std::abs(x(1)-x(0));
                    double dy = std::abs(y(1)-y(0));

                    output_stream.open(output_file_name+".vtk",std::ios::out);

                        // Header of the .vtk file

                        output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                        output_stream << "# 2D slice of a 3D Real Vector"                   << std::endl;

                        // Data type is text

                        output_stream << "ASCII"                                            << std::endl;

                        // Dataset and mesh information

                        output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                        output_stream << "DIMENSIONS " << nxx  << " " << nxy  << " " << "1" << std::endl;
                        output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                        output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                        // Output the vector

                        output_stream << "POINT_DATA " << nxx*nxy                           << std::endl;
                        output_stream << "SCALARS "    << output_vector_name << " double 1" << std::endl;
                        output_stream << "LOOKUP_TABLE default"                             << std::endl;
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << v(i,j,(int) nzv/2) << " ";
                            }
                        }

                    output_stream.close();

                }
            }

            // Slice of v along the x=0 plane

            else if(strcmp(plane,"yz")==0)
            {
                if(nxx != nyv || nxy != nzv)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .vtk file for the plot of a 2D slice of a real  "
                           "Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    double xmin = x(0);
                    double ymin = y(0);
                    double dx = std::abs(x(1)-x(0));
                    double dy = std::abs(y(1)-y(0));

                    output_stream.open(output_file_name+".vtk",std::ios::out);

                        // Header of the .vtk file

                        output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                        output_stream << "# 2D slice of a 3D Real Vector"                   << std::endl;

                        // Data type is text

                        output_stream << "ASCII"                                            << std::endl;

                        // Dataset and mesh information

                        output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                        output_stream << "DIMENSIONS " << nxx  << " " << nxy  << " " << "1" << std::endl;
                        output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                        output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                        // Output the vector

                        output_stream << "POINT_DATA " << nxx*nxy                           << std::endl;
                        output_stream << "SCALARS "    << output_vector_name << " double 1" << std::endl;
                        output_stream << "LOOKUP_TABLE default"                             << std::endl;
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << v((int) nxv/2,i,j) << " ";
                            }
                        }

                    output_stream.close();
                }
            }

            // Slice of v along the y=0 plane

            else if(strcmp(plane,"xz")==0)
            {
                if(nxx != nxv || nxy != nzv)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .vtk file for the plot of a 2D slice of a real  "
                           "Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    double xmin = x(0);
                    double ymin = y(0);
                    double dx = std::abs(x(1)-x(0));
                    double dy = std::abs(y(1)-y(0));

                    output_stream.open(output_file_name+".vtk",std::ios::out);

                        // Header of the .vtk file

                        output_stream << "# vtk DataFile Version 3.0"                         << std::endl;
                        output_stream << "# 2D slice of a 3D real() Vector"                   << std::endl;

                        // Data type is text

                        output_stream << "ASCII"                                            << std::endl;

                        // Dataset and mesh information

                        output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                        output_stream << "DIMENSIONS " << nxx  << " " << nxy  << " " << "1" << std::endl;
                        output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                        output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                        // Output the vector

                        output_stream << "POINT_DATA " << nxx*nxy                           << std::endl;
                        output_stream << "SCALARS "    << output_vector_name << " double 1" << std::endl;
                        output_stream << "LOOKUP_TABLE default"                             << std::endl;
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << v(i,(int) nyv/2,j) << " ";
                            }
                        }

                    output_stream.close();
                }

            }
            else
            {

                // Writing a wrong axis name

                std::cout
                    << "\n\n"
                    << "*****************************************************************************************\n"
                    << "Error found while trying to write a .vtk file for the plot of a 2D slice of a real Vector\n"
                    << "Name of plane invalid. Please use either \"xy\", \"yz\" or \"xz\".\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }

        }

        /**
         * @brief Write an output data file in .vtk format, for 2D slice of complex 3D Vector
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param y *Vector<double>* the second Vector for the output.
         * This is considered as the y-axis for the plot.
         * @param v *Vector<std::complex<double>>* the third Vector for the output.
         * One of its cuts is considered as the z-axis of the plot.
         * @param plane the plane along which the cut is taken. This can either be <code>"xy"</code>,
         * <code>"yz"</code>, or <code>"xz"</code>. If <code>plane="xy"</code>, a slice along
         * the z=0 plane will be taken, and similarly if <code>axis="yz"</code> or <code>axis="xz"</code>.
         * @param output_vector_name *char* the name of the output vector.
         *
         * The output file will be in standard .vtk format for structured data points.
         * It can be readily visualized using programs like
         * <a href="https://www.paraview.org/"> Paraview </a> or
         * <a href="https://visit-dav.github.io/visit-website/index.html">Visit</a>.
         *
         */

        void DataOut::write_slice2d_vtk(Vector<double>& x,
                                        Vector<double>& y,
                                        Vector<std::complex<double>>& v,
                                        const char* plane,
                                        const char* output_vector_name)
        {

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxy = y.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyy = y.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzy = y.extent(2);
            int nzv = v.extent(2);

            // slice of v along the z=0 plane

            if(strcmp(plane,"xy")==0)
            {
                if(nxx != nxv || nxy != nyv || nyx != 0 || nzx != 0 || nyy != 0 || nzy != 0)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .vtk file for the plot of a 2D slice of a complex"
                           " Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // Everything is fine, write the output

                    double xmin = x(0);
                    double ymin = y(0);
                    double dx = std::abs(x(1)-x(0));
                    double dy = std::abs(y(1)-y(0));

                    output_stream.open(output_file_name+".vtk",std::ios::out);

                        // Header of the .vtk file

                        output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                        output_stream << "# 2D slice of a 3D Complex Vector"                << std::endl;

                        // Data type is text

                        output_stream << "ASCII"                                            << std::endl;

                        // Dataset and mesh information

                        output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                        output_stream << "DIMENSIONS " << nxx  << " " << nxy  << " " << "1" << std::endl;
                        output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                        output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                        // Output the real() part of the vector

                        output_stream << "POINT_DATA " << nxx*nxy                                << std::endl;
                        output_stream << "SCALARS real_"    << output_vector_name << " double 1" << std::endl;
                        output_stream << "LOOKUP_TABLE default"                                  << std::endl;
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << v(i,j,(int) nzv/2).real() << " ";
                            }
                        }

                        output_stream << std::endl;
                        output_stream << std::endl;

                        // Output the imaginary part of the vector

                        output_stream << "SCALARS imag_"    << output_vector_name << " double 1" << std::endl;
                        output_stream << "LOOKUP_TABLE default"                                  << std::endl;
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << v(i,j,(int) nzv/2).imag() << " ";
                            }
                        }

                    output_stream.close();
                }
            }

            // Slice of v along the x=0 plane

            else if(strcmp(plane,"yz")==0)
            {
                if(nxx != nyv || nxy != nzv)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .vtk file for the plot of a 2D slice of a complex"
                           "  Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {
                    // Everything is fine, write the output

                    double xmin = x(0);
                    double ymin = y(0);
                    double dx = std::abs(x(1)-x(0));
                    double dy = std::abs(y(1)-y(0));

                    output_stream.open(output_file_name+".vtk",std::ios::out);

                        // Header of the .vtk file

                        output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                        output_stream << "# 2D slice of a 3D Complex Vector"                << std::endl;

                        // Data type is text

                        output_stream << "ASCII"                                            << std::endl;

                        // Dataset and mesh information

                        output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                        output_stream << "DIMENSIONS " << nxx  << " " << nxy  << " " << "1" << std::endl;
                        output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                        output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                        // Output the real() part of the vector

                        output_stream << "POINT_DATA " << nxx*nxy                                << std::endl;
                        output_stream << "SCALARS real_"    << output_vector_name << " double 1" << std::endl;
                        output_stream << "LOOKUP_TABLE default"                                  << std::endl;
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << v((int) nxv/2,i,j).real() << " ";
                            }
                        }

                        output_stream << std::endl;
                        output_stream << std::endl;

                        // Output the imaginary part of the vector

                        output_stream << "SCALARS imag_"    << output_vector_name << " double 1" << std::endl;
                        output_stream << "LOOKUP_TABLE default"                                  << std::endl;
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << v((int) nxv/2,i,j).imag() << " ";
                            }
                        }

                    output_stream.close();
                }
            }

            // Slice of v along the y=0 plane

            else if(strcmp(plane,"xz")==0)
            {
                if(nxx != nxv || nxy != nzv)
                {

                    // The Vectors do not have consistent dimensions

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .vtk file for the plot of a 2D slice of a complex"
                           "  Vector.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {
                    // Everything is fine, write the output

                    double xmin = x(0);
                    double ymin = y(0);
                    double dx = std::abs(x(1)-x(0));
                    double dy = std::abs(y(1)-y(0));

                    output_stream.open(output_file_name+".vtk",std::ios::out);

                        // Header of the .vtk file

                        output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                        output_stream << "# 2D slice of a 3D Complex Vector"                << std::endl;

                        // Data type is text

                        output_stream << "ASCII"                                            << std::endl;

                        // Dataset and mesh information

                        output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                        output_stream << "DIMENSIONS " << nxx  << " " << nxy  << " " << "1" << std::endl;
                        output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                        output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                        // Output the real() part of the vector

                        output_stream << "POINT_DATA "      << nxx*nxy                           << std::endl;
                        output_stream << "SCALARS real_"    << output_vector_name << " double 1" << std::endl;
                        output_stream << "LOOKUP_TABLE default"                                  << std::endl;
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << v(i,(int) nyv/2,j).real() << " ";
                            }
                        }

                        output_stream << std::endl;
                        output_stream << std::endl;

                        // Output the imaginary part of the vector

                        output_stream << "SCALARS imag_"    << output_vector_name << " double 1" << std::endl;
                        output_stream << "LOOKUP_TABLE default"                                  << std::endl;
                        for (size_t j = 0; j < nxy; ++j)
                        {
                            for (size_t i = 0; i < nxx; ++i)
                            {
                                output_stream << v(i,(int) nyv/2,j).imag() << " ";
                            }
                        }

                    output_stream.close();
                }

            }
            else
            {

                // Writing a wrong axis name

                std::cout
                    << "\n\n"
                    << "*****************************************************************************************\n"
                    << "Error found while trying to write a .vtk file for the plot of a 2D slice of a complex  "
                       "Vector.\n"
                    << "Name of plane invalid. Please use either \"xy\", \"yz\" or \"xz\".\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }

        }

    } // namespace RealSpaceOutput

    namespace FourierSpaceOutput
    {

        /**
         * @brief Set the name for the output data file, input as an std::string
         *
         */

        void DataOut::set_output_name(const std::string & name)
        {
            output_file_name = name;
        }

        /**
         * @brief Set the name for the output data file, input as simple text
         *
         */

        void DataOut::set_output_name(const char* name)
        {
            output_file_name = name;
        }

        /**
         * @brief Write an output data file in .csv format, for complex 1D output in Fourier space.
         * @param kx *Vector<double>* the first Vector for the output.
         * This is considered as the kx-axis for the plot. For proper output, it must have
         * been generated with a call to MKLWrappers::generate_mesh_in_Fourier_space.
         * @param v *Vector<std::complex<double>>* the second Vector for the output.
         * This is considered as the data to be plotted on the y-axis of the plot.
         *
         * Notice that we have two possible cases: either <code>v</code> comes from a real-to-complex
         * transform, in which case its size will be one-half (plus one) the size of <code>kx</code>,
         * or it comes from a complex-to-complex transform, in which case the the two vectors will have
         * the same size. In each case, the output format will be similar to the one of the corresponding
         * member functions of the class RealSpaceOutput::DataOut, and can be plotted using, for example,
         * <code>gnuplot</code> in the same way. Check the corresponding documentation for
         * RealSpaceOutput::DataOut for more information.
         *
         */

        void DataOut::write_csv(Vector<double>& x,
                                Vector<std::complex<double>>& v)
        {

            // Set the extension for the file

            output_file_name+=".csv";

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzv = v.extent(2);

            if( ( nxx != nxv && (nxx/2+1) != nxv )  || nyx != 0 || nzx != 0)
            {

                // The two Vectors do not have consistent dimensions

                std::cout
                    << "\n\n"
                    << "*****************************************************************************************\n"
                    << "Error found while trying to write a .csv file for a 1D plot of a complex Vector in "
                       "Fourier space.\n"
                    << "The dimensions of the Vectors provided are not consistent.\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
            else if( nyx != 0 || nyv != 0 || nzx != 0 || nzv != 0)
            {

                // Trying to write a one-dimensional plot for non-one-dimensional Vectors

                std::cout
                    << "\n\n"
                    << "*****************************************************************************************\n"
                    << "Error found while trying to write a .csv file for a 1D plot of a complex Vector in "
                       "Fourier space.\n"
                    << "The Vector(s) provided are not one-dimensional.\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }
            else
            {

                // Everything is fine, write the output.

                output_stream.open(output_file_name,std::ios::out);

                    // If the transform was complex-to-complex

                    if(nxx==nxv)
                    {
                        for (size_t i = nxv/2; i < nxv; ++i)
                        {
                            output_stream << x(i)        << ","
                                          << v(i).real() << ","
                                          << v(i).imag() << std::endl;
                        }
                        for (size_t i = 0; i < nxv/2; ++i)
                        {
                            output_stream << x(i)        << ","
                                          << v(i).real() << ","
                                          << v(i).imag() << std::endl;
                        }
                    }

                    // If the transform was real()-to-complex

                    else
                    {

                        for (size_t i = 0; i < nxv; ++i)
                        {
                            output_stream << x(i)        << ","
                                          << v(i).real() << ","
                                          << v(i).imag() << std::endl;
                        }
                    }

                output_stream.close();
            }

        }


        /**
         * @brief Write an output data file in .csv format, for complex 2D output in Fourier space.
         * @param kx *Vector<double>* the first Vector for the output.
         * This is considered as the kx-axis for the plot. For proper output, it must have
         * been generated with a call to MKLWrappers::generate_mesh_in_Fourier_space.
         * @param ky *Vector<double>* the second Vector for the output.
         * This is considered as the kx-axis for the plot. For proper output, it must have
         * been generated with a call to MKLWrappers::generate_mesh_in_Fourier_space.
         * @param v *Vector<double>* the third Vector for the output.
         * This is considered as the data to be plotted on the z-axis of the plot.
         *
         * Notice that we have two possible cases: either <code>v</code> comes from a real-to-complex
         * transform, in which case its size along the last extent will be one-half (plus one) the
         * size of <code>kx</code>, or it comes from a complex-to-complex transform, in which case
         * they will be equal. In each case, the extent along the second direction of <code>v</code>
         * must be the same as the one of <code>ky</code>. The output format will be similar to the
         * one of the corresponding member functions of the class RealSpaceOutput::DataOut, and can
         * be plotted using, for example, <code>gnuplot</code> in the same way. Check the corresponding
         * documentation for RealSpaceOutput::DataOut for more information.
         *
         */

        void DataOut::write_csv(Vector<double>& x,
                                Vector<double>& y,
                                Vector<std::complex<double>>& v)
        {

            // Set the extension for the file

            output_file_name+=".csv";

            // Check the extents are correct

            int nxx = x.extent(0);
            int nxy = y.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyy = y.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzy = y.extent(2);
            int nzv = v.extent(2);

            if( nxx != nxv || ( nxy != nyv && nxy/2+1 != nyv )  )
            {

                // The Vectors do not have consistent dimensions

                std::cout
                    << "\n\n"
                    << "*****************************************************************************************\n"
                    << "Error found while trying to write a .csv file for a 2D plot of a complex Vector in Fourier"
                       " space.\n"
                    << "The dimensions of the Vectors provided are not consistent.\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
            else if( nyx != 0 || nyy != 0 || nzx != 0 || nzy != 0 || nzv != 0)
            {

                // Trying to write a two-dimensional plot for non-two-dimensional Vector

                std::cout
                    << "\n\n"
                    << "*****************************************************************************************\n"
                    << "Error found while trying to write a .csv file for a 2D plot of a complex Vector in "
                       "Fourier  space.\n"
                    << "The Vector provided is not two-dimensional, or the axis are not one-dimensional.\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }
            else
            {

                // Everything is fine, write the output

                output_stream.open(output_file_name,std::ios::out);

                    // If the transform was complex-to-complex

                    if(nxy==nyv)
                    {
                        for (size_t i = nxv/2; i < nxv; ++i)
                        {
                            for (size_t j = nyv/2; j < nyv; ++j)
                            {
                                output_stream << x(i)          << ","
                                              << y(j)          << ","
                                              << v(i,j).real() << ","
                                              << v(i,j).imag() << std::endl;
                            }
                            for (size_t j = 0; j < nyv/2; ++j)
                            {
                                output_stream << x(i)          << ","
                                              << y(j)          << ","
                                              << v(i,j).real() << ","
                                              << v(i,j).imag() << std::endl;
                            }
                            output_stream << std::endl;
                        }
                        for (size_t i = 0; i < nxv/2; ++i)
                        {
                            for (size_t j = nyv/2; j < nyv; ++j)
                            {
                                output_stream << x(i)          << ","
                                              << y(j)          << ","
                                              << v(i,j).real() << ","
                                              << v(i,j).imag() << std::endl;
                            }
                            for (size_t j = 0; j < nyv/2; ++j)
                            {
                                output_stream << x(i)          << ","
                                              << y(j)          << ","
                                              << v(i,j).real() << ","
                                              << v(i,j).imag() << std::endl;
                            }
                            output_stream << std::endl;
                        }
                    }

                    // If the transform was real-to-complex

                    else
                    {
                        for (size_t i = nxv/2; i < nxv; ++i)
                        {
                            for (size_t j = 0; j < nyv-1; ++j)
                            {
                                output_stream << x(i)          << ","
                                              << y(j)          << ","
                                              << v(i,j).real() << ","
                                              << v(i,j).imag() << std::endl;
                            }
                            output_stream << std::endl;
                        }
                        for (size_t i = 0; i < nxv/2; ++i)
                        {
                            for (size_t j = 0; j < nyv-1; ++j)
                            {
                                output_stream << x(i)          << ","
                                              << y(j)          << ","
                                              << v(i,j).real() << ","
                                              << v(i,j).imag() << std::endl;
                            }
                            output_stream << std::endl;
                        }
                    }

                output_stream.close();
            }

        }

        /**
         * @brief Write an output data file in .csv format, for 1D slice of complex 2D or 3D Vector in Fourier
         * space.
         * @param kx *Vector<double>* the first Vector for the output.
         * This is considered as the kx-axis for the plot. For proper output, it must have been generated with a
         * call to MKLWrappers::generate_mesh_in_Fourier_space.
         * @param v *Vector<std::complex<double>>* the second Vector for the output.
         * One of its cuts is considered as the y-axis of the plot.
         * @param axis the axis along which the cut is taken. This can either be <code>"kx"</code>,
         * <code>"ky"</code>, or <code>"kz"</code>. For a 2D Vector <code>v</code>, if <code>axis="kx"</code> a
         * cut of the 2D Vector <code>v</code> will be taken along the ky=0 axis (and vice-versa).
         * For a 3D Vector <code>v</code>, if <code>axis="kx"</code>, a cut along the intersection of the two planes
         * kz=0 and ky=0 will be taken, ans similarly if <code>axis="ky"</code> or <code>axis="kz"</code>.
         *
         * Notice that we have two possible cases: either <code>v</code> comes from a real-to-complex
         * transform, in which case its size along the last extent will be one-half (plus one) the
         * size of <code>ky</code> (or <code>kz</code>), or it comes from a complex-to-complex transform, in
         * which case
         * they will be equal. The output format will be similar to the
         * one of the corresponding member functions of the class RealSpaceOutput::DataOut, and can
         * be plotted using, for example, <code>gnuplot</code> in the same way. Check the corresponding
         * documentation for RealSpaceOutput::DataOut for more information.
         *
         * \warning This function performs only a few range checking, hence be careful in passing
         * consistent vectors as input. If the extents of the Vectors provided are not consistent,
         * a segmentation fault may arise.
         *
         */

        void DataOut::write_slice1d_csv(Vector<double>& x,
                                        Vector<std::complex<double>>& v,
                                        const char* axis)
        {

            // Set the extension for the file

            output_file_name+=".csv";

            // Get the extents

            int nxx = x.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzv = v.extent(2);

            // Correctness of extents is checked later

            output_stream.open(output_file_name,std::ios::out);

            // If it is requested the cut of v along the kx axis

            if(strcmp(axis,"kx")==0)
            {
                if(nxx != nxv)
                {
                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for a 1D slice along the kx axis"
                        << " of a complex Vector in Fourier space.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {
                    if(nzv == 0)
                    {
                        // 1D slice of 2D vector

                        for (size_t i = nxx/2; i < nxx; ++i)
                        {
                            output_stream << x(i)        << ","
                                        << v(i,0).real() << ","
                                        << v(i,0).imag() << "," << std::endl;
                        }
                        for (size_t i = 0; i < nxx/2; ++i)
                        {
                            output_stream << x(i)        << ","
                                          << v(i,0).real() << ","
                                          << v(i,0).imag() << "," << std::endl;
                        }
                    }
                    else
                    {
                        // 1D slice of 3D vector

                        for (size_t i = nxx/2; i < nxx; ++i)
                        {
                            output_stream << x(i)            << ","
                                          << v(i,0,0).real() << ","
                                          << v(i,0,0).imag() << "," << std::endl;
                        }
                        for (size_t i = 0; i < nxx/2; ++i)
                        {
                            output_stream << x(i)            << ","
                                          << v(i,0,0).real() << ","
                                          << v(i,0,0).imag() << "," << std::endl;
                        }

                    }
                }

            }

            // If it is requested the cut of v along the ky axis

            else if(strcmp(axis,"ky")==0)
            {
                if(nxx != nyv && nxx/2+1 != nyv)
                {

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for a 1D slice along the ky axis"
                        << " of a complex Vector in Fourier space.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {
                    if(nzv == 0)
                    {
                        // 1D slice of 2D vector

                        if(nxx == nyv)
                        {

                            // v from a complex-to-complex transform

                            for (size_t i = nxx/2; i < nxx; ++i)
                            {
                                output_stream << x(i)          << ","
                                              << v(0,i).real() << ","
                                              << v(0,i).imag() << "," << std::endl;
                            }
                            for (size_t i = 0; i < nxx/2; ++i)
                            {
                                output_stream << x(i)          << ","
                                              << v(0,i).real() << ","
                                              << v(0,i).imag() << "," << std::endl;
                            }

                        }
                        else
                        {

                            // v from a real()-to-complex transform

                            for (size_t i = 0; i < nyv-1; ++i)
                            {
                                output_stream << x(i)          << ","
                                              << v(0,i).real() << ","
                                              << v(0,i).imag() << "," << std::endl;
                            }
                        }
                    }
                    else
                    {
                        // 1D slice of 3D vector

                        for (size_t i = nxx/2; i < nxx; ++i)
                        {
                            output_stream << x(i)            << ","
                                          << v(0,i,0).real() << ","
                                          << v(0,i,0).imag() << "," << std::endl;
                        }
                        for (size_t i = 0; i < nxx/2; ++i)
                        {
                            output_stream << x(i)          << ","
                                          << v(0,i,0).real() << ","
                                          << v(0,i,0).imag() << "," << std::endl;
                        }
                    }
                }
            }

            // If it is requested the cut of v along the kz axis

            else if(strcmp(axis,"kz")==0)
            {
                if(nxx != nzv && nxx/2+1 != nzv)
                {

                    std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found while trying to write a .csv file for a 1D slice along the kz axis"
                        << " of a complex Vector in Fourier space.\n"
                        << "The dimensions of the Vectors provided are not consistent.\n"
                        << "Terminating program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
                else
                {

                    // 1D slice of 3D vector

                    if(nxx == nzv)
                    {

                        // v from a complex-to-complex transform

                        for (size_t i = nxx/2; i < nxx; ++i)
                        {
                            output_stream << x(i)            << ","
                                          << v(0,0,i).real() << ","
                                          << v(0,0,i).imag() << "," << std::endl;
                        }
                        for (size_t i = 0; i < nxx/2; ++i)
                        {
                            output_stream << x(i)            << ","
                                          << v(0,0,i).real() << ","
                                          << v(0,0,i).imag() << "," << std::endl;
                        }

                    }
                    else
                    {

                        // v from a real-to-complex transform

                        for (size_t i = 0; i < nzv-1; ++i)
                        {
                            output_stream << x(i)            << ","
                                          << v(0,0,i).real() << ","
                                          << v(0,0,i).imag() << "," << std::endl;
                        }
                    }


                }
            }
            else
            {

                // Writing a wrong axis name

                std::cout
                    << "\n\n"
                    << "*****************************************************************************************\n"
                    << "Error found while trying to write a .csv file for the plot of a 1D slice of a complex "
                       "Vector in Fourier space.\n"
                    << "Name of axis invalid. Please use either \"kx\", \"ky\" or \"kz\".\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }

            output_stream.close();

        }

        /**
         * @brief Write an output data file in .csv format, for 2D slice of complex 3D Vector in Fourier space.
         * @param kx *Vector<double>* the first Vector for the output.
         * This is considered as the kx-axis for the plot. For proper output, it must have
         * been generated with a call to MKLWrappers::generate_mesh_in_Fourier_space.
         * @param ky *Vector<double>* the second Vector for the output.
         * This is considered as the ky-axis for the plot. For proper output, it must have
         * been generated with a call to MKLWrappers::generate_mesh_in_Fourier_space.
         * @param v *Vector<double>* the third Vector for the output.
         * One of its cuts is considered as the z-axis of the plot.
         * @param plane the plane along which the cut is taken. This can either be <code>"kxy"</code>,
         * <code>"kyz"</code>, or <code>"kxz"</code>. If <code>plane="kxy"</code>, a slice along
         * the kz=0 plane will be taken, and similarly if <code>axis="kyz"</code> or <code>axis="kxz"</code>.
         *
         * Notice that we have two possible cases: either <code>v</code> comes from a real-to-complex
         * transform, in which case its size along the last extent will be one-half (plus one) the
         * size of <code>ky</code> (or <code>kz</code>), or it comes from a complex-to-complex transform,
         * in which case they will be equal. The output format will be similar to the
         * one of the corresponding member functions of the class RealSpaceOutput::DataOut, and can
         * be plotted using, for example, <code>gnuplot</code> in the same way. Check the corresponding
         * documentation for RealSpaceOutput::DataOut for more information.
         *
         * \warning This function performs only a few range checking, hence be carefull in passing
         * consistent vectors as input. If the extents of the Vectors provided are not consistent,
         * a segmentation fault may arise.
         *
         */

        void DataOut::write_slice2d_csv(Vector<double>& x,
                                        Vector<double>& y,
                                        Vector<std::complex<double>>& v,
                                        const char* plane)
        {

            // Set the extension for the file

            output_file_name+=".csv";

            // Get the extents

            int nxx = x.extent(0);
            int nxy = y.extent(0);
            int nxv = v.extent(0);

            int nyx = x.extent(1);
            int nyy = y.extent(1);
            int nyv = v.extent(1);

            int nzx = x.extent(2);
            int nzy = y.extent(2);
            int nzv = v.extent(2);

            output_stream.open(output_file_name,std::ios::out);

            // If it is requested the slice of v along the kz=0 plane

            if(strcmp(plane,"kxy")==0)
            {

                for (size_t i = nxv/2; i < nxv; ++i)
                {
                    for (size_t j = nyv/2; j < nyv; ++j)
                    {
                        output_stream << x(i)            << ","
                                      << y(j)            << ","
                                      << v(i,j,0).real() << ","
                                      << v(i,j,0).imag() << std::endl;
                    }
                    for (size_t j = 0; j < nyv/2; ++j)
                    {
                        output_stream << x(i)            << ","
                                      << y(j)            << ","
                                      << v(i,j,0).real() << ","
                                      << v(i,j,0).imag() << std::endl;
                    }
                    output_stream << std::endl;
                }
                for (size_t i = 0; i < nxv/2; ++i)
                {
                    for (size_t j = nyv/2; j < nyv; ++j)
                    {
                        output_stream << x(i)            << ","
                                      << y(j)            << ","
                                      << v(i,j,0).real() << ","
                                      << v(i,j,0).imag() << std::endl;
                    }
                    for (size_t j = 0; j < nyv/2; ++j)
                    {
                        output_stream << x(i)            << ","
                                      << y(j)            << ","
                                      << v(i,j,0).real() << ","
                                      << v(i,j,0).imag() << std::endl;
                    }
                    output_stream << std::endl;
                }
            }

            // If it is requested the slice of v along the kyz plane

            else if(strcmp(plane,"kyz")==0)
            {

                // If the tranform was complex-to-complex

                if(nxy == nzv)
                {
                    for (size_t i = nyv/2; i < nyv; ++i)
                    {
                        for (size_t j = nzv/2; j < nzv; ++j)
                        {
                            output_stream << x(i)            << ","
                                          << y(j)            << ","
                                          << v(0,i,j).real() << ","
                                          << v(0,i,j).imag() << std::endl;
                        }
                        for (size_t j = 0; j < nzv/2; ++j)
                        {
                            output_stream << x(i)            << ","
                                          << y(j)            << ","
                                          << v(0,i,j).real() << ","
                                          << v(0,i,j).imag() << std::endl;
                        }
                        output_stream << std::endl;
                    }
                    for (size_t i = 0; i < nyv/2; ++i)
                    {
                        for (size_t j = nzv/2; j < nzv; ++j)
                        {
                            output_stream << x(i)            << ","
                                          << y(j)            << ","
                                          << v(0,i,j).real() << ","
                                          << v(0,i,j).imag() << std::endl;
                        }
                        for (size_t j = 0; j < nzv/2; ++j)
                        {
                            output_stream << x(i)            << ","
                                          << y(j)            << ","
                                          << v(0,i,j).real() << ","
                                          << v(0,i,j).imag() << std::endl;
                        }
                        output_stream << std::endl;
                    }
                }

                // If the transform was real()-to-complex

                else if (nxy/2+1 == nzv)
                {
                    for (size_t i = nyv/2; i < nyv; ++i)
                    {
                        for (size_t j = 0; j < nzv-1; ++j)
                        {
                            output_stream << x(i)            << ","
                                          << y(j)            << ","
                                          << v(0,i,j).real() << ","
                                          << v(0,i,j).imag() << std::endl;
                        }
                        output_stream << std::endl;
                    }
                    for (size_t i = 0; i < nyv/2; ++i)
                    {
                        for (size_t j = 0; j < nzv-1; ++j)
                        {
                            output_stream << x(i)            << ","
                                          << y(j)            << ","
                                          << v(0,i,j).real() << ","
                                          << v(0,i,j).imag() << std::endl;
                        }
                        output_stream << std::endl;
                    }
                }
            }

            // If it is requested the slice of v along the kxz plane

            else if(strcmp(plane,"kxz")==0)
            {

                // If the transform was complex-to-complex

                if(nxy == nzv)
                {
                    for (size_t i = nxv/2; i < nxv; ++i)
                    {
                        for (size_t j = nzv/2; j < nzv; ++j)
                        {
                            output_stream << x(i)            << ","
                                          << y(j)            << ","
                                          << v(i,0,j).real() << ","
                                          << v(i,0,j).imag() << std::endl;
                        }
                        for (size_t j = 0; j < nzv/2; ++j)
                        {
                            output_stream << x(i)            << ","
                                          << y(j)            << ","
                                          << v(i,0,j).real() << ","
                                          << v(i,0,j).imag() << std::endl;
                        }
                        output_stream << std::endl;
                    }
                    for (size_t i = 0; i < nxv/2; ++i)
                    {
                        for (size_t j = nzv/2; j < nzv; ++j)
                        {
                            output_stream << x(i)            << ","
                                          << y(j)            << ","
                                          << v(i,0,j).real() << ","
                                          << v(i,0,j).imag() << std::endl;
                        }
                        for (size_t j = 0; j < nzv/2; ++j)
                        {
                            output_stream << x(i)            << ","
                                          << y(j)            << ","
                                          << v(i,0,j).real() << ","
                                          << v(i,0,j).imag() << std::endl;
                        }
                        output_stream << std::endl;
                    }
                }

                // If the transform was real-to-complex

                else if (nxy/2+1 == nzv)
                {
                    for (size_t i = nxv/2; i < nxv; ++i)
                    {
                        for (size_t j = 0; j < nzv-1; ++j)
                        {
                            output_stream << x(i)            << ","
                                          << y(j)            << ","
                                          << v(i,0,j).real() << ","
                                          << v(i,0,j).imag() << std::endl;
                        }
                        output_stream << std::endl;
                    }
                    for (size_t i = 0; i < nxv/2; ++i)
                    {
                        for (size_t j = 0; j < nzv-1; ++j)
                        {
                            output_stream << x(i)            << ","
                                          << y(j)            << ","
                                          << v(i,0,j).real() << ","
                                          << v(i,0,j).imag() << std::endl;
                        }
                        output_stream << std::endl;
                    }
                }

            }
            else
            {

                // Writing a wrong axis name

                std::cout
                    << "\n\n"
                    << "*****************************************************************************************\n"
                    << "Error found while trying to write a .csv file for the plot of a 1D slice of a complex "
                       "Vector  in Fourier space.\n"
                    << "Name of axis invalid. Please use either \"kx\", \"ky\" or \"kz\".\n"
                    << "Terminating program now...\n"
                    << "*****************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);

            }

            output_stream.close();

        }

    } // namespace FourierSpaceOutput
} // namespace UltraCold