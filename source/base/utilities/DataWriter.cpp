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

#include "DataWriter.hpp"

namespace UltraCold
{

    namespace GraphicOutput
    {

        /**
         * @brief Set the name for the output data file, input as an std::string
         *
         */

        void DataWriter::set_output_name(const std::string& file_name)
        {
            output_file_name = file_name;
        }

        /**
         * @brief Set the name for the output data file, input as simple text
         *
         */

        void DataWriter::set_output_name(const char* file_name)
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

        void DataWriter::write_csv(Vector<double>& x,
                                   Vector<double>& v)
        {
            int nx = x.extent(0);
            assert(x.order() == 1);
            assert(v.order() == 1);
            assert(v.extent(0) == nx);
            output_stream.open(output_file_name+".csv",std::ios::out);
                for (size_t i = 0; i < nx; ++i)
                    output_stream << x(i) << "," << v(i) << std::endl;
            output_stream.close();

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

        void DataWriter::write_csv(Vector<double>& x,
                                   Vector<std::complex<double>>& v)
        {

            // Check the extents are correct

            int nx = x.extent(0);
            assert(x.order() == 1);
            assert(v.order() == 1);
            assert(v.extent(0) == nx);

            // Everything is fine, write the output

            output_stream.open(output_file_name+".csv",std::ios::out);
                for (size_t i = 0; i < nx; ++i)
                    output_stream << x(i)        << ","
                                  << v(i).real() << ","
                                  << v(i).imag() << std::endl;
            output_stream.close();


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

        void DataWriter::write_csv(Vector<double>& x,
                                   Vector<double>& y,
                                   Vector<double>& v)
        {

            // Check the extents are correct

            int nx = x.extent(0);
            int ny = y.extent(0);

            assert(x.order()==1);
            assert(y.order()==1);
            assert(v.extent(0)==nx);
            assert(v.extent(1)==ny);

            // Everything is fine, write the output

            output_stream.open(output_file_name+".csv",std::ios::out);

                for (size_t i = 0; i < nx; ++i)
                {
                    for (size_t j = 0; j < ny; ++j)
                    {
                        output_stream << x(i) << "," << y(j)
                                              << "," << v(i,j) << std::endl;
                    }
                    output_stream << std::endl;
                }

            output_stream.close();


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

        void DataWriter::write_csv(Vector<double>& x,
                                   Vector<double>& y,
                                   Vector<std::complex<double>>& v)
        {

            // Check the extents are correct

            int nx = x.extent(0);
            int ny = y.extent(0);

            assert(x.order()==1);
            assert(y.order()==1);
            assert(v.extent(0)==nx);
            assert(v.extent(1)==ny);

            // Everything is fine, write the output

            output_stream.open(output_file_name+".csv",std::ios::out);

                for (size_t i = 0; i < nx; ++i)
                {
                    for (size_t j = 0; j < ny; ++j)
                    {
                        output_stream << x(i) << "," << y(j)
                                              << "," << v(i,j).real()
                                              << "," << v(i,j).imag() << std::endl;
                    }
                    output_stream << std::endl;
                }

            output_stream.close();

        }

        /**
         * @brief Write an output data file in .csv format, for 1D slice of real 2D or 3D Vector
         * @param ax *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param v *Vector<double>* the second Vector for the output.
         * One of its cuts is considered as the y-axis of the plot.
         * @param axis the axis along which the cut is taken. This can either be <code>"ax"</code>,
         * <code>"y"</code>, or <code>"z"</code>. For a 2D Vector <code>v</code>, if <code>axis="ax"</code> a cut of
         * the 2D Vector <code>v</code> will be taken along the y=0 axis, in the second along the ax=0 axis.
         * For a 3D Vector <code>v</code>, if <code>axis="ax"</code>, a cut along the intersection of the two planes
         * z=0 and y=0 will be taken, ans similarly if <code>axis="y"</code> or <code>axis="z"</code>.
         *
         * The format of the output data file resembles the one of 1D plots.
         *
         */

        void DataWriter::write_slice1d_csv(Vector<double>& ax,
                                           Vector<double>& v,
                                           const char* axis)
        {

            // Some check on the shape of the vectors

            assert(ax.order() == 1);
            assert(v.order() == 2 || v.order() == 3);
            assert(strcmp(axis,"x") || strcmp(axis,"y") || strcmp(axis,"z"));

            // Open the output stream

            output_stream.open(output_file_name+".csv",std::ios::out);

            // Cut of v along the x-axis

            if(strcmp(axis,"x")==0)
            {
                int nx = ax.extent(0);
                assert(v.extent(0)==nx);
                if(v.order()==2)
                {
                    int ny = v.extent(1);
                    for (size_t i = 0; i < nx; ++i)
                        output_stream << ax(i) << "," << v(i,(int)ny/2) << std::endl;
                }
                if(v.order() == 3)
                {
                    int ny=v.extent(1);
                    int nz=v.extent(2);
                    for (size_t i = 0; i < nx; ++i)
                        output_stream << ax(i) << "," << v(i,(int)ny/2,(int)nz/2) << std::endl;
                }
                output_stream.close();
            }

            // Cut of v along the y-axis

            else if(strcmp(axis,"y")==0)
            {
                int ny=ax.extent(0);
                assert(v.extent(1)==ny);
                output_stream.open(output_file_name,std::ios::out);
                    if(v.order()==2)
                    {
                        int nx=v.extent(0);
                        for (size_t i = 0; i < ny; ++i)
                            output_stream << ax(i) << "," << v((int)nx/2,i) << std::endl;
                    }
                    if(v.order()==3)
                    {
                        int nx=v.extent(0);
                        int nz=v.extent(2);
                        for (size_t i = 0; i < ny; ++i)
                            output_stream << ax(i) << "," << v((int)nx/2,i,(int)nz/2) << std::endl;
                    }

                output_stream.close();
            }

            // Cut of v along the z-axis

            else if(strcmp(axis,"z")==0)
            {
                int nz=ax.extent(0);
                assert(v.order()==3);
                assert(v.extent(2)==nz);
                int nx=v.extent(0);
                int ny=v.extent(1);
                output_stream.open(output_file_name,std::ios::out);
                for (size_t i = 0; i < nz; ++i)
                    output_stream << ax(i) << "," << v((int)nx/2,(int)ny/2,i) << std::endl;
                output_stream.close();
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

        void DataWriter::write_slice1d_csv(Vector<double>& ax,
                                           Vector<std::complex<double>>& v,
                                           const char* axis)
        {

            // Some check on the shape of the vectors

            assert(ax.order() == 1);
            assert(v.order() == 2 || v.order() == 3);
            assert(strcmp(axis,"x") || strcmp(axis,"y") || strcmp(axis,"z"));

            // Open the output stream

            output_stream.open(output_file_name+".csv",std::ios::out);

            // Cut of v along the x-axis

            if(strcmp(axis,"x")==0)
            {
                int nx = ax.extent(0);
                assert(v.extent(0)==nx);
                if(v.order()==2)
                {
                    int ny = v.extent(1);
                    for (size_t i = 0; i < nx; ++i)
                        output_stream << ax(i) << ","
                                      << v(i,(int)ny/2).real() << ","
                                      << v(i,(int)ny/2).imag() << std::endl;
                }
                if(v.order() == 3)
                {
                    int ny=v.extent(1);
                    int nz=v.extent(2);
                    for (size_t i = 0; i < nx; ++i)
                        output_stream << ax(i) << ","
                                      << v(i,(int)ny/2,(int)nz/2).real() << ","
                                      << v(i,(int)ny/2,(int)nz/2).imag() << std::endl;
                }
                output_stream.close();
            }

            // Cut of v along the y-axis

            else if(strcmp(axis,"y")==0)
            {
                int ny=ax.extent(0);
                assert(v.extent(1)==ny);
                output_stream.open(output_file_name,std::ios::out);
                if(v.order()==2)
                {
                    int nx=v.extent(0);
                    for (size_t i = 0; i < ny; ++i)
                        output_stream << ax(i) << ","
                                      << v((int)nx/2,i).real() << ","
                                      << v((int)nx/2,i).imag() << std::endl;
                }
                if(v.order()==3)
                {
                    int nx=v.extent(0);
                    int nz=v.extent(2);
                    for (size_t i = 0; i < ny; ++i)
                        output_stream << ax(i) << ","
                                      << v((int)nx/2,i,(int)nz/2).real() << ","
                                      << v((int)nx/2,i,(int)nz/2).imag() << std::endl;
                }

                output_stream.close();
            }

            // Cut of v along the z-axis

            else if(strcmp(axis,"z")==0)
            {
                int nz=ax.extent(0);
                assert(v.order()==3);
                assert(v.extent(2)==nz);
                int nx=v.extent(0);
                int ny=v.extent(1);
                output_stream.open(output_file_name,std::ios::out);
                for (size_t i = 0; i < nz; ++i)
                    output_stream << ax(i) << ","
                                  << v((int)nx/2,(int)ny/2,i).real() << ","
                                  << v((int)nx/2,(int)ny/2,i).imag() << std::endl;
                output_stream.close();
            }
        }

        /**
         * @brief Write an output data file in .csv format, for 2D slice of real 3D Vector
         * @param ax1 *Vector<double>* the first Vector for the output.
         * This is considered as the ax1-axis for the plot.
         * @param ax2 *Vector<double>* the second Vector for the output.
         * This is considered as the ax2-axis for the plot.
         * @param v *Vector<double>* the third Vector for the output.
         * One of its cuts is considered as the z-axis of the plot.
         * @param plane the plane along which the cut is taken. This can either be <code>"xy"</code>,
         * <code>"yz"</code>, or <code>"xz"</code>. If <code>plane="xy"</code>, a slice along
         * the z=0 plane will be taken, and similarly if <code>axis="yz"</code> or <code>axis="xz"</code>.
         *
         * The format of the output data file resembles the one of 2D plots.
         *
         */

        void DataWriter::write_slice2d_csv(Vector<double>& ax1,
                                           Vector<double>& ax2,
                                           Vector<double>& v,
                                           const char* plane)
        {

            // Check the extents are correct
            assert(ax1.order()==0);
            assert(ax2.order()==0);
            assert(v.order()==3);
            assert(strcmp(plane,"xy") || strcmp(plane,"xz") || strcmp(plane,"yz") );

            // slice of v along the z=0 plane

            if(strcmp(plane,"xy")==0)
            {
                int nx=ax1.extent(0);
                int ny=ax2.extent(0);
                int nz=v.extent(2);
                assert(v.extent(0)==nx);
                assert(v.extent(1)==ny);

                output_stream.open(output_file_name+".csv",std::ios::out);
                    for (size_t i = 0; i < nx; ++i)
                    {
                        for(size_t j = 0; j < ny; ++j)
                        {
                            output_stream << ax1(i) << ","
                                          << ax2(j) << ","
                                          << v(i,j,(int)nz/2) << std::endl;
                        }
                        output_stream << std::endl;
                    }
                output_stream.close();
            }

            // Slice of v along the x=0 plane

            else if(strcmp(plane,"yz")==0)
            {
                int ny=ax1.extent(0);
                int nz=ax2.extent(0);
                int nx=v.extent(0);
                assert(v.extent(1)==ny);
                assert(v.extent(2)==nz);

                output_stream.open(output_file_name+".csv",std::ios::out);

                    for (size_t i = 0; i < ny; ++i)
                    {
                        for(size_t j = 0; j < nz; ++j)
                        {
                            output_stream << ax1(i) << ","
                                          << ax2(j) << ","
                                          << v((int)nx/2, i, j) << std::endl;
                        }
                        output_stream << std::endl;
                    }
                output_stream.close();
            }

            // Slice of v along the ax2=0 plane

            else if(strcmp(plane,"xz")==0)
            {
                int nx=ax1.extent(0);
                int nz=ax2.extent(0);
                int ny=v.extent(1);
                assert(v.extent(0)==nx);
                assert(v.extent(2)==nz);

                output_stream.open(output_file_name+".csv",std::ios::out);
                    for (size_t i = 0; i < nx; ++i)
                    {
                        for(size_t j = 0; j < nz; ++j)
                        {
                            output_stream << ax1(i) << ","
                                          << ax2(j) << ","
                                          << v(i,(int)ny/2,j) << std::endl;
                        }
                        output_stream << std::endl;
                    }
                output_stream.close();
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

        void DataWriter::write_slice2d_csv(Vector<double>& ax1,
                                           Vector<double>& ax2,
                                           Vector<std::complex<double>>& v,
                                           const char* plane)
        {

            // Check the extents are correct
            assert(ax1.order()==0);
            assert(ax2.order()==0);
            assert(v.order()==3);
            assert(strcmp(plane,"xy") || strcmp(plane,"xz") || strcmp(plane,"yz") );

            // slice of v along the z=0 plane

            if(strcmp(plane,"xy")==0)
            {
                int nx=ax1.extent(0);
                int ny=ax2.extent(0);
                int nz=v.extent(2);
                assert(v.extent(0)==nx);
                assert(v.extent(1)==ny);

                output_stream.open(output_file_name+".csv",std::ios::out);
                for (size_t i = 0; i < nx; ++i)
                {
                    for(size_t j = 0; j < ny; ++j)
                    {
                        output_stream << ax1(i) << ","
                                      << ax2(j) << ","
                                      << v(i,j,(int)nz/2).real() << ","
                                      << v(i,j,(int)nz/2).imag() << std::endl;
                    }
                    output_stream << std::endl;
                }
                output_stream.close();
            }

                // Slice of v along the x=0 plane

            else if(strcmp(plane,"yz")==0)
            {
                int ny=ax1.extent(0);
                int nz=ax2.extent(0);
                int nx=v.extent(0);
                assert(v.extent(1)==ny);
                assert(v.extent(2)==nz);

                output_stream.open(output_file_name+".csv",std::ios::out);

                for (size_t i = 0; i < ny; ++i)
                {
                    for(size_t j = 0; j < nz; ++j)
                    {
                        output_stream << ax1(i) << ","
                                      << ax2(j) << ","
                                      << v((int)nx/2,i,j).real() << ","
                                      << v((int)nx/2,i,j).imag() << std::endl;
                    }
                    output_stream << std::endl;
                }
                output_stream.close();
            }

            // Slice of v along the ax2=0 plane

            else if(strcmp(plane,"xz")==0)
            {
                int nx=ax1.extent(0);
                int nz=ax2.extent(0);
                int ny=v.extent(1);
                assert(v.extent(0)==nx);
                assert(v.extent(2)==nz);

                output_stream.open(output_file_name+".csv",std::ios::out);
                for (size_t i = 0; i < nx; ++i)
                {
                    for(size_t j = 0; j < nz; ++j)
                    {
                        output_stream << ax1(i) << ","
                                      << ax2(j) << ","
                                      << v(i,(int)ny/2,j).real() << ","
                                      << v(i,(int)ny/2,j).imag() << std::endl;
                    }
                    output_stream << std::endl;
                }
                output_stream.close();
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

        void DataWriter::stack1d_csv(double time,
                                     Vector<double> & x,
                                     Vector<double> & v)
        {
            // Check the extents are correct

            assert(x.order()==1);
            assert(v.order()==1);
            int nx=x.extent(0);
            assert(v.extent(0) == nx);

            // Everything is fine, write the output

            output_stream.open(output_file_name+".csv",std::ios_base::app);

            output_stream << std::endl;
            for (size_t i = 0; i < nx; ++i)
                output_stream << time << "," << x(i) << "," << v(i) << std::endl;
            output_stream.close();

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

        void DataWriter::stack1d_csv(double time,
                                     Vector<double> & x,
                                     Vector<std::complex<double>> & v)
        {
            // Check the extents are correct

            assert(x.order()==1);
            assert(v.order()==1);
            int nx=x.extent(0);
            assert(v.extent(0) == nx);

            // Everything is fine, write the output

            output_stream.open(output_file_name+".csv",std::ios_base::app);
            output_stream << std::endl;
            for (size_t i = 0; i < nx; ++i)
                output_stream << time << "," << x(i) << "," << v(i).real() << "," << v(i).imag() << std::endl;
            output_stream.close();
        }

        //////////////////////////////////////////////////
        //                  .vtk                        //
        //////////////////////////////////////////////////

        /**
         * @brief Utility that swaps byte ordering from Little Endian to Big Endian. Necessary for writing vtk binary
         * files.
         * @param var The variable whose bytes are going to be swapped
         *
         * */

        template <typename T>
        void SwapEnd(T& var)
        {
            char* varArray = reinterpret_cast<char*>(&var);
            for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
                std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
        }

        /**
         *
         * @brief Write an output data file in .vtk format, for real 2D output
         * @param x *Vector<double>* the first Vector for the output.
         * This is considered as the x-axis for the plot.
         * @param y *Vector<double>* the second Vector for the output.
         * This is considered as the y-axis for the plot.
         * @param v *Vector<double>* the third Vector for the output.
         * This is considered as the data to be plotted on the z-axis of the plot.
         * @param output_vector_name *char* the name of the output vector
         * @param format *char* format of the vtk file, either "ASCII" of "BINARY". Notice that BINARY files occupy less
         * memory and are faster to read/write. However, they may be harder to debug and/or to access via python scripts.
         *
         * The output file will be in standard .vtk format for structured data points.
         * It can be readily visualized using programs like
         * <a href="https://www.paraview.org/"> Paraview </a> or
         * <a href="https://visit-dav.github.io/visit-website/index.html">Visit</a>.
         * \note The data is stored in binary format, and in Big Endian byte ordering.
         *
         */

        void DataWriter::write_vtk(Vector<double>& x,
                                   Vector<double>& y,
                                   Vector<double>& v,
                                   const char* output_vector_name,
                                   const char* format)
        {

            assert(x.order()==1);
            assert(y.order()==1);
            assert(v.order()==2);
            int nx=x.extent(0);
            int ny=y.extent(0);
            assert(v.extent(0)==nx);
            assert(v.extent(1)==ny);
            assert(strcmp(format,"ASCII") || strcmp(format,"BINARY"));

            // Everything is fine, write the output

            double xmin = x(0);
            double ymin = y(0);
            double dx = std::abs(x(1)-x(0));
            double dy = std::abs(y(1)-y(0));

            output_stream.open(output_file_name+".vtk",std::ios_base::out | std::ios_base::binary);

                // Header of the .vtk file

                output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                output_stream << "# 2D plot of a 2D Real Vector"                    << std::endl;

                // Data type, either ASCII or binary

                if(strcmp(format,"ASCII")==0)
                    output_stream << "ASCII"                                            << std::endl;
                else if (strcmp(format,"BINARY")==0)
                    output_stream << "BINARY"                                           << std::endl;


               // Dataset and mesh information

                output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                output_stream << "DIMENSIONS " << nx   << " " << ny   << " " << "1" << std::endl;
                output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                // Output the vector

                output_stream << "POINT_DATA " << nx*ny                             << std::endl;
                output_stream << "SCALARS "    << output_vector_name << " double 1" << std::endl;
                output_stream << "LOOKUP_TABLE default"                             << std::endl;
                if(strcmp(format,"BINARY")==0)
                {
                    for (size_t j = 0; j < ny; ++j)
                        for (size_t i = 0; i < nx; ++i)
                        {
                            double value = v(i,j);
                            SwapEnd(value);
                            output_stream.write(reinterpret_cast<char*>(&value),sizeof(double));
                        }
                }
                else if(strcmp(format,"ASCII") == 0)
                {
                    for (size_t j = 0; j < ny; ++j)
                        for (size_t i = 0; i < nx; ++i)
                            output_stream << v(i,j) << " ";
                }


            output_stream.close();

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
         * @param format *char* format of the vtk file, either "ASCII" of "BINARY". Notice that BINARY files occupy less
         * memory and are faster to read/write. However, they may be harder to debug and/or to access via python scripts.
         *
         * The output file will be in standard .vtk format for structured data points.
         * It can be readily visualized using programs like
         * <a href="https://www.paraview.org/"> Paraview </a> or
         * <a href="https://visit-dav.github.io/visit-website/index.html">Visit</a>.
         * \note The data is stored in binary format, and in Big Endian byte ordering.
         */

        void DataWriter::write_vtk(Vector<double>& x,
                                   Vector<double>& y,
                                   Vector<std::complex<double>>& v,
                                   const char* output_vector_name,
                                   const char* format)
        {

            assert(x.order()==1);
            assert(y.order()==1);
            assert(v.order()==2);
            int nx=x.extent(0);
            int ny=y.extent(0);
            assert(v.extent(0)==nx);
            assert(v.extent(1)==ny);
            assert(strcmp(format,"ASCII") || strcmp(format,"BINARY"));

            // Everything is fine, write the output

            double xmin = x(0);
            double ymin = y(0);
            double dx = std::abs(x(1)-x(0));
            double dy = std::abs(y(1)-y(0));

            output_stream.open(output_file_name+".vtk",std::ios_base::out | std::ios_base::binary);

                // Header of the .vtk file

                output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                output_stream << "# 2D plot of a 2D Complex Vector"                 << std::endl;

                // Data type

                if(strcmp(format,"ASCII")==0)
                    output_stream << "ASCII"                                            << std::endl;
                else if (strcmp(format,"BINARY")==0)
                    output_stream << "BINARY"                                           << std::endl;

                // Dataset and mesh information

                output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                output_stream << "DIMENSIONS " << nx   << " " << ny   << " " << "1" << std::endl;
                output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                // Real part of the input vector

                output_stream << "POINT_DATA "   << nx*ny                             << std::endl;
                output_stream << "SCALARS real_" << output_vector_name << " double 1" << std::endl;
                output_stream << "LOOKUP_TABLE default"                               << std::endl;

                if(strcmp(format,"BINARY")==0)
                {
                    for (size_t j = 0; j < ny; ++j) {
                        for (size_t i = 0; i < nx; ++i)
                        {
                            double value = v(i, j).real();
                            SwapEnd(value);
                            output_stream.write(reinterpret_cast<char *>(&value), sizeof(double));
                        }
                    }
                }
                else if(strcmp(format,"ASCII")==0)
                {
                    for (size_t j = 0; j < ny; ++j)
                        for (size_t i = 0; i < nx; ++i)
                            output_stream << v(i, j).real() << " ";
                }

                output_stream << std::endl;

                // Imaginary part of the input vector

                output_stream << "SCALARS imag_" << output_vector_name << " double 1" << std::endl;
                output_stream << "LOOKUP_TABLE default"                               << std::endl;

                if(strcmp(format,"BINARY")==0)
                {
                    for (size_t j = 0; j < ny; ++j) {
                        for (size_t i = 0; i < nx; ++i)
                        {
                            double value = v(i, j).imag();
                            SwapEnd(value);
                            output_stream.write(reinterpret_cast<char *>(&value), sizeof(double));
                        }
                    }
                }
                else if(strcmp(format,"ASCII")==0)
                {
                    for (size_t j = 0; j < ny; ++j)
                        for (size_t i = 0; i < nx; ++i)
                            output_stream << v(i, j).imag() << " ";
                }

            output_stream.close();

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
         * @param format *char* format of the vtk file, either "ASCII" of "BINARY". Notice that BINARY files occupy less
         * memory and are faster to read/write. However, they may be harder to debug and/or to access via python scripts.
         *
         * The output file will be in standard .vtk format for structured data points.
         * It can be readily visualized using programs like
         * <a href="https://www.paraview.org/"> Paraview </a> or
         * <a href="https://visit-dav.github.io/visit-website/index.html">Visit</a>.
         * \note The data is stored in binary format, and in Big Endian byte ordering.
         */

        void DataWriter::write_vtk(Vector<double>& x,
                                   Vector<double>& y,
                                   Vector<double>& z,
                                   Vector<double>& v,
                                   const char* output_vector_name,
                                   const char* format)
        {

            assert(x.order()==1);
            assert(y.order()==1);
            assert(z.order()==1);
            assert(v.order()==3);
            int nx=x.extent(0);
            int ny=y.extent(0);
            int nz=z.extent(0);
            assert(v.extent(0)==nx);
            assert(v.extent(1)==ny);
            assert(v.extent(2)==nz);
            assert(strcmp(format,"ASCII") || strcmp(format,"BINARY"));

            // Everything is fine, write the output

            double xmin = x(0);
            double ymin = y(0);
            double zmin = z(0);
            double dx = std::abs(x(1)-x(0));
            double dy = std::abs(y(1)-y(0));
            double dz = std::abs(z(1)-z(0));

            output_stream.open(output_file_name+".vtk",std::ios_base::out | std::ios_base::binary);

                // Header of the .vtk file

                output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                output_stream << "# 3D plot of a 3D Real Vector"                    << std::endl;

                // Data type

                if(strcmp(format,"ASCII")==0)
                    output_stream << "ASCII"                                        << std::endl;
                else if (strcmp(format,"BINARY")==0)
                    output_stream << "BINARY"                                       << std::endl;

                // Dataset and mesh information

                output_stream << "DATASET STRUCTURED_POINTS"                         << std::endl;
                output_stream << "DIMENSIONS " << nx   << " " << ny   << " " << nz   << std::endl;
                output_stream << "ORIGIN "     << xmin << " " << ymin << " " << zmin << std::endl;
                output_stream << "SPACING "    << dx   << " " << dy   << " " << dz   << std::endl;

                // Output the vector

                output_stream << "POINT_DATA " << nx*ny*nz                          << std::endl;
                output_stream << "SCALARS "    << output_vector_name << " double 1" << std::endl;
                output_stream << "LOOKUP_TABLE default"                             << std::endl;
                if(strcmp(format,"BINARY")==0)
                {
                    for (size_t k = 0; k < nz; ++k)
                        for (size_t j = 0; j < ny; ++j)
                            for (size_t i = 0; i < nx; ++i)
                            {
                                double value = v(i,j,k);
                                SwapEnd(value);
                                output_stream.write(reinterpret_cast<char*>(&value),sizeof(double));
                            }
                }
                if(strcmp(format,"ASCII")==0)
                {
                    for (size_t k = 0; k < nz; ++k)
                        for (size_t j = 0; j < ny; ++j)
                            for (size_t i = 0; i < nx; ++i)
                            {
                                output_stream << v(i,j,k) << " ";
                            }
                }

            output_stream.close();

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
         * @param format *char* format of the vtk file, either "ASCII" of "BINARY". Notice that BINARY files occupy less
         * memory and are faster to read/write. However, they may be harder to debug and/or to access via python scripts.
         *
         * The output file will be in standard .vtk format for structured data points.
         * It can be readily visualized using programs like
         * <a href="https://www.paraview.org/"> Paraview </a> or
         * <a href="https://visit-dav.github.io/visit-website/index.html">Visit</a>.
         * \note The data is stored in binary format, and in Big Endian byte ordering.
         */

        void DataWriter::write_vtk(Vector<double>& x,
                                   Vector<double>& y,
                                   Vector<double>& z,
                                   Vector<std::complex<double>>& v,
                                   const char* output_vector_name,
                                   const char* format)
        {

            assert(x.order()==1);
            assert(y.order()==1);
            assert(z.order()==1);
            assert(v.order()==3);
            int nx=x.extent(0);
            int ny=y.extent(0);
            int nz=z.extent(0);
            assert(v.extent(0)==nx);
            assert(v.extent(1)==ny);
            assert(v.extent(2)==nz);
            assert(strcmp(format,"ASCII") || strcmp(format,"BINARY"));

            // Everything is fine, write the output

            double xmin = x(0);
            double ymin = y(0);
            double zmin = z(0);
            double dx = std::abs(x(1)-x(0));
            double dy = std::abs(y(1)-y(0));
            double dz = std::abs(z(1)-z(0));

            output_stream.open(output_file_name+".vtk",std::ios_base::out | std::ios_base::binary);

                // Header of the .vtk file

                output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                output_stream << "# 3D plot of a 3D Complex Vector"                 << std::endl;

                // Data type

                if(strcmp(format,"ASCII")==0)
                    output_stream << "ASCII"                                        << std::endl;
                else if (strcmp(format,"BINARY")==0)
                    output_stream << "BINARY"                                       << std::endl;

                // Dataset and mesh information

                output_stream << "DATASET STRUCTURED_POINTS"                         << std::endl;
                output_stream << "DIMENSIONS " << nx   << " " << ny   << " " << nz   << std::endl;
                output_stream << "ORIGIN "     << xmin << " " << ymin << " " << zmin << std::endl;
                output_stream << "SPACING "    << dx   << " " << dy   << " " << dz   << std::endl;

                // Output the real() part of the vector

                output_stream << "POINT_DATA "   << nx*ny*nz                            << std::endl;
                output_stream << "SCALARS real_" << output_vector_name << " double 1"   << std::endl;
                output_stream << "LOOKUP_TABLE default"                                 << std::endl;
                if(strcmp(format,"BINARY")==0)
                {
                    for (size_t k = 0; k < nz; ++k)
                        for (size_t j = 0; j < ny; ++j)
                            for (size_t i = 0; i < nx; ++i)
                            {
                                double value = v(i, j, k).real();
                                SwapEnd(value);
                                output_stream.write(reinterpret_cast<char *>(&value), sizeof(double));
                            }
                }
                if(strcmp(format,"ASCII")==0)
                {
                    for (size_t k = 0; k < nz; ++k)
                        for (size_t j = 0; j < ny; ++j)
                            for (size_t i = 0; i < nx; ++i)
                                output_stream << v(i, j, k).real() << " ";
                }

                // Output the imaginary part of the vector

                output_stream << std::endl;
                output_stream << "SCALARS imag_" << output_vector_name << " double 1"   << std::endl;
                output_stream << "LOOKUP_TABLE default"                                 << std::endl;
                if(strcmp(format,"BINARY")==0)
                {
                    for (size_t k = 0; k < nz; ++k)
                        for (size_t j = 0; j < ny; ++j)
                            for (size_t i = 0; i < nx; ++i)
                            {
                                double value = v(i, j, k).imag();
                                SwapEnd(value);
                                output_stream.write(reinterpret_cast<char *>(&value), sizeof(double));
                            }
                }
                if(strcmp(format,"ASCII")==0)
                {
                    for (size_t k = 0; k < nz; ++k)
                        for (size_t j = 0; j < ny; ++j)
                            for (size_t i = 0; i < nx; ++i)
                                output_stream << v(i, j, k).imag() << " ";
                }
            output_stream.close();

        }

        /**
         * @brief Write an output data file in .vtk format, for 2D slice of real 3D Vector
         * @param ax1 *Vector<double>* the first Vector for the output.
         * This is considered as the ax1-axis for the plot.
         * @param ax2 *Vector<double>* the second Vector for the output.
         * This is considered as the ax2-axis for the plot.
         * @param v *Vector<double>* the third Vector for the output.
         * One of its cuts is considered as the z-axis of the plot.
         * @param plane the plane along which the cut is taken. This can either be <code>"xy"</code>,
         * <code>"yz"</code>, or <code>"xz"</code>. If <code>plane="xy"</code>, a slice along
         * the z=0 plane will be taken, and similarly if <code>axis="yz"</code> or <code>axis="xz"</code>.
         * @param output_vector_name *char* the name of the output
         * @param format *char* format of the vtk file, either "ASCII" of "BINARY". Notice that BINARY files occupy less
         * memory and are faster to read/write. However, they may be harder to debug and/or to access via python scripts.
         *
         * The output file will be in standard .vtk format for structured data points.
         * It can be readily visualized using programs like
         * <a href="https://www.paraview.org/"> Paraview </a> or
         * <a href="https://visit-dav.github.io/visit-website/index.html">Visit</a>.
         * \note The data is stored in binary format, and in Big Endian byte ordering.
         */

        void DataWriter::write_slice2d_vtk(Vector<double>& ax1,
                                           Vector<double>& ax2,
                                           Vector<double>& v,
                                           const char* plane,
                                           const char* output_vector_name,
                                           const char* format)
        {

            // Check the extents are correct
            assert(ax1.order()==1);
            assert(ax2.order()==1);
            assert(v.order()==3);
            assert(strcmp(plane,"xy") || strcmp(plane,"xz") || strcmp(plane,"yz"));
            assert(strcmp(format,"ASCII") || strcmp(format,"BINARY"));

            // slice of v along the z=0 plane

            if(strcmp(plane,"xy")==0)
            {
                int nx = ax1.extent(0);
                int ny = ax2.extent(0);
                int nz = v.extent(2);

                assert(v.extent(0) == nx);
                assert(v.extent(1) == ny);

                // Everything is fine, write the output

                double xmin = ax1(0);
                double ymin = ax2(0);
                double dx = std::abs(ax1(1) - ax1(0));
                double dy = std::abs(ax2(1) - ax2(0));

                output_stream.open(output_file_name+".vtk",std::ios_base::out | std::ios_base::binary);

                    // Header of the .vtk file

                    output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                    output_stream << "# 2D slice of a 3D Real Vector"                   << std::endl;

                    // Data type

                    if(strcmp(format,"ASCII")==0)
                        output_stream << "ASCII"                                        << std::endl;
                    else if (strcmp(format,"BINARY")==0)
                        output_stream << "BINARY"                                       << std::endl;

                    // Dataset and mesh information

                    output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                    output_stream << "DIMENSIONS " << nx   << " " << ny   << " " << "1" << std::endl;
                    output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                    output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                    // Output the vector

                    output_stream << "POINT_DATA " << nx*ny                             << std::endl;
                    output_stream << "SCALARS "    << output_vector_name << " double 1" << std::endl;
                    output_stream << "LOOKUP_TABLE default"                             << std::endl;
                    if(strcmp(format,"BINARY")==0)
                    {
                        for (size_t j = 0; j < ny; ++j)
                            for (size_t i = 0; i < nx; ++i) {
                                double value = v(i, j, (int) nz / 2);
                                SwapEnd(value);
                                output_stream.write(reinterpret_cast<char *>(&value), sizeof(double));
                            }
                    }
                    if(strcmp(format,"ASCII")==0)
                    {
                        for (size_t j = 0; j < ny; ++j)
                            for (size_t i = 0; i < nx; ++i)
                                output_stream << v(i, j, (int) nz / 2) << " ";
                    }

                output_stream.close();

            }

            // Slice of v along the x=0 plane

            else if(strcmp(plane,"yz")==0)
            {
                int nx=v.extent(0);
                int ny=ax1.extent(0);
                int nz=ax2.extent(0);
                assert(v.extent(1)==ny);
                assert(v.extent(2)==nz);

                // Everything is fine, write the output

                double xmin = ax1(0);
                double ymin = ax2(0);
                double dx = std::abs(ax1(1) - ax1(0));
                double dy = std::abs(ax2(1) - ax2(0));

                output_stream.open(output_file_name+".vtk",std::ios_base::out | std::ios_base::binary);

                    // Header of the .vtk file

                    output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                    output_stream << "# 2D slice of a 3D Real Vector"                   << std::endl;

                    // Data type is text

                    if(strcmp(format,"ASCII")==0)
                        output_stream << "ASCII"                                        << std::endl;
                    else if (strcmp(format,"BINARY")==0)
                        output_stream << "BINARY"                                       << std::endl;

                    // Dataset and mesh information

                    output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                    output_stream << "DIMENSIONS " << ny   << " " << nz  << " " << "1" << std::endl;
                    output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                    output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                    // Output the vector

                    output_stream << "POINT_DATA " << ny*nz                             << std::endl;
                    output_stream << "SCALARS "    << output_vector_name << " double 1" << std::endl;
                    output_stream << "LOOKUP_TABLE default"                             << std::endl;
                    if(strcmp(format,"BINARY")==0)
                    {
                        for (size_t j = 0; j < ny; ++j)
                            for (size_t i = 0; i < nz; ++i) {
                                double value = v((int) nx / 2, i, j);
                                SwapEnd(value);
                                output_stream.write(reinterpret_cast<char *>(&value), sizeof(double));
                            }
                    }
                    if(strcmp(format,"ASCII")==0)
                    {
                        for (size_t j = 0; j < ny; ++j)
                            for (size_t i = 0; i < nz; ++i)
                                output_stream << v((int) nx / 2, i, j) << " ";
                    }

                output_stream.close();
            }

            // Slice of v along the y=0 plane

            else if(strcmp(plane,"xz")==0)
            {
                int nx=ax1.extent(0);
                int ny=v.extent(1);
                int nz=ax2.extent(0);
                assert(v.extent(0)==nx);
                assert(v.extent(2)==nz);

                // Everything is fine, write the output

                double xmin = ax1(0);
                double ymin = ax2(0);
                double dx = std::abs(ax1(1) - ax1(0));
                double dy = std::abs(ax2(1) - ax2(0));

                output_stream.open(output_file_name+".vtk",std::ios_base::out | std::ios_base::binary);

                    // Header of the .vtk file

                    output_stream << "# vtk DataFile Version 3.0"                         << std::endl;
                    output_stream << "# 2D slice of a 3D real Vector"                     << std::endl;

                    // Data type

                    if(strcmp(format,"ASCII")==0)
                        output_stream << "ASCII"                                        << std::endl;
                    else if (strcmp(format,"BINARY")==0)
                        output_stream << "BINARY"                                       << std::endl;

                    // Dataset and mesh information

                    output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                    output_stream << "DIMENSIONS " << nx   << " " << nz   << " " << "1" << std::endl;
                    output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                    output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                    // Output the vector

                    output_stream << "POINT_DATA " << nx*nz                           << std::endl;
                    output_stream << "SCALARS "    << output_vector_name << " double 1" << std::endl;
                    output_stream << "LOOKUP_TABLE default"                             << std::endl;

                    if(strcmp(format,"BINARY")==0)
                    {
                        for (size_t j = 0; j < nz; ++j)
                            for (size_t i = 0; i < nx; ++i)
                            {
                                double value = v(i,(int) ny/2,j);
                                SwapEnd(value);
                                output_stream.write(reinterpret_cast<char*>(&value),sizeof(double));
                            }
                    }
                    if(strcmp(format,"ASCII") == 0)
                    {
                        for (size_t j = 0; j < nz; ++j)
                            for (size_t i = 0; i < nx; ++i)
                                output_stream << v(i,(int) ny/2,j) << " ";
                    }

                output_stream.close();

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
         * @param format *char* format of the vtk file, either "ASCII" of "BINARY". Notice that BINARY files occupy less
         * memory and are faster to read/write. However, they may be harder to debug and/or to access via python scripts.
         *
         * The output file will be in standard .vtk format for structured data points.
         * It can be readily visualized using programs like
         * <a href="https://www.paraview.org/"> Paraview </a> or
         * <a href="https://visit-dav.github.io/visit-website/index.html">Visit</a>.
         * \note The data is stored in binary format, and in Big Endian byte ordering.
         */

        void DataWriter::write_slice2d_vtk(Vector<double>& ax1,
                                           Vector<double>& ax2,
                                           Vector<std::complex<double>>& v,
                                           const char* plane,
                                           const char* output_vector_name,
                                           const char* format)
        {

            // Check the extents are correct
            assert(ax1.order()==1);
            assert(ax2.order()==1);
            assert(v.order()==3);
            assert(strcmp(plane,"xy") || strcmp(plane,"xz") || strcmp(plane,"yz"));
            assert(strcmp(format,"ASCII") || strcmp(format,"BINARY"));

            // slice of v along the z=0 plane

            if(strcmp(plane,"xy")==0)
            {
                int nx = ax1.extent(0);
                int ny = ax2.extent(0);
                int nz = v.extent(2);

                assert(v.extent(0) == nx);
                assert(v.extent(1) == ny);

                // Everything is fine, write the output

                double xmin = ax1(0);
                double ymin = ax2(0);
                double dx = std::abs(ax1(1) - ax1(0));
                double dy = std::abs(ax2(1) - ax2(0));

                output_stream.open(output_file_name+".vtk",std::ios_base::out | std::ios_base::binary);

                // Header of the .vtk file

                output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                output_stream << "# 2D slice of a 3D Real Vector"                   << std::endl;

                // Data type

                if(strcmp(format,"ASCII")==0)
                    output_stream << "ASCII"                                        << std::endl;
                else if (strcmp(format,"BINARY")==0)
                    output_stream << "BINARY"                                       << std::endl;

                // Dataset and mesh information

                output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                output_stream << "DIMENSIONS " << nx   << " " << ny   << " " << "1" << std::endl;
                output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                // Output the vector

                output_stream << "POINT_DATA "      << nx*ny                             << std::endl;

                // Real part of the input vector

                output_stream << "SCALARS real_"    << output_vector_name << " double 1" << std::endl;
                output_stream << "LOOKUP_TABLE default"                                  << std::endl;

                if(strcmp(format,"BINARY")==0)
                {
                    for (size_t j = 0; j < ny; ++j)
                        for (size_t i = 0; i < nx; ++i)
                        {
                            double value = v(i,j,(int)nz/2).real();
                            SwapEnd(value);
                            output_stream.write(reinterpret_cast<char*>(&value),sizeof(double));
                        }
                }
                if(strcmp(format,"ASCII")==0)
                {
                    for (size_t j = 0; j < ny; ++j)
                        for (size_t i = 0; i < nx; ++i)
                            output_stream << v(i,j,(int)nz/2).real() << " ";
                }

                output_stream << std::endl;

                // Imaginary part of the input vector
                output_stream << "SCALARS imag_"    << output_vector_name << " double 1" << std::endl;
                output_stream << "LOOKUP_TABLE default"                                  << std::endl;
                if(strcmp(format,"BINARY")==0)
                {
                    for (size_t j = 0; j < ny; ++j)
                        for (size_t i = 0; i < nx; ++i)
                        {
                            double value = v(i,j,(int)nz/2).imag();
                            SwapEnd(value);
                            output_stream.write(reinterpret_cast<char*>(&value),sizeof(double));
                        }
                }
                if(strcmp(format,"ASCII")==0)
                {
                    for (size_t j = 0; j < ny; ++j)
                        for (size_t i = 0; i < nx; ++i)
                            output_stream << v(i,j,(int)nz/2).imag() << " ";
                }

                output_stream.close();

            }

            // Slice of v along the x=0 plane

            else if(strcmp(plane,"yz")==0)
            {
                int nx=v.extent(0);
                int ny=ax1.extent(0);
                int nz=ax2.extent(0);
                assert(v.extent(1)==ny);
                assert(v.extent(2)==nz);

                // Everything is fine, write the output

                double xmin = ax1(0);
                double ymin = ax2(0);
                double dx = std::abs(ax1(1) - ax1(0));
                double dy = std::abs(ax2(1) - ax2(0));

                output_stream.open(output_file_name+".vtk",std::ios_base::out | std::ios_base::binary);

                // Header of the .vtk file

                output_stream << "# vtk DataFile Version 3.0"                       << std::endl;
                output_stream << "# 2D slice of a 3D Real Vector"                   << std::endl;

                // Data type is text

                if(strcmp(format,"ASCII")==0)
                    output_stream << "ASCII"                                        << std::endl;
                else if (strcmp(format,"BINARY")==0)
                    output_stream << "BINARY"                                       << std::endl;

                // Dataset and mesh information

                output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                output_stream << "DIMENSIONS " << ny   << " " << nz  << " " << "1" << std::endl;
                output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                // Output the vector

                output_stream << "POINT_DATA "      << ny*nz                             << std::endl;

                // Real part of the input vector

                output_stream << "SCALARS real_"    << output_vector_name << " double 1" << std::endl;
                output_stream << "LOOKUP_TABLE default"                                  << std::endl;

                if(strcmp(format,"BINARY")==0)
                {
                    for (size_t j = 0; j < ny; ++j)
                        for (size_t i = 0; i < nz; ++i)
                        {
                            double value = v((int) nx / 2, i, j).real();
                            SwapEnd(value);
                            output_stream.write(reinterpret_cast<char *>(&value), sizeof(double));
                        }
                }
                if(strcmp(format,"ASCII")==0)
                {
                    for (size_t j = 0; j < ny; ++j)
                        for (size_t i = 0; i < nz; ++i)
                            output_stream << v((int) nx / 2, i, j).real() << " ";
                }

                output_stream << std::endl;

                // Imaginary part of the input vector

                output_stream << "SCALARS imag_"    << output_vector_name << " double 1" << std::endl;
                output_stream << "LOOKUP_TABLE default"                                  << std::endl;
                if(strcmp(format,"BINARY")==0)
                {
                    for (size_t j = 0; j < ny; ++j)
                        for (size_t i = 0; i < nz; ++i)
                        {
                            double value = v((int) nx / 2, i, j).imag();
                            SwapEnd(value);
                            output_stream.write(reinterpret_cast<char *>(&value), sizeof(double));
                        }
                }
                if(strcmp(format,"ASCII")==0)
                {
                    for (size_t j = 0; j < ny; ++j)
                        for (size_t i = 0; i < nz; ++i)
                            output_stream << v((int) nx / 2, i, j).imag() << " ";
                }

                output_stream.close();
            }

             // Slice of v along the y=0 plane

            else if(strcmp(plane,"xz")==0)
            {
                int nx=ax1.extent(0);
                int ny=v.extent(1);
                int nz=ax2.extent(0);
                assert(v.extent(0)==nx);
                assert(v.extent(2)==nz);

                // Everything is fine, write the output

                double xmin = ax1(0);
                double ymin = ax2(0);
                double dx = std::abs(ax1(1) - ax1(0));
                double dy = std::abs(ax2(1) - ax2(0));

                output_stream.open(output_file_name+".vtk",std::ios_base::out | std::ios_base::binary);

                // Header of the .vtk file

                output_stream << "# vtk DataFile Version 3.0"                         << std::endl;
                output_stream << "# 2D slice of a 3D real Vector"                     << std::endl;

                // Data type is text

                if(strcmp(format,"ASCII")==0)
                    output_stream << "ASCII"                                        << std::endl;
                else if (strcmp(format,"BINARY")==0)
                    output_stream << "BINARY"                                       << std::endl;

                // Dataset and mesh information

                output_stream << "DATASET STRUCTURED_POINTS"                        << std::endl;
                output_stream << "DIMENSIONS " << nx   << " " << nz   << " " << "1" << std::endl;
                output_stream << "ORIGIN "     << xmin << " " << ymin << " " << "0" << std::endl;
                output_stream << "SPACING "    << dx   << " " << dy   << " " << "0" << std::endl;

                // Output the vector

                output_stream << "POINT_DATA "      << nx*nz                             << std::endl;

                // Real part of the input vector

                output_stream << "SCALARS real_"    << output_vector_name << " double 1" << std::endl;
                output_stream << "LOOKUP_TABLE default"                                  << std::endl;
                if(strcmp(format,"BINARY")==0)
                {
                    for (size_t j = 0; j < nz; ++j)
                        for (size_t i = 0; i < nx; ++i)
                        {
                            double value = v(i,(int) ny/2,j).real();
                            SwapEnd(value);
                            output_stream.write(reinterpret_cast<char*>(&value),sizeof(double));
                        }
                }
                if(strcmp(format,"ASCII")==0)
                {
                    for (size_t j = 0; j < nz; ++j)
                        for (size_t i = 0; i < nx; ++i)
                            output_stream << v(i,(int) ny/2,j).real() << " ";
                }

                output_stream << std::endl;

                // Imaginary part of the input vector

                output_stream << "SCALARS imag_"    << output_vector_name << " double 1" << std::endl;
                output_stream << "LOOKUP_TABLE default"                                  << std::endl;
                if(strcmp(format,"BINARY")==0)
                {
                    for (size_t j = 0; j < nz; ++j)
                        for (size_t i = 0; i < nx; ++i)
                        {
                            double value = v(i,(int) ny/2,j).imag();
                            SwapEnd(value);
                            output_stream.write(reinterpret_cast<char*>(&value),sizeof(double));
                        }
                }
                if(strcmp(format,"ASCII")==0)
                {
                    for (size_t j = 0; j < nz; ++j)
                        for (size_t i = 0; i < nx; ++i)
                            output_stream << v(i,(int) ny/2,j).imag() << " ";
                }

                output_stream.close();

            }
        }

    } // namespace GraphicOutput
} // namespace UltraCold