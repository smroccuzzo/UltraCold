/*--------------------------------------------------------------------------------
 *
 *    This file is part of the UltraCold project.
 *
 *    UltraCols is free software: you can redistribute it and/or modify
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

#ifndef ULTRACOLD_INPUT_PARSER
#define ULTRACOLD_INPUT_PARSER

#include <iostream>
#include <fstream>
#include <string>    
#include <algorithm>
#include <map>

namespace UltraCold
{

    namespace Tools
    {

        /**
         * @brief Class to read input parameters from files.
         * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
         *
         * The InputParser class provides an interface to input parameter files. The class can be used as follows. \n
         * Suppose you have a file input_file.prm containing the following text
         * 
             \verbatim

            # Define the mesh parameters

            xmax = 10.0 # Size of the mesh along the x-axis, in micrometers.
            ymax = 10.0 # Size of the mesh along the y-axis, in micrometers.
            zmax = 10.0 # Size of the mesh along the z-axis, in micrometers.

            nx = 100 # number of points along the x-axis
            ny = 100 # number of points along the y-axis
            nz = 100 # number of points along the z-axis

            # Define interaction parameters

            scattering length = 100.0 # Scattering length in micrometers

            # Define the run mode

            calculate ground state = true # If true, calculate a ground state solution of the eGPE

            \endverbatim
        * 
        * Comments are defined by an hashtag (#) and ignored by the class. \n
        * Every line containing an equal (=) sign is interpreted as a line defining an input. Every character on the
        * left of the equal sign is interpreted as the name of the variable (without blanks), while what follows the
         * equal sign on the right is interpreted as the value of the defined variable. \n
        * In order to retrieve the value of these variables and use them in the code, it is necessary 
        * to use the appropriate member function <code><variable type> InputParser::retrieve_<type>("variable name")
        * </code>, which will cast the retrieved variable (which is first interpreted as an
        * <code>std::string</code>)  to the desired <code> <type> </code>. \n
        * So, an example code block that reads the input file defined above and initializes the variables 
        * *xmax*, *ymax*, *zmax*, *scattering length*, and *run in imaginary time*, is the following
        *    
        * \code {.cpp}
        * 
        * InputParser ip("input_file.prm");
        * 
        * ip.read_input_file();
        * 
        * const double xmax = ip.retrieve_double("xmax");
        * const double ymax = ip.retrieve_double("ymax");
        * const double zmax = ip.retrieve_double("zmax");
        * const int nx = ip.retrieve_int("nx");
        * const int ny = ip.retrieve_int("ny");
        * const int nz = ip.retrieve_int("nz");
        * const double scattering_length = ip.retrieve_double("scattering length"));
        * const bool calc_ground_state = ip.retrieve_bool("calculate ground state");
        *
        * \endcode
        * 
        */ 

        class InputParser
        {
            public:  

                InputParser(char* input_file_name);
                ~InputParser();

                // Read input file and initialize the input_list map 
                void read_input_file();

                // Retrieve inputs casting to the appropriate type

                int retrieve_int(const char* variable_name);
                double retrieve_double(const char* variable_name);
                bool retrieve_bool(const char* variable_name);

            private:

                // Retrieve values from the input_list map, cast as string

                std::string retrieve_input(const char*);

                // Internal variable members 

                char* input_file_name;
                std::ifstream ifs;
                std::map<std::string,std::string> input_list;

        };

    } // namespace Tools
} // namespace UltraCold

#endif // ULTRACOLD_INPUT_PARSER
