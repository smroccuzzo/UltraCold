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

#ifndef ULTRACOLD_HARDWARE_INFORMATION
#define ULTRACOLD_HARDWARE_INFORMATION

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/sysinfo.h>

namespace UltraCold
{

    /**
     * @brief Classes and functions of general utility, from input parsers to hardware inspectors.
     *
     */

    namespace Tools
    {

        /**
         * @brief Class to detect hardware capabilities
         * @author Santo Maria Roccuzzo (santom.roccuzzo@gmail.com)
         *
         * The HardwareInspector class provides methods to check hardware information on the current machine, such as
         * the number of available cores, or the total amount of available memory. It also helps in tuning  the
         * performance of UltraCold.
         *
         */

        class HardwareInspector
        {
            public:

                HardwareInspector();
                ~HardwareInspector() = default;

                // Returns the total number of processors present on the machine

                int get_number_of_processors();

                // Returns the actual number of available processors

                int get_number_of_available_processors();

                // Print all this information to the screen

                void print_cpu_information();

            private:

                int number_of_processors;
                int number_of_available_processors;

        };
    } // namespace Tools
} // namespace UltraCold

#endif // ULTRACOLD_HARDWARE_INFORMATION