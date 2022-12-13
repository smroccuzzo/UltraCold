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

#include "HardwareInspector.hpp"

namespace UltraCold{

    namespace Tools
    {
        /**
         * @brief Constructor of class HardwareInspector.
         * 
         * The constructor actually retrieves the hardware information on the machines and initializes the
         * corresponding member variables.
         * 
         */

        HardwareInspector::HardwareInspector()
        {

            number_of_processors = get_nprocs_conf();
            number_of_available_processors = get_nprocs();

        }

        /**
         * @brief Returns the total number of processors present in the machines
         * 
         * Notice that this number may be different from the actual number of processors available for the program.
         * 
         * @return the total number of processors present on the machine
         * 
         */

        int HardwareInspector::get_number_of_processors()
        {

            return number_of_processors;

        }

        /**
         * @brief Returns the number of available processors
         * 
         * Notice that this number may be different from the total number of processors present on the machine.
         * 
         * @return the number of available processors
         * 
         */

        int HardwareInspector::get_number_of_available_processors()
        {

            return number_of_available_processors;
        }

        /**
         * @brief Print all available hardware information on the screen
         * 
         */

        void HardwareInspector::print_cpu_information()
        {

            std::cout 
                << "Number of processors:\t" 
                << number_of_processors    
                <<
            std::endl;

            std::cout
                << "Number of available processors:\t"
                << number_of_available_processors
                <<
            std::endl;
        }
    } // namespace Tools
    
} // namespace UltraCold