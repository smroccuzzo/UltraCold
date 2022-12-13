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

#include "InputParser.hpp"

namespace UltraCold{

    namespace Tools
    {
        /**
         * 
         * @brief Constructor of class InputParser.
         * @param input_file_name *char*  the name of the input file
         * 
         * The constructor tries to open the file <code> input_file_name </code>. 
         * If the file is not found, it will print an error message an terminate the execution of the program. If
         * instead the file is found, the constructor initializes an <code>std::ifstream</code> through which the
         * file can be read.
         * 
        */

        InputParser::InputParser(char* input_file_name)
        {
            
            this->input_file_name = input_file_name;
            ifs.open(this->input_file_name,std::ios::in);
            if (!ifs.is_open())
            {
                std::cout
                    << "\n\n"
                    << "*****************************************************\n"
                    << "Input file " << this->input_file_name << " not found!\n"
                    << "Terminating the program now\n"
                    << "*****************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                std::exit(1);
            }
        }

        /**
         * @brief Destructor
         * 
         * The destructor simply closes the <code>std::ifstream</code> used to read the input file.
         * 
        */

        InputParser::~InputParser()
        {
            ifs.close();
        }


        /**
         * @brief Member function to parse the input file
         * 
         * This function parses the input file, ignoring all comments (i.e., everything coming after an hashtag (#)),
         * and filling an <code>std::map<std::string,std::string></code> with the pairs name=value defined
         * by lines containing an equal (=) sign. \n
         * 
        */

        void InputParser::read_input_file()
        {

            std::string input_line;
            size_t i_line;

            // Reset line counter

            i_line = 1; 
            
            // Read the file line by line

            while(std::getline(ifs,input_line))
            {

                // Remove eventual blanks

                input_line.erase
                (
                    std::remove(input_line.begin(),input_line.end(),' '),
                    input_line.end()
                );

                // If the line is not blank, then we can process it

                if( input_line.length() != 0 )
                {

                    // Remove eventual comments
                    
                    size_t i_hashtag = input_line.find("#");
                    if( i_hashtag != -1 )
                    {
                        input_line = input_line.substr(0,i_hashtag);
                    }

                    // Check if the current line defines an input, i.e. contains an = sign. 
                    // In this case, what comes before the equal sign is the name of the variable, 
                    // and what is on the right is its value. This is always interpreted as a string, 
                    // and must be explicitly cast to the desired type via the appropriate member function.

                    size_t i_equal = input_line.find("=");

                    // If neither an equal sign or an hashtag are found, then the line format is invalid,
                    // as the line does not define a comment nor a parameter. In this case, throw an error
                    // and terminate the execution of the program. 

                    if( (i_equal==-1) && (i_hashtag==-1) )
                    {

                        std::cout
                            << "\n\n"
                            << "*********************************************************************\n"
                            << "Error found in line "
                            << i_line
                            << " while parsing the input file " << input_file_name <<  ". \n"
                            << "Line " 
                            << i_line 
                            << " does not define either an input (does not contain an equal sign)\n" 
                            << "nor a comment (does not contain an hashtag). Please fix this \n"
                            << "in the input file.\n"
                            << "Terminating the program now.\n"
                            << "*********************************************************************\n"
                            << "\n\n"
                            <<
                        std::endl;
                        std::exit(1);
                    }

                    // If an input instead is found, split its name from its value and put them in the input_list map

                    if( i_equal !=-1 )
                    {

                        std::string input_name  = input_line.substr(0,i_equal);
                        std::string input_value = input_line.substr(i_equal+1,input_line.length());

                        input_list.emplace(input_name,input_value);

                    }

                }

                // update line counter

                i_line = i_line+1;

            }

        }

        /**
         * 
         * @brief Retrieve input variables from input file
         * 
         * The input name of the requested element must be provided as an array of <code>char</code>, in the form
         * "name of the variable". Internally, this name is converted to a string, eventual blanks are removed, and
         * then the name of the variable is searched in the input_list map created with
         * <code>InputParser::read_input_file()</code>. If the name of the element is found, its value is returned as
         * an <code>std::string</code>, that must then be cast to the desired type using one of the other member
         * functions.
         * 
         * @param requested_element *char** Name of the element to be retrieved.
         * 
         * @return Value of the requested element, as <code>std::string</code>, that must be explicitly 
         * cast to the desired type.
         * 
        */

        std::string InputParser::retrieve_input(const char* requested_element)
        {

            // Cast requested_element to string and remove eventual blanks from 
            // the requested element name

            std::string req_element = requested_element;

            req_element.erase
            (
                std::remove(req_element.begin(),req_element.end(),' '),
                req_element.end()
            );    

            // Now search for req_element in the input list

            const auto element_iterator = input_list.find(req_element);

            // If the element is not found, print an error message and terminate 
            // the program...

            if(element_iterator==input_list.end())
            {
                std::cout
                    << "\n\n"
                    << "************************************************************\n"
                    << "Requested an element not found in " << input_file_name << ".\n"
                    << "The requested element was\n"
                    << "\t" << requested_element << "\t\t\n"
                    << "Please add\n"
                    << "\t" << requested_element << " = value\n"
                    << "somewhere in " << input_file_name << "\n"
                    << "Terminating the program now.\n"
                    << "************************************************************\n"
                    <<
                std::endl;
                std::exit(1);
            }

            // ... otherwise, return its value as a string.

            return element_iterator->second; 
            
        }        

        /**
         * @brief Cast an element retrieved to int
         * @param requested_element *char** Name of the element to be retrieved.
         * @return the requested element, cast to int
         * 
         */

        int InputParser::retrieve_int(const char* requested_element)
        {
            int element_retrieved;
            element_retrieved = std::stoi(this->retrieve_input(requested_element));
            return element_retrieved;
        }

        /**
         * @brief Cast an element retrieved to double
         * @param requested_element *char** Name of the element to be retrieved.
         * @return the requested element, cast to double
         * 
         */

        double InputParser::retrieve_double(const char* requested_element)
        {
            double element_retrieved;
            element_retrieved = std::stod(this->retrieve_input(requested_element));
            return element_retrieved;
        }

        /**
         * @brief Cast an element retrieved to bool
         * @param requested_element *char** Name of the element to be retrieved.
         * @return the requested element, cast to bool
         * 
         */

        bool InputParser::retrieve_bool(const char* requested_element)
        {
            bool element_retrieved;
            element_retrieved = ("true" == this->retrieve_input(requested_element));
            return element_retrieved;
        }

    } // namespace Tools
} // namespace UltraCold