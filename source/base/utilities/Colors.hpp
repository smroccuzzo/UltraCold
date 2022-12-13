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

#ifndef ULTRACOLD_COLORS
#define ULTRACOLD_COLORS

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <map>
#include <cassert>

#include "Vector.hpp"

namespace UltraCold
{
    namespace GraphicOutput
    {

        /**
         *
         * @brief A simple representation of a pixel as an RGB color. Useful for certain types of visualization outputs.
         *
         * */

        struct Pixel
        {
            unsigned char r;
            unsigned char g;
            unsigned char b;
        };

        /**
         *
         * @brief A class representing a color in an RGB scale
         *
         * This class is useful for creating custom color maps, to be used for example to output two-dimensional
         * density images directly as .ppm images. See for example the function create_linear_color_map()
         *
         **/

        class Color
        {

        public:

            // Default constructor
            Color();

            // Copy constructor
            Color(const Color& other);

            // Copy assignment
            Color& operator=(const Color& other);

            // Move constructor
            Color(const Color&& other) noexcept;

            // Move assignment
            Color& operator=(Color&& other) noexcept;

            // Constructor from the amount of RGB
            Color(double r_in,double g_in,double b_in);

            // Destructor is default
            ~Color()=default;

            // Setters
            void set_color(int r_in,int g_in,int b_in);

            // Sum of two colors
            Color operator+(const Color& other) const;

            // Scale a color
            Color operator*(const double p) const;

            // Getters
            double get_red()   const;
            double get_green() const;
            double get_blue()  const;

            Pixel get_pixel() const;

        private:

            double r;
            double g;
            double b;

            Pixel pixel;
        };

        // Output a color
        std::ostream& operator<<(std::ostream& stream, Color& output_color);


        // Create a linear color map
        std::map<double,Color> create_linear_color_map(double v_min,
                                                       double v_max,
                                                       Color color_min,
                                                       Color color_max,
                                                       int number_of_levels);

        // Create a diverging color map
        std::map<double,Color> create_diverging_color_map(double v_min,
                                                          double v_mid,
                                                          double v_max,
                                                          Color color_min,
                                                          Color color_mid,
                                                          Color color_max,
                                                          int number_of_levels);

        // Write a 2-dimensional double Vector directly as a ppm image
        void write_ppm(const char* filename,
                       Vector<double> values,
                       std::map<double,Color> color_map);

    } // Tools
} // UltraCold

#endif // ULTRACOLD_COLORS