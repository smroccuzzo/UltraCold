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

#include "Colors.hpp"

namespace UltraCold
{
    namespace GraphicOutput
    {

        /**
         * @brief Default constructor. Initialize to white
         *
         */

        Color::Color() :
                    r(1.0),
                    g(1.0),
                    b(1.0),
                    pixel{255,255,255}{}

        /**
         * @brief Copy constructor
         *
         * */

        Color::Color(const Color& other) :
                    r((assert(other.r >= 0.0 && other.r <= 1.0),other.r)),
                    g((assert(other.g >= 0.0 && other.g <= 1.0),other.g)),
                    b((assert(other.b >= 0.0 && other.b <= 1.0),other.b)),
                    pixel{(assert(other.pixel.r >= 0 && other.pixel.r <= 255),other.pixel.r),
                          (assert(other.pixel.g >= 0 && other.pixel.g <= 255),other.pixel.g),
                          (assert(other.pixel.b >= 0 && other.pixel.b <= 255),other.pixel.b)} {}

        /**
         * @brief Copy assignment
         *
         */

        Color& Color::operator=(const Color& other)
        {
            assert(other.r >= 0.0 && other.r <= 1.0);
            assert(other.g >= 0.0 && other.g <= 1.0);
            assert(other.b >= 0.0 && other.b <= 1.0);
            assert(other.pixel.r >= 0 && other.pixel.r <= 255);
            assert(other.pixel.g >= 0 && other.pixel.g <= 255);
            assert(other.pixel.b >= 0 && other.pixel.b <= 255);

            r=other.r;
            g=other.g;
            b=other.b;
            pixel.r=other.pixel.r;
            pixel.g=other.pixel.g;
            pixel.b=other.pixel.b;

            return *this;
        }

        /**
         * @brief Move constructor
         *
         */

        Color::Color(const Color&& other) noexcept
        {
            assert(other.r >= 0.0 && other.r <= 1.0); // Check that the color provided contains red   in the range [0,1]
            assert(other.g >= 0.0 && other.g <= 1.0); // Check that the color provided contains green in the range [0,1]
            assert(other.b >= 0.0 && other.b <= 1.0); // Check that the color provided contains blue  in the range [0,1]

            assert(other.pixel.r >= 0 && other.pixel.r <= 255);
            assert(other.pixel.g >= 0 && other.pixel.g <= 255);
            assert(other.pixel.b >= 0 && other.pixel.b <= 255);

            r=other.r;
            g=other.g;
            b=other.b;

            pixel.r=other.pixel.r;
            pixel.g=other.pixel.g;
            pixel.b=other.pixel.b;
        }

        /**
         * @brief Move assignment
         *
         */

        Color& Color::operator=(Color&& other) noexcept
        {
            assert(other.r >= 0.0 && other.r <= 1.0); // Check the color provided contains red   in the range [0,1]
            assert(other.g >= 0.0 && other.g <= 1.0); // Check the color provided contains green in the range [0,1]
            assert(other.b >= 0.0 && other.b <= 1.0); // Check the color provided contains blue  in the range [0,1]

            assert(other.pixel.r >= 0 && other.pixel.r <= 255);
            assert(other.pixel.g >= 0 && other.pixel.g <= 255);
            assert(other.pixel.b >= 0 && other.pixel.b <= 255);

            r=other.r;
            g=other.g;
            b=other.b;

            pixel.r=other.pixel.r;
            pixel.g=other.pixel.g;
            pixel.b=other.pixel.b;
            return *this;
        }

        /**
         * @brief Constructor from the amount of RGB
         *
         */

        Color::Color(double r_in,double g_in, double b_in) :
                r((assert(r_in <= 1 && r_in >= 0),r_in)),
                g((assert(g_in <= 1 && g_in >= 0),g_in)),
                b((assert(b_in <= 1 && b_in >= 0),b_in)),
                pixel{static_cast<unsigned char>(floor(r_in >= 1.0 ? 255 : r_in * 256.0)),
                      static_cast<unsigned char>(floor(g_in >= 1.0 ? 255 : g_in * 256.0)),
                      static_cast<unsigned char>(floor(b_in >= 1.0 ? 255 : b_in * 256.0))} {}

        /**
         * @brief Set color from the amount of RGB
         *
         */

        void Color::set_color(int r_in,int g_in,int b_in)
        {
            assert(r_in <= 1 && r_in >= 0);
            assert(g_in <= 1 && g_in >= 0);
            assert(b_in <= 1 && b_in >= 0);

            r=r_in;
            g=g_in;
            b=b_in;

        }

        /**
         * @brief Calculate the sum of two colors
         *
         */

        Color Color::operator+(const Color& other) const
        {
            double return_r = ( (r+other.r > 1.0) ? 1.0 : (r+other.r) < 0.0 ? 0.0 : r+other.r);
            double return_g = ( (g+other.g > 1.0) ? 1.0 : (g+other.g) < 0.0 ? 0.0 : g+other.g);
            double return_b = ( (b+other.b > 1.0) ? 1.0 : (b+other.b) < 0.0 ? 0.0 : b+other.b);

            Color return_color(return_r,return_g,return_b);

            return return_color;

        }

        /**
         * @brief Scale a color with a scalar
         *
         */

        Color Color::operator*(const double p) const
        {
            assert(p >= 0.0 && p <= 1.0);

            Color return_color(p*r,p*g,p*b);

            return return_color;
        }

        /**
         * @brief Get the amount of red in a color
         */

        double Color::get_red()   const {return r;}

        /**
         * @brief Get the amount of green in a color
         */

        double Color::get_green() const {return g;}

        /**
         * @brief Get the amount of blue in a color
         */

        double Color::get_blue()  const {return b;};

        /**
         * @brief Get a pixel with this color
         * */

        Pixel Color::get_pixel() const {return pixel;}


        /**
         * @brief Output a pixel with this color
         */

        std::ostream& operator<<(std::ostream& stream, Color& output_color)
        {
            Pixel output_pixel = output_color.get_pixel();
            stream.write(reinterpret_cast<char*>(&output_pixel),sizeof(Pixel));
            return stream;
        }

        /**
         * @brief Create a linear color map
         * @param v_min *double* minimum value to be represented on the color map
         * @param v_max *double* maximum value to be represented on the color map
         * @param color_min *Color* Color at the bottom of the scale
         * @param color_max *Color* Color at the top of the scale
         * @param number_of_levels *int* number of levels between v_min and v_max
         * @return linear_color_map *std::map<double,Color>* a map associating double values in the interval
         * \f$ [v_{min},v_{max}] \f$ to a color linearly interpolated between *color_min* and *color_max*
         *
         * This function creates a linear color map associating the color *color_min* to *v_min* and the color
         * *color_max* to *v_max*, dividing the interval \f$[v_{min},v_{max}]\f$ into *number_of_levels* sub-intervals and
         * associating a color in each sub-interval by interpolating colors linearly between *color_min* and
         * *color_max*.
         *
         */

        std::map<double,Color> create_linear_color_map(double v_min,
                                                       double v_max,
                                                       Color color_min,
                                                       Color color_max,
                                                       int number_of_levels)
        {

            std::map<double,Color> linear_color_map;
            double delta_v=(v_max-v_min)/number_of_levels;
            Color color;
            for(int i = 0; i <= number_of_levels; ++i)
            {
                double p = double(i)/number_of_levels;
                color = color_min*(1-p) + color_max*p;
                linear_color_map.insert(std::make_pair(v_min+i*delta_v,color));
            }

            return linear_color_map;

        }

        /**
          * @brief Create a diverging color map
          * @param v_min *double* minimum value to be represented on the color map
          * @param v_mid *double* middle value to be represented on the color map
          * @param v_max *double* maximum value to be represented on the color map
          * @param color_min *Color* Color at the bottom of the scale
          * @param color_mid *Color* Color at the middle of the scale
          * @param color_max *Color* Color at the top of the scale
          * @param number_of_levels *int* number of levels between v_min and v_max
          * @return diverging_color_map *std::map<double,Color>* a map associating double values in the interval
          * \f$ [v_{min},v_{mid}] \f$ to a color linearly interpolated between *color_min* and *color_mid*, and double values
          * in the interval \f$ [v_{mid},v_{max}] \f$ to a color linearly interpolated between *color_mid* and *color_max*.
          *
          */

        std::map<double,Color> create_diverging_color_map(double v_min,
                                                          double v_mid,
                                                          double v_max,
                                                          Color color_min,
                                                          Color color_mid,
                                                          Color color_max,
                                                          int number_of_levels)
        {

            std::map<double,Color> diverging_color_map;

            Color color;

            for(int i = 0; i <= number_of_levels; ++i)
            {
                double p = double(i)/number_of_levels;
                double delta_v;
                if (p <= 0.5)
                {
                    color = color_mid * p * 2 + color_min * (0.5 - p) * 2;
                    delta_v = std::abs(v_mid-v_min)/number_of_levels;
                }
                else
                {
                    color = color_max * (p - 0.5) * 2.0 + color_mid * (1.0 - p) * 2.0;
                    delta_v = std::abs(v_max-v_mid)/number_of_levels;
                }
                diverging_color_map.insert(std::make_pair(v_min+2*i*delta_v,color));
            }

            return diverging_color_map;

        }

        /**
         * @brief Write a 2-dimensional double Vector directly as a ppm image
         * */

        void write_ppm(const char* filename,
                       Vector<double> output_vector,
                       std::map<double,Color> color_map)
        {

            assert(output_vector.order() == 2);

            std::ofstream ofs(filename, std::ios_base::out | std::ios_base::binary);
            ofs << "P6" << std::endl
                << output_vector.extent(0) << ' '
                << output_vector.extent(1) << std::endl
                << "255" << std::endl;

            for (int i = 0; i < output_vector.extent(0); ++i)
                for (int j = 0; j < output_vector.extent(1); ++j)
                    ofs << color_map.lower_bound(output_vector(i,j))->second;

            ofs.close();
        }

    }
}