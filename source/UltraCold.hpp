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

#ifndef ULTRACOLD
#define ULTRACOLD

// All the includes specific to UltraCold

#include "DataWriter.hpp"
#include "InputParser.hpp"
#include "Vector.hpp"
#include "HardwareInspector.hpp"
#include "DFtCalculator.hpp"
#include "GPSolvers.hpp"
#include "BogolyubovSolvers.hpp"
#include "Colors.hpp"

#ifdef ULTRACOLD_WITH_CUDA
#include "cudaGPSolver.cuh"
#include "cudaDipolarGPSolver.cuh"
#endif

// A few useful macros, defining constants

#ifndef ULTRACOLD_PI
#define PI 3.1415926535897932
#define TWOPI (2*PI)
#endif

// Definition of namespace UltraCold

/**
 * @brief All the classes and functions necessary to work with UltraCold
 * 
 */

namespace UltraCold{};

#endif // ULTRACOLD