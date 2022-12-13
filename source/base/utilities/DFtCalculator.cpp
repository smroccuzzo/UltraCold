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

#include "DFtCalculator.hpp"
#include <iostream>
#include <cstring>
#include <complex>

namespace UltraCold{

    namespace MKLWrappers
    {

        /**
         * @brief Constructor for complex-complex transforms
         * @param forward_domain_vector *Vector<std::complex<double>>* Vector defined on the forward domain, i.e.
         * whose domain lives in real space
         * @param backward_domain_vector *Vector<std::complex<double>>* Vector defined on the backward domain, i.e.
         * whose domain lives in Fourier space
         * 
         * The constructor initializes the appropriate descriptor type.
         * 
         */

        DFtCalculator::DFtCalculator(Vector<std::complex<double>>& forward_domain_vector,
                                     Vector<std::complex<double>>& backward_domain_vector)
        {

            // Initialize the pointers

            forward_domain_pointer  = (std::complex<double>*) forward_domain_vector.data();
            backward_domain_pointer = (std::complex<double>*) backward_domain_vector.data();

            // Get the extents of the Vectors 

            forward_domain_n1=forward_domain_vector.extent(0);
            forward_domain_n2=forward_domain_vector.extent(1);
            forward_domain_n3=forward_domain_vector.extent(2);

            backward_domain_n1=backward_domain_vector.extent(0);
            backward_domain_n2=backward_domain_vector.extent(1);
            backward_domain_n3=backward_domain_vector.extent(2);

            // Check that the dimensions of the Vectors are equal, and rise an error otherwise
            assert(forward_domain_vector.order() == backward_domain_vector.order());
            assert(forward_domain_n1 == backward_domain_n1);
            assert(forward_domain_n2 == backward_domain_n2);
            assert(forward_domain_n3 == backward_domain_n3);

            // Initialize descriptor...

            // ... in 1D

            if(forward_domain_vector.order() == 1)
            {
                status = DftiCreateDescriptor(&dft_descriptor,DFTI_DOUBLE,DFTI_COMPLEX,1,forward_domain_n1);
                if(status != 0)
                {
                    char* error_message = DftiErrorMessage(status);
                    std::cout
                        << "\n\n"
                        << "**********************************************************************************\n"
                        << "Error found in constructor of a DFtCalculator for a complex-complex 1D transform. \n"
                        << "Failed to create descriptor. MKL is raising the following error message: \n"
                        << error_message
                        << "\nTerminating the program now...\n"
                        << "**********************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
            }
            else if(forward_domain_vector.order()==2)
            {

                // ... in 2D

                long N[2] = {forward_domain_n1,forward_domain_n2};
                status = DftiCreateDescriptor(&dft_descriptor,DFTI_DOUBLE,DFTI_COMPLEX,2,N);
                if(status != 0)
                {
                    char* error_message = DftiErrorMessage(status);
                    std::cout
                        << "\n\n"
                        << "**********************************************************************************\n"
                        << "Error found in constructor of a DFtCalculator for a complex-complex 2D transform. \n"
                        << "Failed to create descriptor. MKL is raising the following error message: \n"
                        << error_message
                        << "\nTerminating the program now...\n"
                        << "**********************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
            }
            else
            {
                // ... in 3D

                long N[3] = {forward_domain_n1,forward_domain_n2,forward_domain_n3};
                status = DftiCreateDescriptor(&dft_descriptor,DFTI_DOUBLE,DFTI_COMPLEX,3,N);
                if(status != 0)
                {
                    char* error_message = DftiErrorMessage(status);
                    std::cout
                        << "\n\n"
                        << "**********************************************************************************\n"
                        << "Error found in constructor of a DFtCalculator for a complex-complex 3D transform. \n"
                        << "Failed to create descriptor. MKL is raising the following error message: \n"
                        << error_message
                        << "\nTerminating the program now...\n"
                        << "**********************************************************************************\n"
                        << "\n\n"
                        <<
                    std::endl;
                    exit(1);
                }
            }

            // Specify that the transform is out-of-place.

            status = DftiSetValue(dft_descriptor,DFTI_PLACEMENT,DFTI_NOT_INPLACE);

            // Set the scale for the backward domain. The backward transform will thus be normalized

            status = DftiSetValue(dft_descriptor,DFTI_BACKWARD_SCALE,1./forward_domain_vector.size());

            // Commit descriptor. The object is now ready to be used in calculations

            status = DftiCommitDescriptor(dft_descriptor);

            if(status != 0)
            {
                char* error_message = DftiErrorMessage(status);
                std::cout
                    << "\n\n"
                    << "******************************************************************************\n"
                    << "Error found in constructor of a DFtCalculator for a complex-complex transform.\n"
                    << "Failed to commit the descriptor. MKL is raising the following error message: \n"
                    << error_message
                    << "\nTerminating the program now...\n"
                    << "******************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
        }

        /**
         * @brief Constructor for real-complex transforms.
         * @param forward_domain_vector *Vector<double>* Vector defined on the forward domain, i.e.
         * whose domain lives in real space
         * @param backward_domain_vector *Vector<std::complex<double>>* Vector defined on the backward domain, i.e.
         * whose domain lives in Fourier space
         * 
         * The constructor initializes the appropriate descriptor type.
         *
         * \note Since the Fourier transform of a real function is conjugate symmetric, only half of the values needs
         * to be stored. Using Intel's MKL DFT functions, the halved dimension is the last one.
         * 
         */

        DFtCalculator::DFtCalculator(Vector<double>& forward_domain_vector,
                                     Vector<std::complex<double>>& backward_domain_vector)
        {

            // Set the transform type
            transform_is_real_complex = true;

            // Get the pointers
            forward_domain_pointer  = (double*) forward_domain_vector.data();
            backward_domain_pointer = (std::complex<double>*) backward_domain_vector.data();

            // Get the extents of the Vectors 
            forward_domain_n1=forward_domain_vector.extent(0);
            forward_domain_n2=forward_domain_vector.extent(1);
            forward_domain_n3=forward_domain_vector.extent(2);

            backward_domain_n1=backward_domain_vector.extent(0);
            backward_domain_n2=backward_domain_vector.extent(1);
            backward_domain_n3=backward_domain_vector.extent(2);

            // Check that the dimensions of the Vectors are consistent, and rise an error otherwise.
            // Notice that, for real-complex transform, the extent of the complex (backward_domain_vector) Vector 
            // along the last direction must be one-half (plus one) the extent of the real (forward_domain_vector) Vector 
            // along the same direction.
            if((forward_domain_n2     == 0                   &&
                forward_domain_n3     == 0                   && 
                forward_domain_n1/2+1 != backward_domain_n1 )   
              ||
               (forward_domain_n3     == 0                   && 
                forward_domain_n2     != 0                   &&
                forward_domain_n1     != backward_domain_n1)     
              ||
               (forward_domain_n3     == 0                   && 
                forward_domain_n2     != 0                   &&
                forward_domain_n2/2+1 != backward_domain_n2)     
              || 
               (forward_domain_n3     != 0                   && 
                forward_domain_n2     != backward_domain_n2)         
              ||
               (forward_domain_n3     != 0                   && 
                forward_domain_n3/2+1 != backward_domain_n3)
            )
            {
                std::cout
                    << "\n\n"
                    << "***********************************************************************************************\n"
                    << "Error found in constructor of a DFtCalculator for a real-complex transform.\n"
                    << "The forward_domain_vector and backward_domain_vector Vectors do not have consistent dimensions.\n"
                    << "Terminating program now...\n"
                    << "***********************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
            else
            {

                // Initialize descriptor...

                // ... in 1D
                if(forward_domain_n2 == 0 && forward_domain_n3 == 0)
                {
                    
                    status = DftiCreateDescriptor(&dft_descriptor,DFTI_DOUBLE,DFTI_REAL,1,forward_domain_n1);
                    if(status != 0)
                    {
                        char* error_message = DftiErrorMessage(status);
                        std::cout
                            << "\n\n"
                            << "******************************************************************************\n"
                            << "Error found in constructor of a DFtCalculator for a real-complex 1D transform.\n"
                            << "Failed to create descriptor. MKL is raising the following error message: \n"
                            << error_message
                            << "\nTerminating the program now...\n"
                            << "*******************************************************************************\n"
                            << "\n\n"
                            <<
                        std::endl;
                        exit(1);
                    }
                }
                else if(forward_domain_n3 == 0)
                {

                    // ... in 2D
                    transform_is_2d = true;
                    long N[2] = {forward_domain_n1,forward_domain_n2};
                    status = DftiCreateDescriptor(&dft_descriptor,DFTI_DOUBLE,DFTI_REAL,2,N);
                    if(status != 0)
                    {
                        char* error_message = DftiErrorMessage(status);
                        std::cout
                            << "\n\n"
                            << "*************************************************************************************\n"
                            << "Error found in constructor of a DFtCalculator for a real-complex 2D transform.\n"
                            << "Failed to create descriptor. MKL is raising the following error message: \n"
                            << error_message
                            << "\nTerminating the program now...\n"
                            << "*************************************************************************************\n"
                            << "\n\n"
                            <<
                        std::endl;
                        exit(1);
                    }

                    // Set the strides for the forward and backward domains
                    strides_forward_domain_2d[0] = 0;
                    strides_forward_domain_2d[1] = forward_domain_n2;
                    strides_forward_domain_2d[2] = 1;
                    strides_backward_domain_2d[0] = 0;
                    strides_backward_domain_2d[1] = forward_domain_n2/2+1;
                    strides_backward_domain_2d[2] = 1;

                }
                else
                {
                    // ... in 3D
                    transform_is_3d = true;

                    long N[3] = {forward_domain_n1,forward_domain_n2,forward_domain_n3};
                    status = DftiCreateDescriptor(&dft_descriptor,DFTI_DOUBLE,DFTI_REAL,3,N);

                    if(status != 0)
                    {
                        char* error_message = DftiErrorMessage(status);
                        std::cout
                            << "\n\n"
                            << "*************************************************************************************\n"
                            << "Error found in constructor of a DFtCalculator for a real-complex 3D transform.\n"
                            << "Failed to create descriptor. MKL is raising the following error message: \n"
                            << error_message
                            << "\nTerminating the program now...\n"
                            << "*************************************************************************************\n"
                            << "\n\n"
                            <<
                        std::endl;
                        exit(1);
                    }
                     
                    // Set the strides
                    strides_forward_domain_3d[0] = 0;
                    strides_forward_domain_3d[1] = forward_domain_n2*forward_domain_n3;
                    strides_forward_domain_3d[2] = forward_domain_n3;
                    strides_forward_domain_3d[3] = 1;
                    strides_backward_domain_3d[0] = 0;
                    strides_backward_domain_3d[1] = forward_domain_n2*(forward_domain_n3/2+1);
                    strides_backward_domain_3d[2] = forward_domain_n3/2+1;
                    strides_backward_domain_3d[3] = 1;
                }
            }

            // Specify that the transform is out-of-place.
            status = DftiSetValue(dft_descriptor,DFTI_PLACEMENT,DFTI_NOT_INPLACE);

            // Specify the memory layout for real-complex transform
            status = DftiSetValue(dft_descriptor,DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX);

            // Set the scale for the backward domain. The backward transform will thus be normalized
            status = DftiSetValue(dft_descriptor,DFTI_BACKWARD_SCALE,1./forward_domain_vector.size());

            // Commit descriptor. The object is now ready to be used in calculations
            status = DftiCommitDescriptor(dft_descriptor);

            if(status != 0)
            {
                char* error_message = DftiErrorMessage(status);
                std::cout
                    << "\n\n"
                    << "*************************************************************************************\n"
                    << "Error found in constructor of a DFtCalculator for a real-complex 3D transform.\n"
                    << "Failed to commit the descriptor. MKL is raising the following error message: \n"
                    << error_message
                    << "\nTerminating the program now...\n"
                    << "*************************************************************************************\n"
                    << "\n\n"
                    <<
                std::endl;
                exit(1);
            }
        
        }

        /**
         * @brief Destructor, free the descriptor handler
         * 
         */

        DFtCalculator::~DFtCalculator()
        {
            DftiFreeDescriptor(&dft_descriptor);
        }

        /**
         * @brief Reinit a complex-to-complex transform
         * */

        void DFtCalculator::reinit(Vector<std::complex<double>> &forward_domain_vector,
                                   Vector<std::complex<double>> &backward_domain_vector)
        {

            // Initialize the pointers
            forward_domain_pointer  = (std::complex<double>*) forward_domain_vector.data();
            backward_domain_pointer = (std::complex<double>*) backward_domain_vector.data();

            // Get the extents of the Vectors
            forward_domain_n1=forward_domain_vector.extent(0);
            forward_domain_n2=forward_domain_vector.extent(1);
            forward_domain_n3=forward_domain_vector.extent(2);

            backward_domain_n1=backward_domain_vector.extent(0);
            backward_domain_n2=backward_domain_vector.extent(1);
            backward_domain_n3=backward_domain_vector.extent(2);

            // Check that the dimensions of the Vectors are equal, and rise an error otherwise
            if(forward_domain_n1 != backward_domain_n1 ||
               forward_domain_n2 != backward_domain_n2 ||
               forward_domain_n3 != backward_domain_n3)
            {
                std::cout
                        << "\n\n"
                        << "*********************************************************************************************\n"
                        << "Error found in re-initializer of a DFtCalculator for a complex-complex transform.\n"
                        << "The forward_domain_vector and backward_domain_vector Vectors do not have the same dimensions.\n"
                        << "Terminating program now...\n"
                        << "*********************************************************************************************\n"
                        << "\n\n"
                        <<
                        std::endl;
                exit(1);
            }
            else
            {

                // Initialize descriptor...

                // ... in 1D
                if(forward_domain_n2 == 0 && forward_domain_n3 == 0)
                {
                    status = DftiCreateDescriptor(&dft_descriptor,DFTI_DOUBLE,DFTI_COMPLEX,1,forward_domain_n1);
                    if(status != 0)
                    {
                        char* error_message = DftiErrorMessage(status);
                        std::cout
                                << "\n\n"
                                << "**********************************************************************************\n"
                                << "Error found in re-initializer of a DFtCalculator for a complex-complex 1D transform. \n"
                                << "Failed to create descriptor. MKL is raising the following error message: \n"
                                << error_message
                                << "\nTerminating the program now...\n"
                                << "**********************************************************************************\n"
                                << "\n\n"
                                <<
                                std::endl;
                        exit(1);
                    }
                }
                else if(forward_domain_n3 == 0)
                {

                    // ... in 2D
                    long N[2] = {forward_domain_n1,forward_domain_n2};
                    status = DftiCreateDescriptor(&dft_descriptor,DFTI_DOUBLE,DFTI_COMPLEX,2,N);
                    if(status != 0)
                    {
                        char* error_message = DftiErrorMessage(status);
                        std::cout
                                << "\n\n"
                                << "**********************************************************************************\n"
                                << "Error found in re-initializer of a DFtCalculator for a complex-complex 2D transform. \n"
                                << "Failed to create descriptor. MKL is raising the following error message: \n"
                                << error_message
                                << "\nTerminating the program now...\n"
                                << "**********************************************************************************\n"
                                << "\n\n"
                                <<
                                std::endl;
                        exit(1);
                    }
                }
                else
                {
                    // ... in 3D
                    long N[3] = {forward_domain_n1,forward_domain_n2,forward_domain_n3};
                    status = DftiCreateDescriptor(&dft_descriptor,DFTI_DOUBLE,DFTI_COMPLEX,3,N);
                    if(status != 0)
                    {
                        char* error_message = DftiErrorMessage(status);
                        std::cout
                                << "\n\n"
                                << "**********************************************************************************\n"
                                << "Error found in re-initializer of a DFtCalculator for a complex-complex 3D transform. \n"
                                << "Failed to create descriptor. MKL is raising the following error message: \n"
                                << error_message
                                << "\nTerminating the program now...\n"
                                << "**********************************************************************************\n"
                                << "\n\n"
                                <<
                                std::endl;
                        exit(1);
                    }
                }
            }

            // Specify that the transform is out-of-place.
            status = DftiSetValue(dft_descriptor,DFTI_PLACEMENT,DFTI_NOT_INPLACE);

            // Set the scale for the backward domain. The backward transform will thus be normalized
            status = DftiSetValue(dft_descriptor,DFTI_BACKWARD_SCALE,1./forward_domain_vector.size());

            // Commit descriptor. The object is now ready to be used in calculations
            status = DftiCommitDescriptor(dft_descriptor);

            if(status != 0)
            {
                char* error_message = DftiErrorMessage(status);
                std::cout
                        << "\n\n"
                        << "******************************************************************************\n"
                        << "Error found in re-initializer of a DFtCalculator for a complex-complex transform.\n"
                        << "Failed to commit the descriptor. MKL is raising the following error message: \n"
                        << error_message
                        << "\nTerminating the program now...\n"
                        << "******************************************************************************\n"
                        << "\n\n"
                        <<
                        std::endl;
                exit(1);
            }

        }

        /**
         * @brief Re-initialize a real-to-complex transform
         * */

        void DFtCalculator::reinit(Vector<double> &forward_domain_vector,
                                   Vector<std::complex<double>> &backward_domain_vector)
        {

            // Set the transform type
            transform_is_real_complex = true;

            // Get the pointers
            forward_domain_pointer  = (double*) forward_domain_vector.data();
            backward_domain_pointer = (std::complex<double>*) backward_domain_vector.data();

            // Get the extents of the Vectors
            forward_domain_n1=forward_domain_vector.extent(0);
            forward_domain_n2=forward_domain_vector.extent(1);
            forward_domain_n3=forward_domain_vector.extent(2);

            backward_domain_n1=backward_domain_vector.extent(0);
            backward_domain_n2=backward_domain_vector.extent(1);
            backward_domain_n3=backward_domain_vector.extent(2);

            // Check that the dimensions of the Vectors are consistent, and rise an error otherwise.
            // Notice that, for real-complex transform, the extent of the complex (backward_domain_vector) Vector
            // along the last direction must be one-half (plus one) the extent of the real (forward_domain_vector) Vector
            // along the same direction.
            if((forward_domain_n2     == 0                   &&
                forward_domain_n3     == 0                   &&
                forward_domain_n1/2+1 != backward_domain_n1 )
               ||
               (forward_domain_n3     == 0                   &&
                forward_domain_n2     != 0                   &&
                forward_domain_n1     != backward_domain_n1)
               ||
               (forward_domain_n3     == 0                   &&
                forward_domain_n2     != 0                   &&
                forward_domain_n2/2+1 != backward_domain_n2)
               ||
               (forward_domain_n3     != 0                   &&
                forward_domain_n2     != backward_domain_n2)
               ||
               (forward_domain_n3     != 0                   &&
                forward_domain_n3/2+1 != backward_domain_n3)
                    )
            {
                std::cout
                        << "\n\n"
                        << "***********************************************************************************************\n"
                        << "Error found in re-initializer of a DFtCalculator for a real-complex transform.\n"
                        << "The forward_domain_vector and backward_domain_vector Vectors do not have consistent dimensions.\n"
                        << "Terminating program now...\n"
                        << "***********************************************************************************************\n"
                        << "\n\n"
                        <<
                        std::endl;
                exit(1);
            }
            else
            {

                // Initialize descriptor...

                // ... in 1D
                if(forward_domain_n2 == 0 && forward_domain_n3 == 0)
                {

                    status = DftiCreateDescriptor(&dft_descriptor,DFTI_DOUBLE,DFTI_REAL,1,forward_domain_n1);
                    if(status != 0)
                    {
                        char* error_message = DftiErrorMessage(status);
                        std::cout
                                << "\n\n"
                                << "******************************************************************************\n"
                                << "Error found in re-initializer of a DFtCalculator for a real-complex 1D transform.\n"
                                << "Failed to create descriptor. MKL is raising the following error message: \n"
                                << error_message
                                << "\nTerminating the program now...\n"
                                << "*******************************************************************************\n"
                                << "\n\n"
                                <<
                                std::endl;
                        exit(1);
                    }
                }
                else if(forward_domain_n3 == 0)
                {

                    // ... in 2D
                    transform_is_2d = true;
                    long N[2] = {forward_domain_n1,forward_domain_n2};
                    status = DftiCreateDescriptor(&dft_descriptor,DFTI_DOUBLE,DFTI_REAL,2,N);
                    if(status != 0)
                    {
                        char* error_message = DftiErrorMessage(status);
                        std::cout
                                << "\n\n"
                                << "*************************************************************************************\n"
                                << "Error found in re-initializer of a DFtCalculator for a real-complex 2D transform.\n"
                                << "Failed to create descriptor. MKL is raising the following error message: \n"
                                << error_message
                                << "\nTerminating the program now...\n"
                                << "*************************************************************************************\n"
                                << "\n\n"
                                <<
                                std::endl;
                        exit(1);
                    }

                    // Set the strides for the forward and backward domains
                    strides_forward_domain_2d[0] = 0;
                    strides_forward_domain_2d[1] = forward_domain_n2;
                    strides_forward_domain_2d[2] = 1;
                    strides_backward_domain_2d[0] = 0;
                    strides_backward_domain_2d[1] = forward_domain_n2/2+1;
                    strides_backward_domain_2d[2] = 1;

                }
                else
                {
                    // ... in 3D
                    transform_is_3d = true;

                    long N[3] = {forward_domain_n1,forward_domain_n2,forward_domain_n3};
                    status = DftiCreateDescriptor(&dft_descriptor,DFTI_DOUBLE,DFTI_REAL,3,N);

                    if(status != 0)
                    {
                        char* error_message = DftiErrorMessage(status);
                        std::cout
                                << "\n\n"
                                << "*************************************************************************************\n"
                                << "Error found in re-initializer of a DFtCalculator for a real-complex 3D transform.\n"
                                << "Failed to create descriptor. MKL is raising the following error message: \n"
                                << error_message
                                << "\nTerminating the program now...\n"
                                << "*************************************************************************************\n"
                                << "\n\n"
                                <<
                                std::endl;
                        exit(1);
                    }

                    // Set the strides
                    strides_forward_domain_3d[0] = 0;
                    strides_forward_domain_3d[1] = forward_domain_n2*forward_domain_n3;
                    strides_forward_domain_3d[2] = forward_domain_n3;
                    strides_forward_domain_3d[3] = 1;
                    strides_backward_domain_3d[0] = 0;
                    strides_backward_domain_3d[1] = forward_domain_n2*(forward_domain_n3/2+1);
                    strides_backward_domain_3d[2] = forward_domain_n3/2+1;
                    strides_backward_domain_3d[3] = 1;
                }
            }

            // Specify that the transform is out-of-place.
            status = DftiSetValue(dft_descriptor,DFTI_PLACEMENT,DFTI_NOT_INPLACE);

            // Specify the memory layout for real-complex transform
            status = DftiSetValue(dft_descriptor,DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX);

            // Set the scale for the backward domain. The backward transform will thus be normalized
            status = DftiSetValue(dft_descriptor,DFTI_BACKWARD_SCALE,1./forward_domain_vector.size());

            // Commit descriptor. The object is now ready to be used in calculations
            status = DftiCommitDescriptor(dft_descriptor);

            if(status != 0)
            {
                char* error_message = DftiErrorMessage(status);
                std::cout
                        << "\n\n"
                        << "*************************************************************************************\n"
                        << "Error found in re-initializer of a DFtCalculator for a real-complex 3D transform.\n"
                        << "Failed to commit the descriptor. MKL is raising the following error message: \n"
                        << error_message
                        << "\nTerminating the program now...\n"
                        << "*************************************************************************************\n"
                        << "\n\n"
                        <<
                        std::endl;
                exit(1);
            }

        }

        /**
         * @brief Calculate a forward transform
         * 
         */

        void DFtCalculator::compute_forward()
        {

            // Set the strides correctly for real-complex transforms

            if(transform_is_real_complex)
            {

                // Set the strides  

                if(transform_is_2d)
                {
                    status = DftiSetValue(dft_descriptor, DFTI_INPUT_STRIDES,  strides_forward_domain_2d);
                    status = DftiSetValue(dft_descriptor, DFTI_OUTPUT_STRIDES, strides_backward_domain_2d);    
                }
                else if(transform_is_3d)
                {
                    status = DftiSetValue(dft_descriptor, DFTI_INPUT_STRIDES,  strides_forward_domain_3d);
                    status = DftiSetValue(dft_descriptor, DFTI_OUTPUT_STRIDES, strides_backward_domain_3d);  
                }

                status = DftiCommitDescriptor(dft_descriptor);

            }
            
            status = DftiComputeForward(dft_descriptor,forward_domain_pointer,backward_domain_pointer);

        }


        /**
         * @brief Calculate a backward transform. The transform is already normalized.
         * 
         */

        void DFtCalculator::compute_backward()
        {

            // Set the strides correctly for real-complex transforms

            if(transform_is_real_complex)
            {

                // Set the strides  

                if(transform_is_2d)
                {
                    status = DftiSetValue(dft_descriptor, DFTI_OUTPUT_STRIDES, strides_forward_domain_2d);
                    status = DftiSetValue(dft_descriptor, DFTI_INPUT_STRIDES, strides_backward_domain_2d);
                }
                else if(transform_is_3d)
                {
                    status = DftiSetValue(dft_descriptor, DFTI_OUTPUT_STRIDES, strides_forward_domain_3d);
                    status = DftiSetValue(dft_descriptor, DFTI_INPUT_STRIDES, strides_backward_domain_3d);
                }

                status = DftiCommitDescriptor(dft_descriptor);

            }

            status = DftiComputeBackward(dft_descriptor,backward_domain_pointer,forward_domain_pointer);

        }

    } // namespace MKLWrappers
} // namespace UltraCold