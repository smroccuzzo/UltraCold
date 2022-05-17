# This script configure, compile and install packages bundled to UltraCold

# Configure, build and install arpack-ng

cd arpack-ng
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$1 \
      -DCMAKE_Fortran_COMPILER=ifort \
      -DCMAKE_C_COMPILER=icc \
      -DCMAKE_CXX_COMPILER=icpc \
      -DICB=ON \
      -DMPI=ON \
      -DCMAKE_BUILD_TYPE=Release \
      ..

make -j 4
make test
make install
