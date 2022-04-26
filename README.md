# MCRT_GPU
3D Voxel based Monte Carlo Radiation Transport on the GPU using OpenACC

Original MCRT is taken from [here](http://www-star.st-and.ac.uk/~kw25/research/montecarlo/montecarlo.html) and has been modernized from F7 to Modern Fortran.

GPU computation has been acheived by using OpenACC directives.

This repo is mainly me playing around with OpenACC and possibly other GPU frameworks/libraries (CUDA).

## Building
Requires Nvidia's NVFortran compiler, which can be downloaded from [here](https://developer.nvidia.com/hpc-sdk), as part of Nvidia's HPC SDK.

Can be built with make for [FPM](https://fpm.fortran-lang.org/en/index.html).

Make:

  ```
  cd src
  make
  cd ..
  ./mcgrid
  ```
FPM:

  `fpm @run`
  
When built successfully, program should print the average number of scatterings is ~57

The program swill also output the fluence to data/jmean.dat. This is a 201x201x201 float32 datacube that can easily be read in by Python or other languages.
