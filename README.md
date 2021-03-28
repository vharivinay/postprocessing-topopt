# Post-processing of 2D topology optimized designs and extraction of smooth and crisp geometry
Code to post-process topology optimized 2D results to a usable smooth geometry. TO solver used is 88-line code by DTU, Denmark.

Paper: [Efficient topology optimization in MATLAB using 88 lines of code](https://link.springer.com/article/10.1007%2Fs00158-010-0594-7)

Original TO solver: [top88.m](https://www.topopt.mek.dtu.dk/Apps-and-software/Efficient-topology-optimization-in-MATLAB)

### Abstract
Topology optimization has seen rapid developments in it's field with algorithms getting better and faster all the time. These new algorithms help reduce the lead time from concept development to a finished product. Simulation and post-processing of geometry is one of the major developmental costs. Post-processing of this geometry also takes up a lot of time and is dependent of the quality of the geometry output from the solver to make the product ready for rapid prototyping or final production. 

The work done in this thesis deals with the post-processing of the results obtained from a topology optimization algorithms which outputs the result as a 2D image. A suitable methodology is discussed where this image is processed and converted into an CAD geometry all while minimizing deviation in geometry, compliance and volume fraction. Further on, a validation of the designs is performed as to measure the deviation of the extracted geometry from the post-processed result. 

The workflow is coded using MATLAB and uses an image based post-processing approach. The proposed workflow is tested on several numerical examples to assess the performance, limitations and numerical instabilities.

### What this code does
#### Unsmooth topology optimized design
<img width="600px" height="200px" src="/mbb_unsmooth.png">

#### Smooth and crisp geometry
<img width="600px" height="200px" src="/mbb_smooth.png">

#### Extracted curves
<img width="600px" height="200px" src="/mbb_curves.png">

#### Extruded 3D Geometry
<img width="600px" height="300px" src="/mbb_3dgeom.png">
