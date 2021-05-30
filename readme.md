# BundleAdjustment
## Overview
This project is a version of visual bundle adjustment without using any optimized library. 

The project firstly generated n 3D points and m camera data (**as ground truth**) by simulation, and it is assumed that all of n 3D points can be observed by m cameras. Secondly, it uses n 3D points and m cameras to generate image feature points according to the camera model, and adds Gaussian noise to the image feature points (**as measurement data**). The variance of Gaussian noise can be specified by the user. And then, the measurement data (image feature points) is used to recover the camera pose and 3D  points according to the Structure-from-Motion. Finally, the reprojection error model is constructed, and the LM algorithm is implemented to optimize the reprojection error.

## Dependency
* OpenCV3: **Required at least 3.0**.
* Eigen3: **Required at least 3.0**.
* Sophus
## Code Structure
* `Camera.h/.cpp`: a class which describe the camera model.
* `MapPoint.h/.cpp` : a class which describe the 3d landmark.
* `CostFunction.h/.cpp`: a class which describe the re-projection cost function.
* `Sfm.h/.cpp` : an initialization process for optimization.   
* `BundleAdjustment.h/.cpp`: a Levenberg-Marquardt optimization method is adopted.
* `testBA.cpp`: main function.

## Compile
```
mkdir build
cd build
cmake ..
make
```

## Reference
[1] Triggs B, McLauchlan P F, Hartley R I, et al. Bundle adjustment—a modern synthesis[C]//International workshop on vision algorithms. Springer, Berlin, Heidelberg, 1999: 298-372.

[2] https://zhuanlan.zhihu.com/p/64471565
