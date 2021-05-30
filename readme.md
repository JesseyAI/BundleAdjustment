# BundleAdjustment
## Overview
This project is a version of visual bundle adjustment. 

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
<<<<<<< HEAD

# Compile
```
mkdir build
cd build
cmake ..
make
```
## Result
![](![t](https://github.com/JesseyAI/BundleAdjustment/blob/main/image/Accuracy(a).eps))
## Reference

=======

# Compile
```
mkdir build
cd build
cmake ..
make
```
## Result
![](![t](https://github.com/JesseyAI/BundleAdjustment/blob/main/image/Accuracy(a).eps))
## Reference
>>>>>>> 55c98f725e7ddf94838ca135fe7e822c51cc47a6