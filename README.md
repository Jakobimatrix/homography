# homography
Calculates the 3x3 homography matrix from a given set of corresponding point pairs to transform between two points of view.

Please use clang-tidy if you want to contribute: [easy installation](https://github.com/Jakobimatrix/initRepro)
### Dependencies:
* Eigen 3.3

### Features
* calculate the  3x3 homography matrix useing different solver for linear Ax=b Problems provided by Eigen
* Calculates the ''exact'' homography matrix  when given 4 pairs of points
* Estimates the homography matrix if given more than 4 pairs of points e.g. if the points on one plane are measurements with noise.
* example usage is provided in include/homography/homography.hpp

### Installation
* Install Eigen
* If you use CMake to build your project place this repro inside 'path' and use **add_subdirectory(path/
homography)** and **target_link_libraries(your_lib homography_lib)**
* Or just include the header
