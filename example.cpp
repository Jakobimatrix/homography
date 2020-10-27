#include <iostream>

#include "include/homography/homography.hpp"


int main() {
  const Eigen::Vector2d image_size(300, 300);
  std::array<Eigen::Vector2d, 4> A;
  std::array<Eigen::Vector2d, 4> B;
  std::array<EigenSTL::pair_Vector2d_Vector2d, 4> associated_points;
  Eigen::Matrix3d H;

  // Some points in world...
  A[0] << -2, 2;   // A
  A[1] << 2, 2;    // B
  A[2] << 2, -2;   // C
  A[3] << -2, -2;  // D

  // ...should match these picture coordinates
  B[0] << 0, 0;               // A
  B[1] << image_size.x(), 0;  // B
  B[2] << image_size;         // C
  B[3] << 0, image_size.y();  // D

  for (unsigned int i = 0; i < 4; i++) {
    associated_points[i] = std::make_pair(A[i], B[i]);
  }

  // get the transformation matrix
  H = tool::findHomography<tool::SystemSolverMethode::PARTIAL_PIV_LU>(associated_points);

  std::cout << "The homography matrix for the sample points is:\n"
            << H << std::endl;
  // Transform a picture point to the corresponding world point.
  const Eigen::Vector2d a_picture_point(34, 57);
  const Eigen::Vector2d a_world_point =
      tool::computeTransformation(H.inverse(), a_picture_point);

  std::cout << "The Point in A (" << a_picture_point.transpose() << ") corresponds to the point ("
            << a_world_point.transpose() << ") in B." << std::endl;

  // Transform a world point to the corresponding picture point.
  const Eigen::Vector2d b_world_point(2.224, 3.335);
  const Eigen::Vector2d b_picture_point = tool::computeTransformation(H, b_world_point);
  std::cout << "The Point in A (" << b_picture_point.transpose() << ") corresponds to the point ("
            << b_world_point.transpose() << ") in B." << std::endl;

  return 0;
}
