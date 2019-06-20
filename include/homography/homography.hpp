/*
 * \file: homography.h
 * \brief: Find a find a homography -transformation- matrix to transform between two points of view.
 * \date: 10.11.2017
 * \author: Jakob Wandel
 * \version: 1.0
 * \dependencies: Eigen
 *
Example Usage:
  const Eigen::Vector2d image_size(300, 300);
  std::array<Eigen::Vector2d, 4> A;
  std::array<Eigen::Vector2d, 4> B;
  std::array<EigenSTL::pair_Vector2d_Vector2d, 4> associated_points;
  Eigen::Matrix3d H;

  // Some points in world...
  A[0] << -2, 2;  // A
  A[1] << 2, 2;   // B
  A[2] << 2, -2;  // C
  A[3] << -2, -2; // D

  // ...should match these picture coordinates
  B[0] << 0, 0;              // A
  B[1] << image_size.x(), 0; // B
  B[2] << image_size;        // C
  B[3] << 0, image_size.y(); // D

  for (unsigned int i = 0; i < 4; i++) {
    associated_points[i] = std::make_pair(A[i], B[i]);
  }

  //get the transformation matrix
  H = tool::Homography::findHomography(
      associated_points, tool::Homography::Methode::PARTIAL_PIV_LU);

  // Transform a picture point to the corresponding world point.
  const Eigen::Vector2d a_picture_point(34,57);
  const Eigen::Vector2d a_world_point = tool::Homography::computeTransformation(H.inverse(), a_picture_point);

  // Transform a world point to the corresponding picture point.
  const Eigen::Vector2d b_world_point(2.224,3.335);
  const Eigen::Vector2d b_picture_point = tool::Homography::computeTransformation(H, b_world_point);
*/

#ifndef HOMOGRAPHY_H
#define HOMOGRAPHY_H

#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Householder>
#include <eigen3/Eigen/Geometry>

namespace EigenSTL {
  typedef Eigen::Matrix<double, 8, 1> Vector8d;
  typedef Eigen::Matrix<double, 9, 1> Vector9d;
  typedef std::pair<Eigen::Vector2d,Eigen::Vector2d> pair_Vector2d_Vector2d;
  typedef std::vector<pair_Vector2d_Vector2d, Eigen::aligned_allocator<pair_Vector2d_Vector2d>> vector_pair_Vector2d_Vector2d;
}
namespace tool{

class Homography{
public:
  // clang-format off
  /*
   * Benchmark: https://eigen.tuxfamily.org/dox/group__DenseDecompositionBenchmark.html
   * Tabel origin: https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
   *
   * Decomposition                    Requirements  	Speed             Speed     Accuracy
   *                                  on the matrix   (small-to-medium) (large)
   * -----------------------------------------------------------------------------------
   * PartialPivLU (Cholesky)          Invertible      ++                  ++        +
   * FullPivLU    (Cholesky)          None            -                   --        +++
   * HouseholderQR                    None            ++                  ++        +
   * ColPivHouseholderQR              None            + 	                -         +++
   * FullPivHouseholderQR             None            -                   --        +++
   * CompleteOrthogonalDecomposition 	None            +                   -         +++
   * LLT                              Pos definite    +++                 +++       +
   * LDLT                             Pos or negative	+++                 +         ++
   *                                  semidef
   * BDCSVD                           None            -                   -         +++
   * JacobiSVD                        None            -                   ---       +++
  */
  // clang-format on

  enum Methode{
    FULL_PIV_LU,
    PARTIAL_PIV_LU,
    FULL_PIV_QR,
    PARTIAL_QR,
    QR,
    COMPLETE_ORTHOGONAL_DECOMPOSITION,
    LLT,
    LDLT
  };

  enum DynamicMethode{
    JACOBI_SCD,
    BDCSVD  //! warning: https://eigen.tuxfamily.org/dox/classEigen_1_1BDCSVD.html
  };

  // clang-format off
  /*    |h00  h01 h02|
   * H =|h10  h11 h12|
   *    |h20  h21 h22|
   *
   *                      |Px|
   * ~Px' = [h00 h01 h02]*|Py|
   *                      | 1|
   *
   *                      |Px|
   * ~Py' = [h10 h11 h12]*|Py|
   *                      | 1|
   *
   *                      |Px|
   *  Ph' = [h20 h21 h22]*|Py|
   *                      | 1|
   *
   *  Px' = ~Px'/Ph'
   *  Py' = ~Py'/Ph'
   *
   * -->
   *
   *                      |Px|                 |Px|
   *  Px' * [h20 h21 h22]*|Py| = [h00 h01 h02]*|Py|
   *                      | 1|                 | 1|
   *
   *                      |Px|                 |Px|
   *  Py' * [h20 h21 h22]*|Py| = [h10 h11 h12]*|Py|
   *                      | 1|                 | 1|
   *
   * --> (with a*b = b^T*a^T)
   *
   *              |h00|                 |h20|
   *  [Px Py  1]* |h01| - Px'[Px Py  1]*|h21|  = 0
   *              |h02|                 |h22|
   *
   *
   *              |h10|                 |h20|
   *  [Px Py  1]* |h11| - Py'[Px Py  1]*|h21|  = 0
   *              |h12|                 |h22|
   * -->
   *
   *
   *                                        |h00|
   *                                        |h01|
   * [Px Py  1  -Px'  -Px'Px  -Px'Py  -Px']*|h02| = 0
   *                                        |h20|
   *                                        |h21|
   *                                        |h22|
   *
   *                                        |h10|
   *                                        |h11|
   * [Px Py  1  -Px'  -Py'Px  -Py'Py  -Py']*|h12| = 0
   *                                        |h20|
   *                                        |h21|
   *                                        |h22|
   *
   * --> where P_i and P'_i is the ith- assosiated Point pair.
   *                                                                      |h00|
   *                                                                      |h01|
   *                                                                      |h02|
   *                                                                      |h10|
   * |Px_i  Py_i  1   0     0     0   -P'x_i*Px_i   -P'x_i*Py_i   -P'x_i|*|h11| = |0|
   * |0     0     0   Px_i  Py_i  1   -P'y_i*Px_i   -P'y_i*Py_i   -P'y_i| |h12|   |0|
   *                                                                      |h20|
   *                                                                      |h21|
   *                                                                      |h22|
   * =:P*h=0
   * With n (n>4) Points are given we get an overdetermined system of linear equations
   * which can be solved useing JACOBI_SCD or BDCSVD for bigger matrices.
   *
   * To avoide the trivial solution h=0 we add the constraint h22 = 1.
   */
  // clang-format on

  /*!
   * \brief Calculates the homography matrix (3x3) to convert between two points of view:
   * ~P_dst = H * P_src; useing at least 5 pairs of corresponding Points.
   * \param paired_points_src_dst The paired Points where first is source and second is aim.
   * \param methode Which dynamic methode should be used to calculate Ax=b. Have a look above which are avaiable.
   * \return Eigen::Matrix3d found homography.
   */
  static Eigen::Matrix3d findHomography(const EigenSTL::vector_pair_Vector2d_Vector2d& paired_points_src_dst, DynamicMethode methode){

    assert(paired_points_src_dst.size() > 4 && "tool::Homography::findHomography: There must be at least 5 pairs of points given.");

    Eigen::MatrixXd PH;
    PH.resize(paired_points_src_dst.size()*2 + 1,9);
    Eigen::VectorXd b(paired_points_src_dst.size()*2 + 1);

    for(unsigned int i = 0, j = 0; i < paired_points_src_dst.size(); i++ ){
      const double srcX = paired_points_src_dst[i].first.x();
      const double srcY = paired_points_src_dst[i].first.y();
      const double dstX = paired_points_src_dst[i].second.x();
      const double dstY = paired_points_src_dst[i].second.y();

      // clang-format off
      b(j) = 0.;
      PH.row(j++) << srcX, srcY,  1.,   0.,   0.,  0., -srcX*dstX, -srcY*dstX, -dstX;
      b(j) = 0.;
      PH.row(j++) <<   0.,   0.,  0., srcX, srcY,  1., -srcX*dstY, -srcY*dstY, -dstY;
      // clang-format on
    }
    // avoide the trivial solution h=0
    PH.row(paired_points_src_dst.size()*2) << 0., 0., 0., 0., 0., 0., 0., 0., 1;
    b(paired_points_src_dst.size()*2) = 1;

    // solv PH*x=b
    EigenSTL::Vector9d x;
    Eigen::Matrix3d result;

    switch (methode) {
    case DynamicMethode::JACOBI_SCD:
    {
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(PH, Eigen::ComputeThinU | Eigen::ComputeThinV);
      x = svd.solve(b);
      break;
    }
    case DynamicMethode::BDCSVD:
    {
      Eigen::BDCSVD<Eigen::MatrixXd> bdcsvd(PH, Eigen::ComputeThinU | Eigen::ComputeThinV);
      x = bdcsvd.solve(b);
      break;
    }
    }

    //const double relative_error = (PH*x - b).norm() / b.norm(); // norm() is L2 norm

    Eigen::Vector3d x1,x2,x3;
    x1 << x.head(3);
    x2 << x.segment(3,3);
    x3 << x.tail(3);

    //result << x; //doesnt work for reasons
    result << x1, x2, x3;
    return result.transpose();
  };


  // clang-format off
  /*
   * Useing the general solution from above:
   * Since P'h is just a factor which we will get rid of in the converstaion anyway
   * it is possible to restrict H by setting h22 = 1
   *
   * If we do this we can rewrite the general solution:
   *
   *                                                              |h00|
   *                                                              |h01|
   *                                                              |h02|
   * |Px_i  Py_i  1   0     0     0   -P'x_i*Px_i   -P'x_i*Py_i|* |h10| =  |P'x_i|
   * |0     0     0   Px_i  Py_i  1   -P'y_i*Px_i   -P'y_i*Py_i|  |h11|    |P'y_i|
   *                                                              |h12|
   *                                                              |h20|
   *                                                              |h21|
   *
   * With 4 Points we have a system of 8 linear equations with 8 unknown.
   * Now we can use use less expensive methodes to solve this system.
   * If there is no noise on the 4 given point pairs the solution will
   * be optimal apart from numerical problems.
   *
   * If you have no ground truth for the source / target Points but rather
   * measurements, consider useing the general methode to filter the noise.
   */
  // clang-format on

  /*!
   * \brief Calculates the homography matrix (3x3) to convert between two points of view:
   * ~P_dst = H * P_src; useing 4 pairs of corresponding Points.
   * \param paired_points_src_dst The 4 paired Points where first is source and second is aim.
   * \param methode Which methode should be used to calculate Ax=b. Have a look above which are avaiable.
   * \return Eigen::Matrix3d found homography.
   */
  typedef Eigen::Matrix<double, 8,8> HomograpyMatrix;
  static Eigen::Matrix3d findHomography(const std::array<EigenSTL::pair_Vector2d_Vector2d,4>& paired_points_src_dst, Methode methode){

      HomograpyMatrix PH;
      EigenSTL::Vector8d b;
      for(unsigned int i = 0, j = 0; i < 4; i++ ){

        const double srcX = paired_points_src_dst[i].first.x();
        const double srcY = paired_points_src_dst[i].first.y();
        const double dstX = paired_points_src_dst[i].second.x();
        const double dstY = paired_points_src_dst[i].second.y();

        // clang-format off
        b(j) = dstX;
        PH.row(j++) << srcX, srcY,  1.,   0.,   0.,  0., -srcX*dstX, -srcY*dstX;
        b(j) = dstY;
        PH.row(j++) <<   0.,   0.,  0., srcX, srcY,  1., -srcX*dstY, -srcY*dstY;
        // clang-format on
      }

      // solv PH*x=b
      EigenSTL::Vector8d x;
      Eigen::Matrix3d result;

      switch (methode) {
      case Methode::FULL_PIV_LU:
      {
        Eigen::FullPivLU<HomograpyMatrix> lu(PH);
        x = lu.solve(b);
        //x = PH.fullPivLu().solve(b);
        break;
      }
      case Methode::PARTIAL_PIV_LU:
      {
        Eigen::PartialPivLU<HomograpyMatrix> plu(PH);
        x = plu.solve(b);
        //x = PH.partialPivLu().solve(b);
        break;
      }
      case Methode::FULL_PIV_QR:
      {
        Eigen::FullPivHouseholderQR<HomograpyMatrix> fqr(PH);
        x = fqr.solve(b);
        //x = PH.fullPivHouseholderQr().solve(b);
        break;
      }
      case Methode::PARTIAL_QR:
      {
        Eigen::ColPivHouseholderQR<HomograpyMatrix> cqr(PH);
        x = cqr.solve(b);
        //x = PH.colPivHouseholderQr().solve(b);
        break;
      }
      case Methode::QR:
      {
        Eigen::HouseholderQR<HomograpyMatrix> qr(PH);
        x = qr.solve(b);
        //x = PH.householderQr().solve(b);
        break;
      }
      case Methode::LLT:
      {
        //Eigen::LLT<HomograpyMatrix> llt(PH);
        //x = llt.solve(b);
        x = PH.llt().solve(b);
        x = (PH.transpose() * PH).llt().solve(PH.transpose() * b);
        break;
      }
      case Methode::LDLT:
      {
        //Eigen::LDLT<HomograpyMatrix> ldlt(PH);
        //x = ldlt.solve(b);
        //x = PH.ldlt().solve(b);
        x = (PH.transpose() * PH).ldlt().solve(PH.transpose() * b);
        break;
      }
      case Methode::COMPLETE_ORTHOGONAL_DECOMPOSITION:
      {
        Eigen::CompleteOrthogonalDecomposition<HomograpyMatrix> cod(PH);
        x = cod.solve(b);
        //x = PH.completeOrthogonalDecomposition().solve(b);
        break;
      }
      }

      //const double relative_error = (PH*x - b).norm() / b.norm(); // norm() is L2 norm

      Eigen::Vector3d x1,x2,x3;
      x1 << x.head<3>();
      x2 << x.segment<3>(3);
      x3 << x.tail<2>(), 1.;

      //result << x, 1.; //doesnt work for reasons
      result << x1, x2, x3;
      return result.transpose();
  }

  static Eigen::Vector2d computeTransformation(const Eigen::Matrix3d& homography, const Eigen::Vector2d& src){
    const Eigen::Vector3d dst_homogeneous =
        homography * src.homogeneous();
    assert(dst_homogeneous(2) != 0.0 &&
           "A 3D point from the horizon was asked to be transformed!");

    return dst_homogeneous.hnormalized().head(2);
  };

};

}
#endif
