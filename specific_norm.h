/*! @file specific_norm.h
 *  @brief calculate some specific norm
 */
#include <Eigen/Dense>
#include <cmath>

/*! @brief Calculate difference bewteen two rotation matrices
 *  @details calculate the minimum angle to align the two matrices
 *  @return difference
 *  @param src rotation matrix 1
 *  @param src2 rotation matrix 2
 */
template<typename T>
T CalcDifferenceBetweenTwoRotationMatrices(const Eigen::Matrix<T, 3, 3, 0, 3, 3> &src, const Eigen::Matrix<T, 3, 3, 0, 3, 3> &src2){
	Eigen::Quaternion<T> qt(src.transpose() * src2);
	return acos(qt.w());
}