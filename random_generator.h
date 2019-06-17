/*! @file random_generator.h
 *  @brief random generator for s1, s2 space, uniform in sphere, uniform in SO(3) 
 *  @author Jun Takamatsu <j-taka@is.naist.jp>
 */
#include <Eigen/Dense>

/*! @brief random generator (double type only) */
class RandomGenerator
{
public:
	/*! @brief constructor */
	RandomGenerator();
	/*! @brief random number on s1-space 
	 *  @retval random number 
	 */
	Eigen::Vector2d RandInS1() const;
	/*! @brief random number in the unit circle (include the inside)
	 *  @retval random number
	 */
	Eigen::Vector2d RandInUnitCircle() const;
	/*! @brief random number on s2-space 
	 *  @retval random number
	 */
	Eigen::Vector3d RandInS2() const;
	/*! @brief random number in the unit sphere (include the inside) 
	 *  @retval random number
	 */
	Eigen::Vector3d RandInUnitSphere() const;
	/*! @brief random number in rotation matrix
	 *  @details use the implementation in gamegems
	 *  @retval rotation matrix
	 */
	Eigen::Matrix3d RandInRotationMatrix() const;
	double UniformRandom() const{
		return (double) rand() / (double) RAND_MAX;
	}
};