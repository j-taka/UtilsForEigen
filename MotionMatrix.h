/*! homogeneous matrix for representing 6-DOF motion */
// actually I would like to use MotionMatrixd ...
#include <Eigen/Dense>

template<typename Type>
class MotionMatrix
{
private:
    /*! rotation */
	Eigen::Matrix<Type, 3, 3, 0, 3, 3> r;
    /*! translation */
	Eigen::Matrix<Type, 3, 1, 0, 3, 1> t;
    
public:
    /*! default constructor */
	MotionMatrix() { 
		r = Eigen::Matrix<Type, 3, 3, 0, 3, 3>::Identity(); 
		t = Eigen::Matrix<Type, 3, 1, 0, 3, 1>::Zero(); 
	}
    /*! constructor 
	 *  \param location location
	 *  \param rpy representation of orientation
	 */
	MotionMatrix(const Eigen::Matrix<Type, 3, 1, 0, 3, 1> &loc, const Eigen::Matrix<Type, 3, 3, 0, 3, 3> &ori){
		r = ori;
		t = loc;
	}
	MotionMatrix(const MotionMatrix<Type> &src){
		r = src.r;
		t = src.t;
	}
    /*! access to rotational matrix 
	 *  \return rotational matrix
	 */
	Eigen::Matrix<Type, 3, 3, 0, 3, 3>& R(){ return r; }
    /*! get rotational matrix 
	 *  \return rotational matrix
	 */
	const Eigen::Matrix<Type, 3, 3, 0, 3, 3>& R() const{ return r; }
    /*! access to location vector 
	 *  \return location vector
	 */
	Eigen::Matrix<Type, 3, 1, 0, 3, 1>& T(){ return t; }
	/*! get location vector
	 *  \return location vector
	 */
	const Eigen::Matrix<Type, 3, 1, 0, 3, 1>& T() const{ return t; }
    /*! substitution 
	 *  \return this
	 *  \param  src  homogeneous matrix
	 */
	MotionMatrix<Type>& operator=(const MotionMatrix<Type> &src){
		r = src.r;
		t = src.t;
		return *this;
	}
    /*! multiple to homogeneous matrix (apply some displacement)
	 *  \return result of the multiplication
	 *  \param  src multiplicand
	 */
	MotionMatrix<Type> operator*(const MotionMatrix<Type> &src) const{
		MotionMatrix<Type> dest;
		dest.r = r * src.r;
		dest.t = r * src.t + t;
		return dest;
	}
	/*! multiple to vector (the conversion of the coordinate system)
	 *  \return result of the multiplication
	 *  \param  src multiplicand
	 */
	Eigen::Matrix<Type, 3, 1, 0, 3, 1> operator*(const Eigen::Matrix<Type, 3, 1, 0, 3, 1> &src) const{
		Eigen::Matrix<Type, 3, 1, 0, 3, 1> dest;
		dest = r * src + t;
		return dest;
	}
    /*! get inverse matrix 
	 *  \return inverse matrix
	 */
	MotionMatrix<Type> Inv() const{
		MotionMatrix<Type> dest;
		dest.r = r.Trans();
		dest.t = -dest.r * t;
		return dest;
	}
	/* inverse matrix */
	MotionMatrix<Type> operator~() const{
		MotionMatrix<Type> dest;
		dest.r = r.Trans();
		dest.t = -dest.r * t;
		return dest;
	}
    /*! convert homogeneous matrix into XYZ-RPY representation 
	 *  \retval xyz xyz value
	 *  \retval rpy rpy value
	 */
	void MM2XYZRPY(Eigen::Matrix<Type, 3, 1, 0, 3, 1> &xyz,Eigen::Matrix<Type, 3, 1, 0, 3, 1> &rpy) const{
	    rpy = Mat2RPY(r);
	    xyz = t;
	}
    /*! convert XYZ-RPY representation to homogeneous matrix 
	 *  \param xyz xyz value
	 *  \param rpy rpy value
	 */
	void XYZRPY2MM(const Eigen::Matrix<Type, 3, 1, 0, 3, 1> &xyz,const Eigen::Matrix<Type, 3, 1, 0, 3, 1> &rpy){
		r = RPY2Mat(rpy);
		t = xyz;
	}
    /*! print */
	friend std::ostream& operator<<(std::ostream &s,const MotionMatrix<Type> &mm){
	    s << mm.r[0][0] << " " << mm.r[0][1] << " " << " " << mm.r[0][2]
	      << "  " << mm.t[0] << std::endl;
	    s << mm.r[1][0] << " " << mm.r[1][1] << " " << " " << mm.r[1][2]
	      << "  " << mm.t[1] << std::endl;
	    s << mm.r[2][0] << " " << mm.r[2][1] << " " << " " << mm.r[2][2]
	      << "  " << mm.t[2] << std::endl;
	    s << "0 0 0 1" << std::endl;
		return s;
	}
};
