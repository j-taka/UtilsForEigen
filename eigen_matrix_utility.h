/*!
 *  \brief  utilities of eigen matrix
 *  \file   eigen_matrix_utility.h
 *  \author Jun Takamatsu <j-taka@cvl.iis.u-tokyo.ac.jp>
 *  \date   11 Jul. 2004
 */

#pragma once

#include <Eigen/Dense>

const double TRITRUNC = 1.0e-6;

/// a too small number is truncated to zero
template<typename T>
static T check_tol(T num)
{
    if (fabs(num) < (T) TRITRUNC) return 0.0;
    return num;
}

template<typename T>
static inline T ProperTriFunc(T src)
{
	if (src < -1.0){ return (T) -1.0; }
	if (src >  1.0){ return (T) 1.0; }
	return src;
}

template<typename T>
Eigen::Matrix<T, 3, 3, 0, 3, 3> dRotX2Mat(const T &src)
{
	Eigen::Matrix<T, 3, 3, 0, 3, 3> mat;
	T c = (T) cos(src);
	T s = (T) sin(src);
	mat(0, 0) =  0;
	mat(0, 1) =  0;
	mat(0, 2) =  0;
	mat(1, 0) =  0;
	mat(1, 1) = -s;
	mat(1, 2) = -c;
	mat(2, 0) =  0;
	mat(2, 1) =  c;
	mat(2, 2) = -s;
	return mat;
}

template<typename T>
Eigen::Matrix<T, 3, 3, 0, 3, 3> dRotY2Mat(const T &src)
{
	Eigen::Matrix<T, 3, 3, 0, 3, 3> mat;
	T c = (T) cos(src);
	T s = (T) sin(src);
	mat(0, 0) = -s;
	mat(0, 1) =  0;
	mat(0, 2) =  c;
	mat(1, 0) =  0;
	mat(1, 1) =  0;
	mat(1, 2) =  0;
	mat(2, 0) = -c;
	mat(2, 1) =  0;
	mat(2, 2) = -s;
	return mat;
}

template<typename T>
Eigen::Matrix<T, 3, 3, 0, 3, 3> dRotZ2Mat(const T &src)
{
	Eigen::Matrix<T, 3, 3, 0, 3, 3> mat;
	T c = (T) cos(src);
	T s = (T) sin(src);
	mat(0, 0) = -s;
	mat(0, 1) = -c;
	mat(0, 2) =  0;
	mat(1, 0) =  c;
	mat(1, 1) = -s;
	mat(1, 2) =  0;
	mat(2, 0) =  0;
	mat(2, 1) =  0;
	mat(2, 2) =  0;
	return mat;
}

/*! rpy representation to rotation matrix 
 *  \return rotational matrix
 *  \param  rpy rpy representation
 */
template<typename T>
Eigen::Matrix<T, 3, 3, 0, 3, 3> RPY2Mat(const Eigen::Matrix<T, 3, 1, 0, 3, 1> &rpy)
{
	Eigen::Matrix<T, 3, 3, 0, 3, 3> mat;
    T cr = (T) cos(rpy[0]);
    T cp = (T) cos(rpy[1]);
    T cy = (T) cos(rpy[2]);
    T sr = (T) sin(rpy[0]);
    T sp = (T) sin(rpy[1]);
    T sy = (T) sin(rpy[2]);
    mat(0, 0) = check_tol(cr * cp);
    mat(0, 1) = check_tol((cr * sp * sy) - (sr * cy));
    mat(0, 2) = check_tol((cr * sp * cy) + (sr * sy));
    mat(1, 0) = check_tol(sr * cp);
    mat(1, 1) = check_tol((sr * sp * sy) + (cr * cy));
    mat(1, 2) = check_tol((sr * sp * cy) - (cr * sy));
    mat(2, 0) = check_tol(-sp);
    mat(2, 1) = check_tol(cp * sy);
    mat(2, 2) = check_tol(cp * cy);
    return mat;
}

/*! the matrix to be derivatived by r
 *  \return derivatived matrix
 *  \param  rpy rpy representation
 */
template<typename T>
Eigen::Matrix<T, 3, 3, 0, 3, 3> RPYdR2Mat(const Eigen::Matrix<T, 3, 1, 0, 3, 1> &rpy)
{
	Eigen::Matrix<T, 3, 3, 0, 3, 3> mat;
    T cr = (T) cos(rpy[0]);
    T cp = (T) cos(rpy[1]);
    T cy = (T) cos(rpy[2]);
    T sr = (T) sin(rpy[0]);
    T sp = (T) sin(rpy[1]);
    T sy = (T) sin(rpy[2]);
    mat(0, 0) = check_tol(-sr * cp);
    mat(0, 1) = check_tol((-sr * sp * sy) - (cr * cy));
    mat(0, 2) = check_tol((-sr * sp * cy) + (cr * sy));
    mat(1, 0) = check_tol(cr * cp);
    mat(1, 1) = check_tol((cr * sp * sy) - (sr * cy));
    mat(1, 2) = check_tol((cr * sp * cy) + (sr * sy));
    mat(2, 0) = 0;
    mat(2, 1) = 0;
    mat(2, 2) = 0;
    return mat;
}

/*! the matrix to be derivatived by p 
 *  \return derivatived matrix
 *  \param  rpy rpy representation
 */
template<typename T> 
Eigen::Matrix<T, 3, 3, 0, 3, 3> RPYdP2Mat(const Eigen::Matrix<T, 3, 1, 0, 3, 1> &rpy)
{
    Eigen::Matrix<T, 3, 3, 0, 3, 3> mat;
    T cr = (T) cos(rpy[0]);
    T cp = (T) cos(rpy[1]);
    T cy = (T) cos(rpy[2]);
    T sr = (T) sin(rpy[0]);
    T sp = (T) sin(rpy[1]);
    T sy = (T) sin(rpy[2]);
    mat(0, 0) = check_tol(cr * -sp);
    mat(0, 1) = check_tol(cr * cp * sy);
    mat(0, 2) = check_tol(cr * cp * cy);
    mat(1, 0) = check_tol(sr * -sp);
    mat(1, 1) = check_tol(sr * cp * sy); 
    mat(1, 2) = check_tol(sr * cp * cy);
    mat(2, 0) = check_tol(-cp);
    mat(2, 1) = check_tol(-sp * sy);
    mat(2, 2) = check_tol(-sp * cy);
    return mat;
}

/*! the matrix to be derivatived by y 
 *  \return derivatived matrix
 *  \param  rpy rpy representation
 */
template<typename T>
Eigen::Matrix<T, 3, 3, 0, 3, 3> RPYdY2Mat(const Eigen::Matrix<T, 3, 1, 0, 3, 1> &rpy)
{
    Eigen::Matrix<T, 3, 3, 0, 3, 3> mat;
    T cr = (T) cos(rpy[0]);
    T cp = (T) cos(rpy[1]);
    T cy = (T) cos(rpy[2]);
    T sr = (T) sin(rpy[0]);
    T sp = (T) sin(rpy[1]);
    T sy = (T) sin(rpy[2]);
    mat(0, 0) = 0;
    mat(0, 1) = check_tol((cr * sp * cy) + (sr * sy));
    mat(0, 2) = check_tol((cr * sp * -sy) + (sr * cy));
    mat(1, 0) = 0;
    mat(1, 1) = check_tol((sr * sp * cy) - (cr * sy));
    mat(1, 2) = check_tol((sr * sp * -sy) - (cr * cy));
    mat(2, 0) = 0;
    mat(2, 1) = check_tol(cp * cy);
    mat(2, 2) = check_tol(cp * -sy);
    return mat;
}

/*! rotation matrix to rpy representation 
 *  \return rpy representation
 *  \param  rotational matrix
 */
template<typename T>
Eigen::Matrix<T, 3, 1, 0, 3, 1> Mat2RPY(const Eigen::Matrix<T, 3, 3, 0, 3, 3> &mat)
{
	Eigen::Matrix<T, 3, 1, 0, 3, 1> rpy;
    T n1 = mat(0, 0);
    T n2 = mat(1, 0);
    T n3 = mat(2, 0);
    T o1 = mat(0, 1);
    T o2 = mat(1, 1);
    T a1 = mat(0, 2);
    T a2 = mat(1, 2);
    
    rpy[0] = (T) atan2(n2, n1);
    rpy[1] = (T) atan2(-n3, (cos(rpy[0]) * n1) + (sin(rpy[0]) * n2));
    rpy[2] = (T) atan2((sin(rpy[0]) * a1) - (cos(rpy[0]) * a2),
                         (cos(rpy[0]) * o2) - (sin(rpy[0]) * o1));
    return rpy;
}

template<typename T>
void Mat2AngleAxis(T &angle, Eigen::Matrix<T, 3, 1, 0, 3, 1> &axis, const Eigen::Matrix<T, 3, 3, 0, 3, 3> &mat)
{
	// convert to quaternion
	Eigen::Quaternion<T> q(mat);
	// get angle axis
	const T thr = (T) (1 - 1.0e-8);
	if (q.w() >= thr || q.w() <= -thr){ // no rotation
		angle = 0.0;
		axis = Eigen::Matrix<T, 3, 1, 0, 3, 1>(0, 0, 1); // z-axiz (no meaning)
	}
	else{
		angle = 2.0 * acos(q.w());
		axis[0] = q.x();
		axis[1] = q.y();
		axis[2] = q.z();
		axis /= axis.norm();
	}
}

template<typename T>
Eigen::Matrix<T, 3, 3, 0, 3, 3> dAngleAxis2Mat(const Eigen::Matrix<T, 3, 1, 0, 3, 1> &axis, const T &angle)
{
	const double c = cos(angle);
	const double s = sin(angle);
	Eigen::Matrix<T, 3, 3, 0, 3, 3> res;

	res(0, 0) = -s + axis[0] * axis[0] * s;
	res(0, 1) = axis[0] * axis[1] * s - axis[2] * c;
	res(0, 2) = axis[0] * axis[2] * s + axis[1] * c;
	res(1, 0) = axis[0] * axis[1] * s + axis[2] * c;
	res(1, 1) = -s + axis[1] * axis[1] * s;
	res(1, 2) = axis[1] * axis[2] * s -	axis[0] * c;
	res(2, 0) = axis[0] + axis[2] * s - axis[1] * c;
	res(2, 1) = axis[1] * axis[2] * s + axis[0] * c;
	res(2, 2) = -s + axis[2] * axis[2] * s;
	return res;
}

#if 0
// define only one axis
// assign the other axes of orthogonal coordinates.
template<int dim, typename T>
Eigen::Matrix<T, dim, dim, 0, dim, dim> SetCoord(const Eigen::Matrix<T, dim, 1, 0, dim, 1> &axis)
{
	

	CVLVector<dim, double> workspace;
	CVLMatrix<dim, dim, T> coord;
	T swap;
	// 共分散行列の決定
	for (i = 0; i < dim; i++){
		coord[i][i] = axis[i] * axis[i];
		for (j = 0; j < dim; j++){
			coord[j][i] = coord[i][j] = axis[i] * axis[j];
		}
	}
	// trick
	v = CVLMatrixTo2DArray(coord);
	eigen(dim, v, d);
	// メモリの解放
	FreeConverted2DArray(v);
	// ソーティング
	wmax = d[0]; pos = 0;
	for (i = 1; i < dim; i++){
		if (wmax < d[i]){ pos = i; wmax = d[pos]; }
	}
	// 大きいものを適切な軸方向に
	for (i = 0; i < dim; i++){
		swap = coord[pos][i];
		coord[pos][i] = coord[num][i];
		coord[num][i] = swap;
	}
	return coord;
}

template<class T>
Eigen::Matrix<T, 3, 3, 0, 3, 3> SetCoord3(unsigned int num, const Eigen::Matrix<T, 3, 1, 0, 3, 1> &axis, bool rightHandCoord = true)
{
	assert(num < 3);
	unsigned int i, j, pos;
	double **v, d[3], wmax;
	CVLVector<3, double> workspace;
	Eigen::Matrix<T, 3, 3, 0, 3, 3> coord;
	// 共分散行列の決定
	for (i = 0; i < 3; i++){
		coord[i][i] = axis[i] * axis[i];
		for (j = 0; j < 3; j++){
			coord[j][i] = coord[i][j] = axis[i] * axis[j];
		}
	}
	// trick
	v = CVLMatrixTo2DArray(coord);
	eigen(3, v, d);
	// メモリの解放
	FreeConverted2DArray(v);
	// ソーティング
	wmax = d[0]; pos = 0;
	for (i = 1; i < 3; i++){
		if (wmax < d[i]){ pos = i; wmax = d[pos]; }
	}
	// 大きいものを適切な軸方向に
	coord[pos][0] = coord[num][0]; coord[pos][1] = coord[num][1]; coord[pos][2] = coord[num][2];
	workspace = coord.Row(0).Cross(coord.Row(1));
	if (rightHandCoord){
		coord[2][0] = workspace[0]; coord[2][1] = workspace[1]; coord[2][2] = workspace[2];
	}
	else{
		coord[2][0] = -workspace[0]; coord[2][1] = -workspace[1]; coord[2][2] = -workspace[2];
	}
	return coord;
}

/*!
 *	target and org must be normalized.
 */
template<class T>
Eigen::Matrix<T, 3, 3, 0, 3, 3> AxisAlign(const Eigen::Matrix<T, 3, 1, 0, 3, 1> &target, const Eigen::Matrix<T, 3, 1, 0, 3, 1> &org)
{
	Eigen::Matrix<T, 3, 1, 0, 3, 1> rotAxis;
	Eigen::Matrix<T, 3, 3, 0, 3, 3> mat;
	T angle; 
	rotAxis = org.Cross(target);
	angle = atan2(rotAxis.Length(), target.Dot(org));
	if (rotAxis.Length() == 0.0 && target.Dot(org) < 0){		
		CVLVector<3, double> temp;
		temp = org;
		rotAxis = SetCoord(0, temp).Row(0);
	}
	mat = AA2M(rotAxis, angle);
	return mat;
}
#endif
