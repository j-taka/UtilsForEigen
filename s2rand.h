#ifndef _S2RAND_H
#define _S2RAND_H

#include <cmath>
#include <cstdlib>

static const double X = .525731112119133606;
static const double Z = .850650808352039932;

/* vertex of Icosahedron */
static const double vp[12][3]=
{
    {-X, 0, Z}, { X, 0, Z}, {-X, 0,-Z}, { X, 0,-Z},
    { 0, Z, X}, { 0, Z,-X}, { 0,-Z, X}, { 0,-Z,-X},
    { Z, X, 0}, {-Z, X, 0}, { Z,-X, 0}, {-Z,-X, 0}
};      

/* connectivity of Icosahedron */
static const int fn[20][3]=
{
    { 0, 4, 1}, { 0, 9, 4}, { 9, 5, 4}, { 4, 5, 8}, { 4, 8, 1},
    { 8,10, 1}, { 8, 3,10}, { 5, 3, 8}, { 5, 2, 3}, { 2, 7, 3},
    { 7,10, 3}, { 7, 6,10}, { 7,11, 6}, {11, 0, 6}, { 0, 1, 6},
    { 6, 1,10}, { 9, 0,11}, { 9,11, 2}, { 9, 2, 5}, { 7, 2,11}
};

// convert degree to radian
static double DegToRad(double deg)
{
    return deg * M_PI / 180.0;
}

void s2rand(double dest[3])
{
    double prop0, prop1, prop2;
    int triNum = (int) ((rand() / ((double) RAND_MAX + 1.0)) * 20);
    for (;;){
        prop1 = rand() / (double) RAND_MAX;
        prop2 = rand() / (double) RAND_MAX;
        if (prop1 + prop2 <= 1) break;
    }
    prop0 = 1 - prop1 - prop2;
	for (int i(0); i < 3; ++i){
        dest[i] = prop0 * vp[fn[triNum][0]][i] + prop1 * vp[fn[triNum][1]][i]
				+ prop2 * vp[fn[triNum][2]][i];
	}
    double length = 0.0;
	for (int i(0); i < 3; ++i){
		length += dest[i] * dest[i];
	}
	length = sqrt(length);
	for (int i(0); i < 3; i++){
		dest[i] /= length;
	}
    return;
}

template<class T> void s2rand(T &dest)
{
    typename T::VectorType prop0, prop1, prop2, length;
    int triNum = rand() % 20, i;
    for (;;){
		prop1 = rand() / (typename T::VectorType) RAND_MAX;
		prop2 = rand() / (typename T::VectorType) RAND_MAX;
        if (prop1 + prop2 <= 1) break;
    }
    prop0 = 1 - prop1 - prop2;
	for (i = 0; i < 3; i++){
		dest[i] = prop0 * (typename T::VectorType)vp[fn[triNum][0]][i] 
				+ prop1 * (typename T::VectorType)vp[fn[triNum][1]][i]
				+ prop2 * (typename T::VectorType)vp[fn[triNum][2]][i];
	}
    length = 0.0;
	for (i = 0; i < 3; i++){
		length += dest[i] * dest[i];
	}
	length = sqrt(length);
	for (i = 0; i < 3; i++){
		dest[i] /= length;
	}
    return;
}

#endif /* _S2RAND_H */
