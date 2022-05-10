

#include <iostream>
#include <cmath>
using namespace std;
struct Kepl {
	double M, i, w, W, e, a;
};
struct Dec {
	double x, y, z, Vx, Vy,Vz;
};
struct Vector{
	double x, y, z;
};
Vector VectorMultiply(Vector & R, Vector & V) {
	Vector H;
	H.x = R.y * V.z - R.z * V.y;
	H.y = R.z * V.x - R.x * V.z;
	H.z = R.x * V.y - V.x * R.y;
	return H;
}
Vector operator - (Vector A, Vector B) {
	Vector C;
	C.x = A.x - B.x;
	C.y = A.y - B.y;
	C.z = A.z - B.z;
	return C;
}
Vector operator + (Vector A, Vector B) {
	Vector C;
	C.x = A.x + B.x;
	C.y = A.y + B.y;
	C.z = A.z + B.z;
	return C;
}
double operator * (Vector A, Vector B) {
	double C;
	C = A.x * B.x + A.y * B.y + A.z * B.z;
	return C;
}
Vector operator * (double a, Vector& B) {
	Vector C;
	C.x = a*B.x;
	C.y = a*B.y;
	C.z = a*B.z;
	return C;
}

Vector operator * ( Vector B,double a) {
	Vector C;
	C.x = a * B.x;
	C.y = a * B.y;
	C.z = a * B.z;
	return C;
}

Vector operator / (Vector B, double a) {
	Vector C;
	C.x = B.x/a;
	C.y = B.y/a;
	C.z = B.z/a;
	return C;
}

double NormVector(Vector V) {
	double Norm;
	Norm = sqrt(V.x * V.x + V.y * V.y + V.z * V.z);
	return Norm;
}
Kepl DecToKepl(Vector R, Vector V,double mu) {
	Vector h,e,n;
	double nu,E;
	Kepl K = {0,0,0};
	h = VectorMultiply(R, V);
	e = VectorMultiply(V, h)/mu - R/NormVector(R);
	Vector N = { 0,0,1 };
	n = VectorMultiply(N, h);
	if (R * V >= 0.0) {
		nu = acos(e * R / (NormVector(e) * NormVector(R)));
	}
	else {
		nu = 2 * 3.14 - acos(e * R / (NormVector(e) * NormVector(R)));;
	}
	K.i = acos(h.z / NormVector(h));
	K.e = NormVector(e);
	
	E = 2 * atan(tan(nu / 2) / sqrt((1 + K.e) / (1 - K.e)));
	if (n.y >= 0) {
		K.W = acos(n.x / NormVector(n));
	}
	else {
		K.W = 2 * 3.14 - acos(n.x / NormVector(n));
	}
	
	if (e.z >= 0) {
		if (n * e / (NormVector(e) * NormVector(n)) > 1) {
			K.w = acos(1);
		}
		else {
			K.w = acos(n * e / (NormVector(e) * NormVector(n)));
		}
	}
	else {
		if (n * e / (NormVector(e) * NormVector(n)) > 1) {
			K.w = 2 * 3.14159 - acos(1);
		}
		else {
			K.w = 2 * 3.14159 - acos(n * e / (NormVector(e) * NormVector(n)));
		}
	}
	K.w = fmod(K.w,6.28317);
	K.M = E - K.e * sin(E);
	if (nu > 3.14159265) {
		K.M = 2 * 3.14159265 - K.M;
	}
	K.a = 1 / (2/NormVector(R)-NormVector(V)*NormVector(V)/mu);
	
	return K;
}
Dec KeplToDec(Kepl & K, double mu) {
	Dec D;
	double E, E_, eps, nu, r,u,cosnu,Vr,Vn,p;
	E_ = 0;
	E = K.M;
	eps = 0.000001;
	while (E - E_ > eps) {
			E_ = E;
			E = K.M + K.e * sin(E_);
	}
	cosnu = (1.0 - K.e * K.e) / (K.e * (1.0 - K.e * cos(E))) - 1.0 / K.e;
	nu = acos(cosnu);
	r = K.a * (1 - K.e * cos(E));
	if (K.M > 3.14) {
		u = K.w - nu+2*3.14;
	}
	else {
		u = K.w + nu;
	}
	D.x = r * (cos(K.W) * cos(u) - sin(K.W) * sin(u) * cos(K.i));
	D.y = r * (sin(K.W) * cos(u) + cos(K.W) * sin(u) * cos(K.i));
	D.z = r * (sin(u) * sin(K.i));
	p = K.a * (1 - K.e * K.e);
	Vr = sqrt(mu / p) * K.e * sin(nu);
	Vn = sqrt(mu / p) * (1 + K.e * cos(nu));
	D.Vx = cos(K.W)*(Vr*cos(u)-Vn*sin(u))- cos(K.i)* sin(K.W)*(Vr*sin(u)+Vn*cos(u));
	D.Vy = sin(K.W)*(Vr * cos(u) - Vn * sin(u)) + cos(K.W)*cos(K.i)*(Vr * sin(u) + Vn * cos(u)) ;
	D.Vz =sin(K.i)* (Vr * sin(u) + Vn * cos(u));
	return D;
}



