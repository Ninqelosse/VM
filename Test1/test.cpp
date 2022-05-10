#include "pch.h"
#include "../Vychmeh/Vychmeh.h"
TEST(KeplToDec, ok) {
	double mu = 100000;
	Kepl K, k;
	for (int i = 1; i < 50; i++) {
		K = { i * 0.0628,0.5,0,1.7,0.3,5 };
		Dec d = KeplToDec(K, mu);
		//cout << d.x << " " << d.y << " " << d.z << " "<< d.Vx << " " << d.Vy << " " << d.Vz<< endl;
		Vector R = { d.x, d.y, d.z };
		Vector V = { d.Vx, d.Vy, d.Vz };
		k = DecToKepl(R, V, mu);
		//cout << k.M << " " << k.i << " " << k.w << " " << k.W << " " << k.e << " " << k.a << endl;
		EXPECT_NEAR(K.M, k.M, 0.1);
		EXPECT_NEAR(K.i, k.i, 0.0001);
		EXPECT_NEAR(K.w, k.w, 0.0001);
		EXPECT_NEAR(K.W, k.W, 0.0001);
		EXPECT_NEAR(K.e, k.e, 0.0001);
		EXPECT_NEAR(K.a, k.a, 0.0001);


	}
}
TEST(KeplToDec1, ok1) {
		double mu = 100000;
		Kepl K, k;
		for (int i = 1; i < 100; i++) {
			K = { i * 0.0628,0.5,0,1.7,0.5,6 };
			Dec d = KeplToDec(K, mu);
			//cout << d.x << " " << d.y << " " << d.z << " "<< d.Vx << " " << d.Vy << " " << d.Vz<< endl;
			Vector R = { d.x, d.y, d.z };
			Vector V = { d.Vx, d.Vy, d.Vz };
			k = DecToKepl(R, V, mu);
			//cout << k.M << " " << k.i << " " << k.w << " " << k.W << " " << k.e << " " << k.a << endl;
			//EXPECT_NEAR(K.M, k.M, 0.0001);
			EXPECT_NEAR(K.i, k.i, 0.0001);
			//EXPECT_NEAR(K.w, k.w, 0.0001);
			EXPECT_NEAR(K.W, k.W, 0.0001);
			EXPECT_NEAR(K.e, k.e, 0.0001);
			EXPECT_NEAR(K.a, k.a, 0.0001);


		}

}
TEST(KeplToDec2, ok2) {
	double mu = 100000;
	Kepl K, k;
	for (int i = 1; i < 100; i++) {
		K = { i * 0.0628,0.8,1,1.7,0.5,6 };
		Dec d = KeplToDec(K, mu);
		//cout << d.x << " " << d.y << " " << d.z << " "<< d.Vx << " " << d.Vy << " " << d.Vz<< endl;
		Vector R = { d.x, d.y, d.z };
		Vector V = { d.Vx, d.Vy, d.Vz };
		k = DecToKepl(R, V, mu);
		//cout << k.M << " " << k.i << " " << k.w << " " << k.W << " " << k.e << " " << k.a << endl;
		//EXPECT_NEAR(K.M, k.M, 0.0001);
		EXPECT_NEAR(K.i, k.i, 0.0001);
		//EXPECT_NEAR(K.w, k.w, 0.0001);
		EXPECT_NEAR(K.W, k.W, 0.0001);
		EXPECT_NEAR(K.e, k.e, 0.0001);
		EXPECT_NEAR(K.a, k.a, 0.0001);


	}

}



