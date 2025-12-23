#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define LHX 2.0E+00
#define LHZ 3.0E+00
#define RAD 4.0E-02  
#define NSX 20
#define NSZ 30

//  -------------------------------------------------------------------
//  Mesh Generator 001
//  Rectangular Regular Mesh
//  Date:    2025-08-30
//  -------------------------------------------------------------------
//  LHX - Length of the region in X direction
//  LHZ - Length of the region in Z direction
//  RAD - Radius of the spherical coal particle
//  NSX  - Number of spherical particles in X direction
//  NSZ  - Number of spherical particles in Z direction
//  -------------------------------------------------------------------

int main(void)
{
	int j, k;
	double CTX, CTZ;
	double DX1, DX2, DZ1, DZ2;
	FILE* fp;
	DX1 = LHX / NSX;
	DX2 = 0.5 * DX1;
	DZ1 = LHZ / NSZ;
	DZ2 = 0.5 * DZ1;
	//fp = fopen("Coal-Part-Type-I.txt", "w");
	//if (fp == NULL) { printf("\n Cannot make Coal-Part-Type-I.txt."); return 1; }
	errno_t err = fopen_s(&fp, "Coal-Part-Type-I.txt", "w");
	if (err != 0 || fp == NULL) {
		printf("\n Cannot make Coal-Part-Type-I.txt.");
		return 1;
	}
	for (j = 0; j < NSX; ++j)
	{
		CTX = DX2 + j * DX1;
		for (k = 0; k < NSZ; ++k)
		{
			CTZ = DZ2 + k * DZ1;
			fprintf(fp, " %6.3f %6.3f %6.3f \n", CTZ, CTX, RAD);
		}
	}
	fclose(fp);
	printf("\n Rectangular regular mesh generated.");
	printf("\n LHX=%6.3f  LHZ=%6.3f  RAD=%6.3f", LHX, LHZ, RAD);
	printf("\n NSX=%4d    NSZ=%4d    NST=%4d", NSX, NSZ, NSX * NSZ);

	printf("\n Ready.");
	return 0;
}
