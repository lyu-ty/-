#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//  -------------------------------------------------------------------
//  LHX - Length of the region in X direction
//  LHZ - Length of the region in Z direction
//  NST - Total number of spherical particles
//  PX  - Number of mesh points in X direction
//  PZ  - Number of mesh points in Z direction
//  -------------------------------------------------------------------

#define LHX 2.0E+00
#define LHZ 3.0E+00
#define NST 750
#define PX  201
#define PZ  301

//  -------------------------------------------------------------------
//  Model of the Shale Gas Transport 
//  Version 1.00
//  Date:   2025-09-10

//  Discontinuous diffusion coefficient
//  Conservative finite difference scheme
//  Dirichlet boundary condition at the left edge
//  Robin boundary condition at the right edge
//
//  2D Diffusion equation in (Z,X) with D = D(Z,X)
//  Square rectangular region [0:LZ,0:LX] 
//  Spherical carbon particles without inner boundaries
//
//  T = 0  - Initial condition       : U = 0 
//  Z = 0  - 1-st boundary condition : U = U0 
//  Z = LZ - 3-rd boundary condition : D*DU/DZ + U = 0
//  X = 0  - 2-nd boundary condition : DU/DX = 0
//  X = LX - 2-nd boundary condition : DU/DX = 0
//
//  Explicit finite difference scheme
//
//  -------------------------------------------------------------------
double DT = 1.0E-03;
double D0 = 1.0E-06;
double D1 = 1.0E-03;
double U0 = 1.0E+00;

int NStep = 200000;
int NSave = 10000;

double CZ[NST], CX[NST], CR[NST];

int    NX, MX, NZ, MZ, NS;
double HX, HZ, CFI, CFE, DLZ, Time;
double X[PX], Z[PZ], D[PZ][PX];
double DTXP[PZ][PX], DTXM[PZ][PX];
double DTZP[PZ][PX], DTZM[PZ][PX];
double Uold[PZ][PX], Unew[PZ][PX];

int Kaishi(void)
{
	int j, k;
	FILE* fp;

	//fp = fopen("Coal-Part-Type-I.txt", "r");
	//if (fp == NULL) { printf("\n Coal-Part-Type-I.txt not found."); return 1; }

	errno_t err;
	err = fopen_s(&fp, "Coal-Part-Type-I.txt", "r");
	if (err != 0 || fp == NULL) {
		printf("\n Coal-Part-Type-I.txt not found.");
		return 1;
	}
	NS = 0;
	while (NS < NST && fscanf_s(fp, "%lf %lf %lf", &CZ[NS], &CX[NS], &CR[NS]) == 3) {
		printf("\n NS=%d CZ=%6.3f CX=%6.3f CR=%6.3f", NS, CZ[NS], CX[NS], CR[NS]);
		NS++;
	}
	//while (fscanf(fp, "%lf %lf %lf", &CZ[NS], &CX[NS], &CR[NS]) != EOF)
	//{
	//	printf("\n NS=%d CZ=%6.3f CX=%6.3f CR=%6.3f", NS, CZ[NS], CX[NS], CR[NS]); NS++;
	//}
	fclose(fp);

	NX = PX - 1; MX = NX - 1;
	NZ = PZ - 1; MZ = NZ - 1;
	NZ = PZ - 1; MZ = NZ - 1;

	HX = LHX / NX;
	HZ = LHZ / NZ;
	CFI = HX / 3;
	CFE = D1 / (D1 + HZ);
	DLZ = 1.0E+00 / (D1 + LHZ);

	for (j = 0; j < PX; ++j) X[j] = j * HX;
	for (k = 0; k < PZ; ++k) Z[k] = k * HZ;

	for (j = 0; j < PX; ++j)
		for (k = 1; k < NZ; ++k)
			Uold[k][j] = 0.0E+00;

	for (j = 0; j < PX; ++j)
	{
		Uold[0][j] = U0; Uold[NZ][j] = 0.0E+00;
	}
	return 0;
}

void WanGe(void)
{
	int    j, k, m;
	double dx, dz, dr;
	double DTZZ, DTXX;
	for (j = 0; j < PX; ++j)
		for (k = 0; k < PZ; ++k)
		{
			D[k][j] = D1;
			for (m = 0; m < NS; ++m)
			{
				dx = X[j] - CX[m];
				dz = Z[k] - CZ[m];
				dr = sqrt(dx * dx + dz * dz);
				if (dr > CR[m]) continue;
				D[k][j] = D0;
			}
		}
	DTZZ = 0.5 * DT / (HZ * HZ);
	DTXX = 0.5 * DT / (HX * HX);
	for (j = 0; j < PX; ++j)
		for (k = 0; k < PX; ++k)
		{
			DTZP[k][j] = (D[k][j] + D[k + 1][j]) * DTZZ;
			DTZM[k][j] = (D[k][j] + D[k - 1][j]) * DTZZ;
			DTXP[k][j] = (D[k][j] + D[k][j + 1]) * DTXX;
			DTXM[k][j] = (D[k][j] + D[k][j - 1]) * DTXX;
		}
	return;
}

void Bu(void)
{
	int j, k;
	double SMZP, SMZM, SMXP, SMXM;
	for (j = 1; j < NX; ++j)
		for (k = 1; k < NZ; ++k)
		{
			SMZP = DTZP[k][j] * (Uold[k + 1][j] - Uold[k][j]);
			SMZM = DTZM[k][j] * (Uold[k][j] - Uold[k - 1][j]);
			SMXP = DTXP[k][j] * (Uold[k][j + 1] - Uold[k][j]);
			SMXM = DTXM[k][j] * (Uold[k][j] - Uold[k][j - 1]);
			Unew[k][j] = Uold[k][j] + (SMZP - SMZM) + (SMXP - SMXM);
		}
	for (j = 0; j < PX; ++j) Unew[0][j] = U0;
	for (j = 0; j < PX; ++j) Unew[NZ][j] = CFE * Unew[MZ][j];
	for (k = 0; k < PZ; ++k) { Unew[k][0] = Unew[k][1]; Unew[k][NX] = Unew[k][MX]; }

	for (j = 0; j < PX; ++j)
		for (k = 0; k < PZ; ++k)
			Uold[k][j] = Unew[k][j];
	Time += DT;
	return;
}

int Baocun(int num)
{
	FILE* fp;
	FILE* fq;
	int  j, k;
	char TimeLayer[32];
	char FileName1[32];
	char FileName2[32];
	errno_t err;


	double Uexact, LocDif;

	/*
	for (k = 0; k < PZ; k++) {
		double z = Z[k];
		double u_exact = 1.0 - z / (D_coeff + L);


		for (j = 0; j < PX; j++) {
			double diff = fabs(Uold[k][j] - u_exact);
			if (diff > max_diff) {
				max_diff = diff;
			}
		}
	}
	*/
	printf("\n Time = %6.3e", Time);

	//itoa(num, TimeLayer, 10);
	err = _itoa_s(num, TimeLayer, sizeof(TimeLayer), 10);
	if (err != 0) {
		printf("\n Integer conversion failed.");
		return 1;
	}
	//strcpy(FileName, "res-"); strcat(FileName, TimeLayer); strcat(FileName, ".txt");
	err = strcpy_s(FileName1, sizeof(FileName1), "res-");
	if (err != 0) {
		printf("\n String copy failed.");
		return 1;
	}

	err = strcat_s(FileName1, sizeof(FileName1), TimeLayer);
	if (err != 0) {
		printf("\n String concatenation failed.");
		return 1;
	}

	err = strcat_s(FileName1, sizeof(FileName1), ".txt");
	if (err != 0) {
		printf("\n String concatenation failed.");
		return 1;
	}
	//fp = fopen(FileName, "w"); if (fp == NULL) { printf("\n FAIL"); return 1; }
	err = fopen_s(&fp, FileName1, "w");
	if (err != 0 || fp == NULL) {
		printf("\n FAIL to open file for writing.");
		return 1;
	}

	//strcpy(FileName, "res-"); strcat(FileName, TimeLayer); strcat(FileName, ".txt");
	err = strcpy_s(FileName2, sizeof(FileName2), "dif-");
	if (err != 0) {
		printf("\n String copy failed.");
		return 1;
	}

	err = strcat_s(FileName2, sizeof(FileName2), TimeLayer);
	if (err != 0) {
		printf("\n String concatenation failed.");
		return 1;
	}

	err = strcat_s(FileName2, sizeof(FileName2), ".txt");
	if (err != 0) {
		printf("\n String concatenation failed.");
		return 1;
	}
	//fp = fopen(FileName, "w"); if (fp == NULL) { printf("\n FAIL"); return 1; }
	err = fopen_s(&fq, FileName2, "w");
	if (err != 0 || fq == NULL) {
		printf("\n FAIL to open file for writing.");
		return 1;
	}


	for (k = 0; k < PZ; k++)
	{

		Uexact = 1.0E+00 - Z[k] * DLZ;

		for (j = 0; j < PX; j++)
		{
			LocDif = fabs(Uexact - Uold[k][j]);
			fprintf(fp, " %6.3f %6.3f %6.3f \n", Z[k], X[j], Uold[k][j]);
			fprintf(fq, " %6.3f %6.3f %6.3e \n", Z[k], X[j], LocDif);

		}
		if (k != NZ) fprintf(fp, "\n"); fprintf(fq, "\n");
	}
	fclose(fp);
	fclose(fq);
	return 0;
}

void main(void)
{
	int nst, nsv;
	if (Kaishi()) return; else WanGe();
	Time = 0.0E+00;
	nsv = 0; nst = 0;
	do
	{
		Bu(); nsv++; nst++;
		if (nsv == NSave) { Baocun(nst); nsv = 0; }
	} while (nst < NStep);
	return;
}
