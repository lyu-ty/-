//
// Nonlinear Multicomponent Diffusion
// Explicit finite difference  scheme
// Zero  flux   boundary   conditions
// Skew-symmetric  initial condotions
// Version 0.05        2025 - 03 - 06
// 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define NP 1001
#define NC 1000
#define NM  999
#define NSTP 200000
#define NOUT 2000
#define eps 1E-03
const double
Pi = 3.14159265358979E+00;
double
L0 = 1.00E+00,
C0 = 5.00E+03,
D0 = 1.00E-09,
D1 = 1.50E+00,
D2 = 1.00E+00;
double
Miu = 1.00E+00;
double
XL, XR, XD, XS, HST, TAU, THS, TIM, D11, D22;
double
X[NP], C1old[NP], C2old[NP], C1new[NP], C2new[NP];

//---------------------------------------------------------------------

void Start(void)
{
	int    jc;
	double A, B, C, D, Z;

	TIM = 0.00E+00;
	TAU = 1.00E-06;
	HST = L0 / NC;
	THS = TAU / (HST * HST);
	D11 = D0 * D1 * (C0 / L0) * (C0 / L0) * THS;
	D22 = D0 * D2 * (C0 / L0) * (C0 / L0) * THS;
	XL = 0.49 * L0;
	XR = 0.50 * L0;
	XD = XR - XL;
	XS = L0 - XL;
	A = 1.00E+00;
	B = 0.00E+00;
	C = -3.0 / (XD * XD);
	D = 2.0 / (XD * XD * XD);

	for (jc = 0; jc <= NC; jc++)
	{
		X[jc] = jc * HST;
		//if (X[jc] < XL) { C1old[jc] = 0.8E+00; continue; }
		//if (X[jc] > XR) { C1old[jc] = 0.00E+00; continue; }
		//Z = X[jc] - XL;  
		//C1old[jc] = 0.8*(A + B * Z + C * Z * Z + D * Z * Z * Z);
		double xx = jc*HST;
		C1old[jc] = (0.5)*(1-0.9*cos(3.14159265358979 *xx));
		//printf("ss= %f \n", xx);
	    C2old[jc] = 1 -(C1old[jc])+ 0.001*cos(3.14159265358979 * xx);
	}
	//for (jc = 0; jc <= NC; jc++) 
	return;
}

//---------------------------------------------------------------------

void Step(void)
{
	int jm, jc, jp;
	double C11, C12, C21, C22, DF1, DF2;
	for (jc = 1; jc <= NM; jc++)
	{
		jm = jc - 1;
		jp = jc + 1;
		C11 = D11 * (1.00E+00 - Miu * C2old[jc]);
		C12 = D11 * Miu * C1old[jc];
		C21 = D22 * Miu * C2old[jc];
		C22 = D22 * (1.00E+00 - Miu * C1old[jc]);
		DF1 = C1old[jm] - 2 * C1old[jc] + C1old[jp];
		DF2 = C2old[jm] - 2 * C2old[jc] + C2old[jp];
		C1new[jc] = C1old[jc] + C11 * DF1 + C12 * DF2;
		C2new[jc] = C2old[jc] + C21 * DF1 + C22 * DF2;
	}
	C1new[0] = C1new[1]; C2new[0] = C2new[1];
	C1new[NC] = C1new[NM]; C2new[NC] = C2new[NM];
	for (jc = 0; jc <= NC; jc++)
	{
		C1old[jc] = C1new[jc]; C2old[jc] = C2new[jc];
	}
	TIM += TAU;
	return;
}

//---------------------------------------------------------------------

int Result(int num)
{
	FILE* fp;
	int   jc;
	char TimeLayer[32];
	char FileName[32];
	printf("\n Layer Number =  %d   Time = %e", num, TIM);
	_itoa_s(num, TimeLayer, 10);
	strcpy_s(FileName, "res-");
	strcat_s(FileName, TimeLayer);
	strcat_s(FileName, ".txt");
	//fp = fopen_s(FileName, "w");
	//if (fp == NULL) { printf("\n FAIL"); return 1; }
	if (fopen_s(&fp, FileName, "w") != 0) {
		printf("\n FAIL");
		return 1;
	}
	for (jc = 0; jc < NP; jc++)
		fprintf(fp, " %e %e %e \n", X[jc], C1old[jc], C2old[jc]);
	fclose(fp);
	return 0;
}

//---------------------------------------------------------------------

int main(void)
{
	int m = 0;
	int n = 0;
	Start(); Result(n);
	while (n <= NSTP)
	{
		Step(); if (m == NOUT) { m = 0; Result(n); }  m++; n++;
	}
	printf("\n Ready. Press ENTER to quit."); getchar();
	return 0;
}
