//The mass header file

#include "mass.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <armadillo>
#include <array>     // to be able to use this library complie using this g++ -std=c++0x FileNames

using namespace std;

int main()
{
	FILE * hessian;
	hessian=fopen ("h2o_hessian.txt","r");
	int natom;
	fscanf(hessian,"%d",&natom);
	double **H=new double* [natom*3];
	for(int i=0; i<natom*3;i++)
	{
		H[i]=new double[natom*3];
	}
	for(int i=0; i<natom*3; i++)
	{
		for(int j=0;j<natom;j++)
		{
			fscanf(hessian,"%lf %lf %lf",&H[i][3*j], &H[i][3*j+1], &H[i][3*j+2]);
		}
	}
	
	for(int i=0; i<natom*3; i++)
	{
		for(int j=0;j<natom*3;j++)
		{
		printf("%10.7f ", H[i][j]);
		}
		cout<<endl;
	}
}
