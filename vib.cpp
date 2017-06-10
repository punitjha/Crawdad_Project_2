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
#include <complex>

#include <array>   // to be able to use this library complie using this g++ -std=c++0x FileNames

using namespace std;

int main()
{
	//double **geom;
	//std::vector<std::vector<double>> geom;
	arma::mat geom;
	arma::mat hessian;
	arma::vec Zval;
	int  Natom=0;
	int  Natom_h=0;


/**************************************************************/	
// Opening the coordinate file
/**************************************************************/	


	ifstream coord ("h2o_geom.txt", ios::in);
	if (coord.is_open())
	{
		cout<<"The coordinate file is now open\n";
		//first read the number of atoms
		coord>>Natom;
		geom=arma::zeros(Natom,3);
		Zval=arma::zeros(Natom);
		for(unsigned int i=0; i<Natom; i++)
		{
			coord >> Zval(i)>> geom(i,0)>> geom(i,1)>>geom(i,2);
		}
	}
	else cout<<"Sorry! The file could not be opened\n";
	geom.print("\nThis is the coordinate file\n");

/**************************************************************/	
// Opening the Hessian file
/**************************************************************/


	ifstream hess ("h2o_hessian.txt", ios::in);
	if (hess.is_open())
	{
		cout<<"\nThe Hessian file is \n";
		//first read the number of atoms
		hess>>Natom_h;
		if(Natom_h == Natom)
			cout<<"\nTotal no. of atoms in both files consistent\n";
		hessian=arma::zeros(Natom*3,Natom*3);
		for(unsigned int i=0; i<3*Natom; i++)
		{
			for(int j=0; j<Natom*3;j++)
			{
			hess >> hessian(i,j);
			}
		}
	}
	else cout<<"Sorry! The file could not be opened\n";
	hessian.print("\nThis is the Hessian file\n");


/************************************************************************************/	
//divide each element of the Hessian by the product of SR of the masses of the atoms 
/************************************************************************************/

	cout<<"\nThe elements generated for the mass-weighted Hessian\n\n";
	int first=0;
	
	for (int i=0; i<Natom*3; i++)
	{
		if (i !=0 && i % 3  == 0)
			first++;
		int second=0;
		for (int j=0; j<Natom*3;j++)
		{
			cout<<first<<second<<" ";
			hessian(i,j)/=(sqrt(mass[int(Zval[first])]*mass[int(Zval[second])]));
		//	hessian(i,j)*=100000*15.82;
			if (j != 0 && (j+1) % 3 == 0)
				second++;
		}
		cout<<endl;
	}
	hessian.print("\nThis is the mass weighted Hessian\n");

/************************************************************************************/
//diagonalizing the matrix elements to get the eigenvalues
/************************************************************************************/

	arma::cx_vec eigval;
	arma::cx_mat eigvec;
	arma::eig_gen(eigval,eigvec,hessian);
	arma::vec eigval1=arma::sort(real(eigval));
	eigval1.print("\nThe eigenvalues of the Hessian are: (in Hartree/amu-bohr^2)");
	
/************************************************************************************/
//Computing the Harmonic Vibrational Frequencies
/************************************************************************************/
	
	
	cout<<"\n The Harmonic vibrational frequencies are (in cm-1)";
	eigval1=sqrt(eigval1)/(2.0*3.14);
	cout<<endl;
	eigval1.print();




	cout<<"\nNote that there are only three zero frequencies in this case when there \n"
		"should be six. This is because the structure used in the computation is \n"
		"not a stationary point on the potential energy surface, and thus the three \n"
		"rotational frequencis show up.\n"<<endl;



}
