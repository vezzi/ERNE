/*
 *  This file is part of ERNE.
 *  Copyright (c) 2011 by Cristian Del Fabbro <delfabbro@appliedgenomics.org>,
 *  Francesco Vezzi <vezzi@appliedgenomics.org>,
 *  Alexandru Tomescu <alexandru.tomescu@uniud.it>
 *  Nicola Prezza <nicolapr@gmail.com>, and
 *  Alberto Policriti <policriti@uniud.it>
 *
 *   ERNE is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   ERNE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with ERNE.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "MismatchCombination.h"
#include "stdio.h"

void MismatchCombination::swap(int i, int j){

	char temp1;
	double temp2;

	temp1 = numberOfErrors[i*2];
	numberOfErrors[i*2] = numberOfErrors[j*2];
	numberOfErrors[j*2] = temp1;

	temp1 = numberOfErrors[i*2+1];
	numberOfErrors[i*2+1] = numberOfErrors[j*2+1];
	numberOfErrors[j*2+1] = temp1;;

	temp2 = errorValues[i];
	errorValues[i] = errorValues[j];
	errorValues[j] = temp2;
}//swap

void MismatchCombination::order(){//order with respect to d

	for(int j=size;j>1;j--)
		for(int i=1;i<j;i++)
			if( errorValues[i-1] > errorValues[i] ) swap(i-1,i);

}//orderDynamicArray

void MismatchCombination::insert(char a, char g, double d){//a=number of alpha-mismatch, g=number of 1-mismatch, d=a*alpha+g

	size++;
	numberOfErrors[2*(size-1)] = a;
	numberOfErrors[2*(size-1)+1] = g;
	errorValues[size-1] = d;

}//insert

char MismatchCombination::numberOfAlphaAt(int i){
	return numberOfErrors[2*i];
}
char MismatchCombination::numberOfOneAt(int i){
	return numberOfErrors[2*i+1];
}
double MismatchCombination::errorValueAt(int i){
	return errorValues[i];
}

void MismatchCombination::initMismatchCombination(double alpha, double max_d, int k){

	size = 0;
	int max_number_of_alpha = 5;

	for(int g=0; g<=(int)max_d; g++)
		for( int a=0; a<=(max_d-g)/alpha+0.0000001 && a <= max_number_of_alpha ;a++ )
				if(a+g<k) size++;

	numberOfErrors = new char[size*2];
	errorValues = new double[size];

	size = 0;

	for(int g=0; g<=(int)max_d; g++)
		for( int a=0; a<=(max_d-g)/alpha+0.000000001 && a <= max_number_of_alpha;a++ )
				if(a+g<k) insert( a,g,a*alpha+g );

	order();
	//for(int i=0;i<size;i++) printf("a=%i g=%i d=%f \n",numberOfAlphaAt(i),numberOfOneAt(i),errorValueAt(i));

}

MismatchCombination::~MismatchCombination(){

	delete [] numberOfErrors;
	delete [] errorValues;

}//destructor

