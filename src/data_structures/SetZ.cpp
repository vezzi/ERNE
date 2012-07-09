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

#include "SetZ.h"
#include "stdio.h"


//---------------------------------TEST--------------


void printEl(unsigned long z){

	int i=0;

	for(i=16-1;i>=0;i--)
		VERBOSE_CHANNEL << (unsigned int)((z >> (i*2))&3);

	VERBOSE_CHANNEL << endl;

}

//---------------------------------END TEST------------------------------------------------------------



void SetZ::initZ(int max_d, int k, bool methyl_hash){

	ZSize=0;

	int bit;

	if(methyl_hash)
		bit = 1;
	else
		bit = 2;

	this->methyl_hash = methyl_hash;
	this->max_d = max_d;

	Zarray = new unsigned long int*[max_d+1];

	for(int partition = 0; partition <= max_d; partition++){

		if(methyl_hash){
			if(partition==0)
				partitionSize=1;
			else if(partition==1)
				partitionSize = k;
			else if(partition==2)
				partitionSize = (k*(k-1))/2;
			else if(partition==3)
				partitionSize = (k*(k-1)*(k-2))/6;
		}else{
			if(partition==0)
				partitionSize=1;
			else if(partition==1)
				partitionSize = k*3;
			else if(partition==2)
				partitionSize = ((k*(k-1))/2)*9;
			else if(partition==3)
				partitionSize = ((k*(k-1)*(k-2))/6)*27;
		}

		Zarray[partition] = new unsigned long int[partitionSize+1];
		Zarray[partition][0] = 0;//first element of each partition is the size of the partition

		buildRecursive((unsigned long int)0,partition, k, partition);

	}//for

}//constructor


SetZ::~SetZ(){

	for(int partition = 0; partition <= this->max_d; partition++) delete [] Zarray[partition];
	delete [] Zarray;

}//destructor

void SetZ::buildRecursive(unsigned long int z, int g, int toInsert,int partition){

	int bit;
	if(methyl_hash)
		bit = 1;
	else
		bit = 2;

	if(toInsert==0 && g == 0) insert(z,partition);
	else if(g == 0 && toInsert > 0){

		z = z<<(toInsert*bit);
		insert(z,partition);

	}else if(g <= toInsert && toInsert>0){

		buildRecursive(z<<bit,g,toInsert-1,partition);//'0' cipher inserted
		buildRecursive((z<<bit)^1,g-1,toInsert-1,partition);//'1' cipher inserted

		if(!methyl_hash){
			buildRecursive((z<<2)^2,g-1,toInsert-1,partition);//'2' cipher inserted
			buildRecursive((z<<2)^3,g-1,toInsert-1,partition);//'3' cipher inserted
		}

	}//else

}//buildrecursive

int SetZ::size(unsigned int partition){

	return Zarray[partition][0];

}//size

unsigned long int SetZ::getElement(unsigned int partition,unsigned int element){

	return Zarray[partition][element+1];

}//getElement

void SetZ::updateIndex(int *partition,int *idxZ,int error,bool *get_z_element){

	if(*idxZ >= size(*partition)-1){//it was the last element in the partition

		if(*partition >= this->max_d ) *get_z_element = false; //it was the last partition
		else{//not the last partition

			if(*partition+1>error) *get_z_element = false;
			else{
				*partition = *partition + 1;
				*idxZ = 0;
			}

		}
	}else *idxZ = *idxZ + 1;	//not the last element in the partition

}//update index

void SetZ::insert(unsigned long int z,int partition){

	//printEl(z);

	ZSize++;
	Zarray[partition][0]++;//increment size
	Zarray[partition][ Zarray[partition][0] ] = z;

}//insert
