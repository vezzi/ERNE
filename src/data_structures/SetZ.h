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

#ifndef SETZ_H_
#define SETZ_H_

#include <common.h>
#include <iostream>

using namespace std;

class SetZ{
public:

	SetZ() {}
	~SetZ();
	unsigned long int **Zarray;//partitions (errors)-> a set of Z elements for every partition
	void initZ(int max_d,int k,bool methyl_hash);
	unsigned long int getElement(unsigned int, unsigned int);//get next element in Z[indexZ]
	int size(unsigned int);
	unsigned long int ZSize;
	bool methyl_hash;
	void updateIndex(int *partition,int *idxZ,int error,bool *get_z_element);

private:

	void insert(unsigned long int z, int partition);
	int max_d;//max distance between blocks and text
	unsigned int partitionSize;
	void buildRecursive(unsigned long int z, int g, int toInsert,int partition);

};

#endif /* SETZ_H_ */
