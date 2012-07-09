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

#ifndef GLOBALTOLOCAL_H_
#define GLOBALTOLOCAL_H_

#include "common.h"

#include <string>
#include <iostream>
#include <cstdlib>
#include <set>
using namespace std;

class GlobalToLocal {
public:
	~GlobalToLocal() {
		if (contigs > 0) {
			delete [] startPositions;
			delete [] endPositions;
			delete [] contigsNames;
			contigs = 0;
		}
	}
	GlobalToLocal(const GlobalToLocal &GTL); // copy constructor
	GlobalToLocal & operator=(const GlobalToLocal & GTL);
	GlobalToLocal();
	GlobalToLocal(int numContig);
	int searchContig(unsigned int GlobalCoordinate);
	bool detectDuplicates();

	// TODO rename Contigs -> contigs
	int contigs; //Number of contigs
	unsigned int *startPositions;
	unsigned int *endPositions;
	string *contigsNames;
};

#endif /* GLOBALTOLOCAL_H_ */
