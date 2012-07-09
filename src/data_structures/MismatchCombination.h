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

#ifndef MISMATCHCOMBINATION_H_
#define MISMATCHCOMBINATION_H_

using namespace std;

#include <stdlib.h>

class MismatchCombination{

	public:

		~MismatchCombination();
		MismatchCombination() {}

		void initMismatchCombination(double alpha, double max_d,int k);
		char numberOfAlphaAt(int i);
		char numberOfOneAt(int i);
		double errorValueAt(int i);

		int size;
		char *numberOfErrors;
		double *errorValues;

	private:

		void swap(int i, int j);
		void order();//orders structs with respect to d
		void insert(char a, char g, double d);

};

#endif /* MISMATCHCOMBINATION_H_ */
