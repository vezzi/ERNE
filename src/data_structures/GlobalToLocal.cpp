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

#include "GlobalToLocal.h"

GlobalToLocal::GlobalToLocal() {
	contigs = 0;
}

GlobalToLocal::GlobalToLocal(int numContig) {
	contigs = numContig;
	startPositions = new unsigned int[numContig];
	endPositions = new unsigned int[numContig];
	contigsNames = new string[numContig];
}

GlobalToLocal::GlobalToLocal(const GlobalToLocal &GTL) {
	contigs = GTL.contigs;
	if (contigs > 0) {
		startPositions = new unsigned int[contigs];
		endPositions = new unsigned int[contigs];
		contigsNames = new string[contigs];
		for (int i = 0; i < contigs; i++) {
			startPositions = GTL.startPositions;
			endPositions = GTL.endPositions;
			contigsNames = GTL.contigsNames;
		}
	}
}

GlobalToLocal & GlobalToLocal::operator=(const GlobalToLocal & GTL) {
	if (this != &GTL) {
		if (contigs > 0) {
			delete [] startPositions;
			delete [] endPositions;
			delete [] contigsNames;
		}
		contigs = GTL.contigs;
		if (contigs > 0) {
			startPositions = new unsigned int[contigs];
			endPositions = new unsigned int[contigs];
			contigsNames = new string[contigs];
			for (int i = 0; i < contigs; i++) {
				startPositions[i] = GTL.startPositions[i];
				endPositions[i] = GTL.endPositions[i];
				contigsNames[i] = GTL.contigsNames[i];
			}
		}
	}
	return *this;
}

int GlobalToLocal::searchContig(unsigned int GlobalCoordinate) {
	int min = 0;
	int max = contigs-1;
	int i;
	while (min <= max) {
		i = (min+max)/2;
		if (GlobalCoordinate > endPositions[i])
				min = i+1;
		else if (GlobalCoordinate < startPositions[i])
			max = i-1;
		else
			return i;
	}
	ERROR_CHANNEL << "Position \"" << GlobalCoordinate << "\" not found for GlobalToLocal::searchContig(int position)";
	std::exit(2);
	return 0;
}

bool GlobalToLocal::detectDuplicates() {
	set<string> temp;
	pair<set<string>::iterator,bool> ret;
	bool all_ok = true;
	for (int i = 0; i < contigs; i++) {
		ret = temp.insert(string(contigsNames[i]));
		if (ret.second == false) {
			ERROR_CHANNEL << "Error: name \"" << contigsNames[i] << "\" already exists!" << endl;
			all_ok = false;
		}
	}
	return all_ok;
}
