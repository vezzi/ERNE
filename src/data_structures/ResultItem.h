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

#ifndef RESULTITEM_H_
#define RESULTITEM_H_



class ResultItem {
public:

	struct less {
		bool operator()(const ResultItem & hit1, const ResultItem & hit2) {
			int p1 = hit1.globalPosition;
			int p2 = hit2.globalPosition;
			int e1 = hit1.errors;
			int e2 = hit2.errors;

			return ( e1<e2 or (e1==e2 and p1<p2) or (e1==e2 and p1==p2 and hit1.strand < hit2.strand)) ? true:false;
		}
	};

	struct equal {
		bool operator()(const ResultItem & hit1, const ResultItem & hit2) {
			return ( hit1.globalPosition==hit2.globalPosition and hit1.strand == hit2.strand) ? true:false;
		}
	};

	unsigned int globalPosition; // position in the concatenated reference
	int errors; // number of errors
	bool strand; // strand : true=+, false=-
	bool indels;
	//	vector<pair<char, unsigned int> > cigar;


};

#endif /* RESULTITEM_H_ */
