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

#ifndef RESULTITEMGAPPED_H_
#define RESULTITEMGAPPED_H_

class ResultItemGapped {
public:
	int globalPosition1; // position in the concatenated reference
	int globalPosition2; // position in the concatenated reference
	int length1;
	int length2;
	int errors1; // number of errors
	int errors2; // number of errors
	bool strand; // strand : true=+, false=-
	int contig;
	ResultItemGapped(int GP1, int GP2, bool s, int l1, int l2, int e1, int e2, int c) {
		globalPosition1 = GP1;
		globalPosition2 = GP2;
		strand = s;
		length1 = l1;
		length2 = l2;
		errors1 = e1;
		errors2 = e2;
		contig = c;
	}

	struct less {
		bool operator()(const ResultItemGapped & hit1, const ResultItemGapped & hit2) {
			int p1 = hit1.globalPosition1 < hit1.globalPosition1 ? hit1.globalPosition1 : hit1.globalPosition2;
			int p2 = hit2.globalPosition1 < hit2.globalPosition2 ? hit2.globalPosition1 : hit2.globalPosition2;
			int e1 = hit1.errors1 + hit1.errors2;
			int e2 = hit2.errors1 + hit2.errors2;

			return ( e1<e2 or (e1==e2 and p1<p2) or (e1==e2 and p1==p2 and hit1.strand < hit2.strand)) ? true:false;
		}
	};

	struct equal {
		bool operator()(const ResultItemGapped & hit1, const ResultItemGapped & hit2) {
			return ( hit1.globalPosition1 == hit2.globalPosition1
					and hit1.globalPosition2 == hit2.globalPosition2
					and hit1.strand == hit2.strand) ? true : false;
		}
	};
};


#endif /* RESULTITEMGAPPED_H_ */
