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

#ifndef GENERICITEM_H_
#define GENERICITEM_H_

//#include "ResultItem.h"
//#include "ResultItemGapped.h"

class GenericItem {
public:
	GenericItem() { gapped = false; }
	GenericItem(const ResultItem & item) { set(item); }
	GenericItem(const ResultItemGapped & item) { set(item); }
	virtual ~GenericItem() { }

	bool strand; // strand : true=+, false=-
	bool indels;

	bool gapped;

	unsigned int globalPosition1;
	unsigned int globalPosition2;
	unsigned int length1;
	unsigned int length2;
	int errors1;
	int errors2;
	int contig;

	int errors() const { return (gapped ? (errors1 + errors2) : errors1); }
	unsigned int min_position() const { return (gapped ? (globalPosition1 < globalPosition2 ? globalPosition1 : globalPosition2) : globalPosition1); }

	void set(const ResultItem & item) {
		gapped = false;
		globalPosition1 = globalPosition2 = item.globalPosition;
		indels = item.indels;
		strand = item.strand;
		errors1 = item.errors;
		errors2 = length1 = length2 = 0;
	}
	void set(const ResultItemGapped & item) {
		gapped = true;
		globalPosition1 = item.globalPosition2;
		globalPosition2 = item.globalPosition2;
		strand = item.strand;
		errors1 = item.errors1;
		errors2 = item.errors2;
		length1 = item.length1;
		length2 = item.length2;
		contig = item.contig;
	}

	struct less {
		bool operator()(const GenericItem & hit1, const GenericItem & hit2) {
			int p1 = hit1.globalPosition1 < hit1.globalPosition1 ? hit1.globalPosition1 : hit1.globalPosition2;
			int p2 = hit2.globalPosition1 < hit2.globalPosition2 ? hit2.globalPosition1 : hit2.globalPosition2;
			int e1 = hit1.errors1 + hit1.errors2;
			int e2 = hit2.errors1 + hit2.errors2;

			return ( e1<e2 or (e1==e2 and p1<p2) or (e1==e2 and p1==p2 and hit1.strand < hit2.strand)) ? true:false;
		}
	};

	struct equal {
		bool operator()(const GenericItem & hit1, const GenericItem & hit2) {
			return ( hit1.globalPosition1 == hit2.globalPosition1
					and hit1.globalPosition2 == hit2.globalPosition2
					and hit1.strand == hit2.strand) ? true:false;
		}
	};
};

#endif /* GENERICITEM_H_ */
