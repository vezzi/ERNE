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

#ifndef SAM_SEARCH_H_
#define SAM_SEARCH_H_

#include "SAM.h"
#include "BAM1.h"

#include <list>
using namespace std;

namespace samobjects {

class SAM_search : public SAM {
public:
	typedef list<BAM1> bam_set_t;

	SAM_search();
	SAM_search(const char * name);
	virtual ~SAM_search() throw (SAM_IO_Error);

	void open(const char * filename = NULL, bam_header_t * header = NULL) throw (SAM_IO_Error);
	void close();

	BAM1 * mate(const BAM1 *bam) { return mate(bam->pointer()); }
	BAM1 * mate(const BAM1 &bam) { return mate(bam.pointer()); }
	BAM1 * mate(const bam1_t &bam) { return mate(&bam); }
	BAM1 * mate(const bam1_t *bam);

	bam_set_t * overlap(BAM1 * bam);
	void overlap(bam_set_t & l, BAM1 * bam);
	bam_set_t * overlap(int tid, int beg, int end);
	void overlap(bam_set_t & l, int tid, int beg, int end);

protected:

	bam_index_t * index;

	static bool file_exists(string & strFilename);

	void create_index() throw (SAM_IO_Error);
	void load_index() throw (SAM_IO_Error);

	int compare_function(const bam1_t *b, void *data);

};

}

#endif /* SAM_SEARCH_H_ */
