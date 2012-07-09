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

#ifndef SAM_OSTREAM_H_
#define SAM_OSTREAM_H_

#include "SAM.h"

namespace samobjects {

class SAM_ostream : public SAM {
public:
	SAM_ostream();
	SAM_ostream(const char * filename, bam_header_t * header) throw (SAM_IO_Error);
	virtual ~SAM_ostream() throw (SAM_IO_Error);
	virtual void open(const char * filename = NULL, bam_header_t * header = NULL) throw (SAM_IO_Error);

	void write(bam1_t * bam) throw (SAM_IO_Error);
};

}

#endif /* SAM_OSTREAM_H_ */
