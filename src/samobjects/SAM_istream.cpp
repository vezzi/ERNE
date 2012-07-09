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

#include "SAM_istream.h"

namespace samobjects {

SAM_istream::SAM_istream() : SAM() {
	end_of_file = false;
}

SAM_istream::SAM_istream(const char * filename) throw (SAM_IO_Error) : SAM() {
	open(filename);
}

SAM_istream::~SAM_istream() throw (SAM_IO_Error) {
	// Nothing to do
}

void SAM_istream::open(const char * _filename, bam_header_t * header) throw (SAM_IO_Error) {
	close();
	if (_filename != NULL)
		filename.assign(_filename);

	if (filename.size() == 0)
		throw SAM_IO_Error(SAM_IO_Error::empty_file_name);

	switch (detect_file_type_by_name(filename.c_str())) {
	case unknown_sam_format:
		throw SAM_IO_Error(SAM_IO_Error::incorrect_filename_format,filename.c_str());
		break;
	case bam_format:
		sam_file = samopen(filename.c_str(),"rb",NULL);
		break;
	case tam_format:
		sam_file = samopen(filename.c_str(),"r",NULL);
		break;
	case tam_compressed_format:
		throw SAM_IO_Error(SAM_IO_Error::not_supported_method,"tam compressed file");
		break;
	}
	if (sam_file == NULL) {
		throw SAM_IO_Error(SAM_IO_Error::cannot_open_file_r,filename.c_str());

	}
	popolate_name_list(sam_file->header);
	end_of_file = false;
	samfile_is_open = true;
}

void SAM_istream::close() {
	end_of_file = true;
	SAM::close();
}

bam1_t * SAM_istream::read() throw (SAM_IO_Error) {
	if (sam_file == NULL) {
		throw SAM_IO_Error(SAM_IO_Error::file_not_opened, "tried to read from a not previously opened file");
	}
	bam1_t * b = bam_init1();
	int bytes = samread(sam_file,b);
	if (bytes == -1)
		end_of_file = true;
	return b;
}

}
