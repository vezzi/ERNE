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

#include "SAM_IO_Error.h"

namespace samobjects {

SAM_IO_Error::SAM_IO_Error() {
	e = generic_error;
	prepare();
}

SAM_IO_Error::~SAM_IO_Error() throw() {
	// Nothing to do
}

SAM_IO_Error::SAM_IO_Error(error_description error_type, const char *comment_description) {
	e = error_type;
	if (comment_description != NULL)
		comment.assign(comment_description);
	prepare();
}

void SAM_IO_Error::prepare() {
	switch (e) {
	case generic_error:
		description.assign("Generic error");
		break;
	case file_not_found:
		description.assign("File not found");
		break;
	case incorrect_filename_format:
		description.assign("Incorrect filename format");
		break;
	case empty_file_name:
		description.assign("Empty file name");
		break;
	case not_supported_method:
		description.assign("Not supported method");
		break;
	case incorrect_bam1_t:
		description.assign("Incorrect BAM struct");
		break;
	case file_not_opened:
		description.assign("File not opened");
		break;
	case cannot_open_file_r:
		description.assign("Cannot open file for input");
		break;
	case cannot_open_file_rw:
		description.assign("Cannot open file for output");
		break;
	default:
		description.assign("Unknown error");
		break;
	}
	if (comment.size() > 0) {
		description.append(" : ");
		description.append(comment);
	}
}

const char* SAM_IO_Error::what() const throw() {
	return description.c_str();
}

}

