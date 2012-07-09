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
#ifndef SAM_IO_ERROR_H_
#define SAM_IO_ERROR_H_

#include <exception>
#include <string>
#include <sstream>
using namespace std;

namespace samobjects {

class SAM_IO_Error : public std::exception {
public:
	enum error_description { generic_error, file_not_found, incorrect_filename_format,
		empty_file_name, not_supported_method, incorrect_bam1_t, file_not_opened, cannot_open_file_r, cannot_open_file_rw };
	SAM_IO_Error();
	SAM_IO_Error(error_description error_type,const char * comment_description = NULL);
	virtual ~SAM_IO_Error() throw();
    const virtual char *what() const throw ();

    const char * getComment() const { return comment.c_str(); }
    error_description getError() const { return e; }

protected:
	string comment;
	error_description e;
	string description;

	void prepare();
};

}

#endif /* SAM_IO_ERROR_H_ */
