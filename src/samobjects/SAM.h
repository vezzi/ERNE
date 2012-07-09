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
#ifndef SAM_H_
#define SAM_H_

#include <map>
#include <vector>
#include <string>
using namespace std;

#include <boost/regex.hpp>
using namespace boost;

#include "SAM_IO_Error.h"

#include "samtools/sam.h"

namespace samobjects {

class SAM {
public:
	enum sam_format_t { unknown_sam_format, tam_format, tam_compressed_format, bam_format };
	typedef map<const string, int32_t> name_to_tid_t;

	SAM();
	virtual ~SAM() throw (SAM_IO_Error);

	virtual void open(const char * filename = NULL, bam_header_t * header = NULL) throw (SAM_IO_Error) = 0;
	void revert() throw (SAM_IO_Error);
	void close();

	static sam_format_t detect_file_type_by_name(const char * filename);

	bool is_open() const { return samfile_is_open; }

	void popolate_name_list(bam_header_t * header);
	const string * tid_to_name(int32_t i) const;
	int32_t name_to_tid(const string * name) const;

	bam_header_t * header() const { return sam_file->header; }

	typedef vector<pair < string, uint32_t > > names_list_t;
	static bam_header_t * update_header_from_list(bam_header_t *header, names_list_t & list);

	void change_to_unsorted();
	bool check_ordered() const;

	static void regex_header(bam_header_t * header, const char * pattern, const char * replace);

protected:
	typedef vector<string> name_list_t;
	samfile_t * sam_file;
	string filename;
	bool samfile_is_open;

	name_list_t names_list;
	name_to_tid_t names_to_tid;


	void unset();
	void clear() throw (SAM_IO_Error);


};

}

#endif /* SAM_H_ */
