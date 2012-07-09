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

#include "SAM_search.h"

#include <sys/stat.h>
#include <utility>

namespace samobjects {

SAM_search::SAM_search() {
	index = NULL;
}

SAM_search::SAM_search(const char * name) {
	index = NULL;
	open(name);
}

SAM_search::~SAM_search() throw (SAM_IO_Error) {
	close();
}

void SAM_search::close() {
	if (index != NULL)
		bam_index_destroy(index);
	SAM::close();
}

/* static */ bool SAM_search::file_exists(string & strFilename) {
	struct stat stFileInfo;
	// Attempt to get the file attributes
	if (stat(strFilename.c_str(),&stFileInfo) == 0)
		return true;
	else
		return false;
}

void SAM_search::open(const char * _filename, bam_header_t * header) throw (SAM_IO_Error) {
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
		throw SAM_IO_Error(SAM_IO_Error::not_supported_method,"tam file");
		//sam_file = samopen(filename.c_str(),"r",NULL);
		break;
	case tam_compressed_format:
		throw SAM_IO_Error(SAM_IO_Error::not_supported_method,"tam compressed file");
		break;
	}

	if (sam_file == NULL) {
		throw SAM_IO_Error(SAM_IO_Error::cannot_open_file_r,filename.c_str());

	}

	//TODO: controllare apertura file

	//TODO: controllare!!!
	/*
	// check for ordered field
	if (not check_ordered()) {
		close();
		throw SAM_IO_Error(SAM_IO_Error::not_supported_method,"bam compressed file is not sorted");
	}
	*/

	load_index();

	popolate_name_list(sam_file->header);
	samfile_is_open = true;

}

void SAM_search::create_index() throw (SAM_IO_Error) {
	 if (bam_index_build(filename.c_str()) != 0)
		 throw SAM_IO_Error(SAM_IO_Error::generic_error,string("problem creating index (.bai) file for ").append(filename).c_str());
}

void SAM_search::load_index() throw (SAM_IO_Error) {
	if (not file_exists(string(filename).append(".bai")))
		create_index();

	index = bam_index_load(filename.c_str());
	if (index == NULL)
		throw SAM_IO_Error(SAM_IO_Error::generic_error,string("problem loading index (.bai) file for ").append(filename).c_str());
}

typedef pair<const bam1_t *,BAM1 *> bam_input_e_output_t;
typedef pair<bool,bam_input_e_output_t> bool_and_bams_t;

int search_mate_function(const bam1_t *bam, void *data) {
	bool_and_bams_t * info = (bool_and_bams_t *)data;
	if (info->first == false) {
		uint32_t bam_flag = bam->core.flag;
		const bam1_t * search = info->second.first;
		uint32_t search_flag = search->core.flag;
		if ( (bam_flag & BAM_FPAIRED)
			and
				( ((search_flag & BAM_FREAD1) and (bam_flag & BAM_FREAD2))
						or
				  ((search_flag & BAM_FREAD2) and (bam_flag & BAM_FREAD1)) )
			and
			(strcmp(bam1_qname(search),bam1_qname(bam)) == 0) ) {
				//bam1_t * mate = info->second.second->pointer();
				// FOUND!
			info->second.second->copy(bam);
			info->first = true;
		}
	}
	return 0;
}

BAM1 * SAM_search::mate(const bam1_t *bam) {
	BAM1 * mate = new BAM1();
	if (bam != NULL and (bam->core.flag & BAM_FPAIRED) and not (bam->core.flag & BAM_FMUNMAP)) {
		bool_and_bams_t info(false,bam_input_e_output_t(bam,mate));

		bam_fetch(sam_file->x.bam, index, bam->core.mtid, bam->core.mpos, bam->core.mpos+1,
			&info,&search_mate_function);
	}
	return mate;
}

SAM_search::bam_set_t * SAM_search::overlap(BAM1 * bam) {
	bam_set_t * l = new bam_set_t;
	overlap(*l, bam->core()->tid, bam->get_pos(), bam->get_end());
	return l;
}

void SAM_search::overlap(bam_set_t & l, BAM1 * bam) {
	overlap(l, bam->core()->tid, bam->core()->pos, bam->core()->pos+bam->aligned_region_size());
}

SAM_search::bam_set_t * SAM_search::overlap(int tid, int beg, int end) {
	bam_set_t * l = new bam_set_t;
	overlap(*l, tid, beg, end);
	return l;
}

int overlap_function(const bam1_t *bam, void *data) {
	SAM_search::bam_set_t * l = (SAM_search::bam_set_t *)data;
	l->push_back(BAM1());
	l->rbegin()->copy(bam);
	return 0;
}

void SAM_search::overlap(bam_set_t & l, int tid, int beg, int end) {
	bam_fetch(sam_file->x.bam, index, tid, beg, end, &l,&overlap_function);
}


}
