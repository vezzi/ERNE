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

#include "Hash_To_Samfile.h"

#include <sys/stat.h>

Hash_To_Samfile::Hash_To_Samfile() {
	index = NULL;
	max = 0;
}

Hash_To_Samfile::~Hash_To_Samfile() {
	clear();
}

void Hash_To_Samfile::clear() {
	if (index != NULL) {
		delete [] index;
		index = NULL;
	}
}

void Hash_To_Samfile::link(const GlobalToLocal & gtl, const string header_file) {
	if (not file_exists(header_file.c_str())) {
		ERROR_CHANNEL << "=== Cannot open file " << header_file << endl;
		ERROR_CHANNEL << "=== Check parameters! " << endl;
		exit(3);
	}

	bam_header_t * sam_header = sam_header_read2(header_file.c_str());

	clear();
	max = gtl.contigs;
	if (max == 0) {
		ERROR_CHANNEL << "There are no contigs in Hash_To_Samfile::link(const GlobalToLocal & gtl, const string header_file)" << endl;
		exit(4);
	}
	index = new index_t[max];

	map<string,int32_t> conversion;

	for (int32_t i = 0; i < sam_header->n_targets; i++)
		conversion.insert(make_pair(string(sam_header->target_name[i]),i));

	for (int i = 0; i < max; i++) {
		map<string,int32_t>::iterator iter = conversion.find(gtl.contigsNames[i]);
		if (iter == conversion.end()) {
			ERROR_CHANNEL << "=== reference " << gtl.contigsNames[i] << " not found in file " << header_file << endl;
			exit(4);
		}
		index[i] = iter->second;
	}

	bam_header_destroy(sam_header);
}

int32_t Hash_To_Samfile::convert(int i) {
	if (i < 0 or i >= max) {
		ERROR_CHANNEL << "Value " << i << " out of range: [0," << (max-1) << "] for Hash_To_Samfile::convert(int i)" << endl;
		std::exit(2);
	}
	return index[i];
}


/* static */
bool Hash_To_Samfile::file_exists(const char * filename) {
	struct stat stFileInfo;
	return stat(filename, &stFileInfo) == 0;
}
