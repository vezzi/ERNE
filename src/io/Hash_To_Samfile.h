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

#ifndef HASH_TO_SAMFILE_H_
#define HASH_TO_SAMFILE_H_

#include "data_structures/GlobalToLocal.h"
#include "samtools/sam.h"
#include "common.h"
#include <map>

class Hash_To_Samfile {
public:

	typedef int32_t index_t;

	Hash_To_Samfile();
	virtual ~Hash_To_Samfile();
	void link(const GlobalToLocal & gtl, const string header_file);

	int32_t convert(int i);
protected:
	index_t * index;
	int max;

	void clear();
	static bool file_exists(const char * filename);
};

#endif /* HASH_TO_SAMFILE_H_ */
