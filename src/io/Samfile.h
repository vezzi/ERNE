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

#ifndef SAMFILE_H_
#define SAMFILE_H_

#include "data_structures/Mask.h"
#include "data_structures/Hash.h"

#include "options/Options.h"
using namespace options;

#include "samtools/sam.h"
#include <boost/thread/pthread/mutex.hpp>

class Samfile {
public:
	enum print_method_t { NO_ONE = 0, FIRST_ONLY = 1, SECOND_ONLY = 2, BOTH = 3};

	Samfile() { sam_file = NULL; CR = NULL; gui_output = false; insert_size_check = false;}
	virtual ~Samfile() { close_file(); }

#ifdef HAVE_MPI
	void open_file_mpi(const Options & options);
#endif
	void open_file(const Options & options, const Hash & H, const Hash & CR, int index = -1);
	void close_file();

	void print_output(bam1_t* b);
	void print_output(Mask &read);
	void print_output_paired(Mask &first, Mask &second, print_method_t print_method = BOTH);

	static bool file_exists(const char * filename);

	void set_insert_size(unsigned int min, unsigned int max) { insert_size_check = true; insert_size_min = min; insert_size_max = max; }

	void set_reading_group_id(string rgi) { reading_group_id = rgi; }
	string get_reading_group_id() const { return reading_group_id; }

private:

	mutex write_mutex;
	mutex gui_output_lock;

	samfile_t * sam_file;
	bam_header_t * sam_header;

	const Hash * CR;
	bool gui_output;

	string reading_group_id;

	bool insert_size_check;
	unsigned int insert_size_min;
	unsigned int insert_size_max;

	void print_bam_debug(bam1_t * b, const Mask & read);
	string header_description(const Options & options);
	static size_t aligned_region_size(bam1_t * bam);

};

#endif /* SAMFILE_H_ */
