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

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <string>
using namespace std;

#include <common.h>

#include "io/Fasta.h"
using namespace fasta;

#include "data_structures/Mask.h"

namespace options {

class Options {
public:
	Options() { set_defaults(); }
	virtual ~Options() { }
	Options(const Options & orig);

	virtual bool process(int argc, char *argv[]) = 0;

	enum program_mode_t { program_unknown, program_create, program_dcreate, program_dmap, program_map, program_bs5, program_filter};

	program_mode_t program_mode;

	int argc;
	char **argv;

	// input options
	string reference_file;
	vector<string> validation_files;
	bool auto_errors;
	t_errors errors_rate;
	t_errors common_errors_allowed;

	vector<string> input_files;

	Fasta::FASTQ_encoding fastqformat;
	bool force_fastqformat;

	bool contamination_check;
	string contamination_file;

	bool quality_check;
	bool trim;
	Mask::t_min_phred_value_MOTT min_phred_value_MOTT;
	Mask::t_min_phred_value_MOTT min_mean_quality;
	t_pattern_length min_size;

	t_pattern_length seed_sizes;
	t_errors seed_errors;
	t_length max_gap;
	bool gap;
	bool transcriptome;

	bool insert_size_check;
	unsigned int insert_size_min;
	unsigned int insert_size_max;

	string sample;
	bool bam_format;

	int threads_number;

	string query1;
	string query2;
	bool paired_ends;

	bool methyl_hash;
	bool skip_alignment;
	bool alignment_only;
	unsigned int coverage_threshold;
	bool only_unique;
	unsigned int min_C_cov;
	unsigned int max_C_cov;
	bool verbose;
	double error_threshold;

	string vectors_file;

	bool auto_k;
	int k;
	int blockLength;

	unsigned int delta;
	bool indels;
	int indels_max_value;

	// output options
	bool printAll;
	unsigned int toBePrinted;

	string output_file;

	bool gui_output;

	bool preserve_encoding;

#ifdef HAVE_MPI
	bool balancing;
#endif

protected:
	void set_defaults();

};

}

#endif /* OPTIONS_H_ */
