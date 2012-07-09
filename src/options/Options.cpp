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

#include "Options.h"

namespace options {

Options::Options(const Options & orig) {
	program_mode = orig.program_mode;

	argc = orig.argc;
	argv = orig.argv;

	// input options
	reference_file = orig.reference_file;
	validation_files = orig.validation_files;
	auto_errors = orig.auto_errors;
	errors_rate = orig.errors_rate;
	common_errors_allowed = orig.common_errors_allowed;

	input_files = orig.input_files;

	fastqformat = orig.fastqformat;
	force_fastqformat = orig.fastqformat;

	contamination_check = orig.contamination_check;
	contamination_file = orig.contamination_file;

	quality_check = orig.quality_check;
	trim = orig.trim;
	min_phred_value_MOTT = orig.min_phred_value_MOTT;
	min_mean_quality = orig.min_mean_quality;
	min_size = orig.min_size;

	seed_sizes = orig.seed_sizes;
	seed_errors = orig.seed_errors;
	max_gap = orig.max_gap;
	gap = orig.gap;
	transcriptome = orig.transcriptome;

	insert_size_check = orig.insert_size_check;
	insert_size_min = orig.insert_size_min;
	insert_size_max = orig.insert_size_max;

	sample = orig.sample;
	bam_format = orig.bam_format;

	threads_number = orig.threads_number;

	query1 = orig.query1;
	query2 = orig.query2;
	paired_ends = orig.paired_ends;

	vectors_file = orig.vectors_file;

	auto_k = orig.auto_k;
	k = orig.k;
	blockLength = orig.blockLength;

	delta = orig.delta;
	indels = orig.indels;
	indels_max_value = orig.indels_max_value;

	// output options
	printAll = orig.printAll;
	toBePrinted = orig.toBePrinted;

	output_file = orig.output_file;

	gui_output = orig.gui_output;

	preserve_encoding = orig.preserve_encoding;
}

void Options::set_defaults() {
	program_mode = program_unknown;

	auto_errors = false;
	errors_rate = 15;
	common_errors_allowed = 0;

	coverage_threshold = 1;
	skip_alignment = false;
	alignment_only = false;
	only_unique = false;
	min_C_cov = 10;
	max_C_cov = 20;
	verbose=false;
	error_threshold = 1;

	force_fastqformat = false;

	contamination_check = false;

	quality_check = true;
	trim = true;
	min_phred_value_MOTT = 20;
	min_mean_quality = 20;
	min_size = 25;

	seed_sizes = 20;
	seed_errors = 1;
	max_gap = 100;
	gap = false;
	transcriptome = false;

	sample = string("no_sample_specified");
	bam_format = true;

	threads_number = 1;

	paired_ends = false;

	insert_size_check = false;

	auto_k = true;
	k = 12;
	blockLength = 15;

	delta = 0;
	indels = false;
	indels_max_value = 5;

	// output options
	printAll = false;
	toBePrinted = UINT_MAX;

	gui_output = false;

	preserve_encoding = false;

}

}
