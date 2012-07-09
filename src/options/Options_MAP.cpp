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

#include "Options_MAP.h"

namespace options {



bool Options_MAP::process(int argc, char *argv[]) {

	this->argc = argc;
	this->argv = argv;

	// PROCESS PARAMETERS
	stringstream ss;
	ss << package_description() << endl << endl << "This is ERNE-MAP. Allowed options";
	po::options_description desc(ss.str().c_str());
	desc.add_options() ("help", "produce help message")
		("version", "print version and exit")
		("validation", po::value< vector < string > >(), "fasta validation files for methylation check")
		("methyl-hash", "in --create create reference for methylation search. This reference can be used also for standard search. The opposite is not true: if a reference is created without this option, it cannot be used for methylation search.")
		("only-unique", "during methylation reconstruction use only unique alignments")
		("use-all-alignments", "during methylation reconstruction with skip-alignment use all the alignments from the input bam (ignore IH,HI,NH)")
		("skip-alignment", "during methylation reconstruction skip alignment phase; a pre-created bam with search must be present with extension _0")
		("max-methyl-error", po::value<double>(), "max methylation error expressed in % allowed for an alignment. Default: 1%")
		("reference", po::value<string>(), "reference file to use (in our format)")
		("contamination-reference", po::value<string>(), "reference file to use for contamination check (in ERNE format)")
		("query1", po::value<string>(), "query1 file (can be compressed with gzip or bzip2, or a pipe)")
		("query2", po::value<string>(), "query2 file (can be compressed with gzip or bzip2, or a pipe)")
//		("vectors-file",po::value<string>(),"vector file for contamination trimming")
		("output", po::value<string>(), "SAM/BAM output file (default: SAM format)")
		("sam", "output file in SAM format instead of BAM format")
		("gui-output","enable output information for GUI")
		("force-illumina","force ILLUMINA 1.3+ FASTQ format (default: auto-detect)")
		("force-standard","force standard SANGER FASTQ format (default: auto-detect)")
		("threads", po::value<unsigned int>(), "maximum number of allowed threads (default 1)")
		("auto-errors","use automatically one error every ~15bp")
		("errors-rate",po::value<t_errors>(),"change automatically error rate (default 15)")
		("errors", po::value<t_errors>(), "errors allowed (>= 0, default 0)")
		("delta", po::value<unsigned int>(), "DELTA value (default 0)")
		("indels", "allow indels in read alignment")
		("indels-max",po::value<int>(),"max base pairs indels value (default: 5)")
		("insert-size-min", po::value<unsigned int>(), "minimum insertion size for proper pair (default: none, if --insert-size-max is defined, it is optional and default is 0)")
		("insert-size-max", po::value<unsigned int>(), "maximum insertion size for proper pair (default: none, required if --insert-size-min is defined)")
		("sample",po::value<string>(), "sample name")
		("no-auto-trim","disable automatic trim")
		("min-phred-value-mott",po::value<Mask::t_min_phred_value_MOTT>(),"minimum value used by Mott-like trimming (default 20)")
		("min-mean-phred-quality",po::value<Mask::t_min_phred_value_MOTT>(),"minimum mean value to accept a (trimmed) sequence (default 20)")
		("min-size", po::value<t_pattern_length>(), "min length for a sequence (default 25)")
		("gap",po::value<t_length>(),"Efficiently search 1 gap inside an interval of arg bases")
		("transcriptome","detect correct intron-exons junctions in transcriptome-seq (require --gap)")
//		("seed-sizes", po::value<t_pattern_length>(), "seed sizes for a sequence (default 10)")
//		("seed-errors", po::value<t_errors>(), "seed errors for a sequence (default 1)")
//		("no-quality-check","sequences may not be dropped using quality information")
		("print-all","print all possible alignments [only for single reads]")
		("print-first", po::value<int>(), "print the number of specified alignments [only for single reads]")
		;

	po::variables_map vm;
	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	} catch (boost::program_options::error & error) {
		ERROR_CHANNEL <<  error.what() << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		exit(2);
	}

	if (vm.count("help")) {
		DEFAULT_CHANNEL << desc << endl;
		exit(0);
	}

	if (vm.count("version")) {
		DEFAULT_CHANNEL << package_description() << endl;
		exit(0);
	}

	// Process common parameters
	program_mode = program_map;

	if (not(vm.count("reference") )) {
		ERROR_CHANNEL << "'--reference' is required" << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		return false;
	}

	if (!(vm.count("query1")) and !(vm.count("query1") and vm.count("query2")) ) {
		ERROR_CHANNEL << "At least one --query1 (single read) or a pair --query1 --query2 (pair ends) are required!" << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		return false;
	}

	if (vm.count("no-auto-trim"))
		trim = false;

	if (vm.count("sam"))
		bam_format = false;

	if (vm.count("sample"))
		sample = vm["sample"].as<string>();

	if (vm.count("query1"))
		query1= vm["query1"].as<string>();

	if (vm.count("query2"))
		query2 = vm["query2"].as<string>();

	if (vm.count("gui-output"))
		gui_output = true;

	if (!vm.count("output")) {
		ERROR_CHANNEL << "--output parameter is required!" << endl;
		return false;
	}

	if (vm.count("force-illumina") + vm.count("force-standard") > 1) {
		ERROR_CHANNEL << "At most one between --force-illumina --force-standard is allowed!" << endl;
		return false;
	}

	if (vm.count("force-illumina")) {
		force_fastqformat = true;
		fastqformat = Fasta::illumina_fastq_encoding;
	}
	if (vm.count("force-standard")) {
		force_fastqformat = true;
		fastqformat = Fasta::standard_fastq_encoding;
	}

	if (vm.count("threads"))
		threads_number = vm["threads"].as<unsigned int>();

	if (vm.count("auto-errors"))
		auto_errors = true;

	if (vm.count("errors-rate"))
		errors_rate = vm["errors-rate"].as<t_errors>();

	if (vm.count("errors"))
		common_errors_allowed = vm["errors"].as<t_errors>();

	if (vm.count("no-seed") and vm.count("improved")) {
		ERROR_CHANNEL << "At most one between '--no-seed' adn '--improved' is allowed!" << endl;
		ERROR_CHANNEL  << "Try \"--help\" for help" << endl;
		return false;
	}

	if (vm.count("insert-size-max")) {
		insert_size_check = true;
		insert_size_min = 0;
		insert_size_max = vm["insert-size-max"].as<unsigned int>();
	}

	if (vm.count("insert-size-min")) {
		if (not insert_size_check) {
			ERROR_CHANNEL << "--insert-size-max is required if --insert-size-min is specified!" << endl;
			exit(1);
		}
		insert_size_min = vm["insert-size-min"].as<unsigned int>();
	}

	if (vm.count("min-size"))
		min_size = vm["min-size"].as<t_pattern_length>();


	if (vm.count("seed-sizes"))
		seed_sizes = vm["seed-sizes"].as<t_pattern_length>();

	if (vm.count("seed-errors"))
		seed_errors = vm["seed-errors"].as<t_errors>();

	if (vm.count("no-quality-check"))
		quality_check = false;

	if (vm.count("min-phred-value-mott"))
		min_phred_value_MOTT = vm["min-phred-value-mott"].as<Mask::t_min_phred_value_MOTT>();

	if (vm.count("min-mean-phread-quality"))
		min_mean_quality = vm["min-mean-phread-quality"].as<Mask::t_min_phred_value_MOTT>();

	if (vm.count("print-all"))
		printAll = true;

	if (vm.count("print-first")) {
		printAll = true;
		toBePrinted = vm["print-first"].as<int>();
	}

	reference_file = vm["reference"].as<string>();

	if (vm.count("contamination-reference")) {
		contamination_check = true;
		contamination_file = vm["contamination-reference"].as<string>();
		if (reference_file.compare(contamination_file.c_str()) == 0) {
			ERROR_CHANNEL << "--reference and --contamination-reference cannot be equal" << endl;
			ERROR_CHANNEL << "Try \"--help\" for help" << endl;
			exit(1);
		}
	}

	if (vm.count("gap")) {
		gap = true;
		max_gap = vm["gap"].as<t_length>();
	}

	if (vm.count("transcriptome")) {
		if (not gap) {
			ERROR_CHANNEL << "--transcriptome require --gap parameter!!" << endl;
			ERROR_CHANNEL << "Try \"--help\" for help" << endl;
			exit(1);
		}
		transcriptome  = true;
	}

	if (vm.count("delta"))
		delta = vm["delta"].as<unsigned int>();

	if (vm.count("indels")) {
		indels = true;
		if (vm.count("indels-max")) {
			indels_max_value = vm["indels-max"].as<int>();
		}
	}

	output_file = vm["output"].as<string>();

	if (vm.count("query1") and vm.count("query2"))
		paired_ends = true;

	return true;
}

}
