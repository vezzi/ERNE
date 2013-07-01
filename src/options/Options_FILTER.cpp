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

#include "Options_FILTER.h"

namespace options {



bool Options_FILTER::process(int argc, char *argv[]) {

	this->argc = argc;
	this->argv = argv;

	// PROCESS PARAMETERS
	stringstream ss;
	ss << package_description() << endl << endl << "This is ERNE-FILTER. Allowed options";
	po::options_description desc(ss.str().c_str());
	desc.add_options() ("help", "produce help message")
		("version", "print version and exit")
		("reference", po::value<string>(), "reference file to use (in our format)")
		("query1", po::value<string>(), "query1 file (can be compressed with gzip or bzip2, or a pipe)")
		("query2", po::value<string>(), "query2 file (can be compressed with gzip or bzip2, or a pipe)")
		("output", po::value<string>(), "output prefix")
		("force-illumina","force ILLUMINA 1.3+ FASTQ format (default: auto-detect)")
		("force-standard","force standard SANGER FASTQ format (default: auto-detect)")
		("threads", po::value<unsigned int>(), "maximum number of allowed threads (default 1)")
		("auto-errors","use automatically one error every ~15bp")
		("errors", po::value<t_errors>(), "fixed number of errors allowed (>= 0, default 0) --- overwrites auto-errors ")
		("no-auto-trim","disable automatic trimming (read are trimmed and aligned but only original reads are returned)")
		("min-phred-value-mott",po::value<Mask::t_min_phred_value_MOTT>(),"minimum value used by Mott-like trimming (default 20)")
		("min-mean-phred-quality",po::value<Mask::t_min_phred_value_MOTT>(),"minimum mean value to accept a (trimmed) sequence (default 20)")
		("min-size", po::value<t_pattern_length>(), "minimum sequence length after trimming (default 25)")
		("preserve-encoding","preserve input encoding")
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
	program_mode = program_filter;
	if (not vm.count("output")) {
		ERROR_CHANNEL << "--output is required" << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		exit(1);
	}
	if (not(vm.count("query1") and vm.count("query2")) ) {
		ERROR_CHANNEL << "both --query1 and --query2 are required!" << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		exit(1);
	}
	query1 = vm["query1"].as<string>();
	query2 = vm["query2"].as<string>();

	if (vm.count("min-phred-value-CLC"))
		min_phred_value_MOTT = vm["min-phred-value-CLC"].as<Mask::t_min_phred_value_MOTT>();

	if (vm.count("min-mean-phred-quality"))
		min_mean_quality = vm["min-mean-phred-quality"].as<Mask::t_min_phred_value_MOTT>();

	if (vm.count("min-size"))
		min_size = vm["min-size"].as<t_pattern_length>();

	if (vm.count("no-auto-trim"))
		trim = false;

	if (vm.count("force-illumina")) {
		force_fastqformat = true;
		fastqformat = Fasta::illumina_fastq_encoding;
	}

	if (vm.count("force-standard")) {
		force_fastqformat = true;
		fastqformat = Fasta::standard_fastq_encoding;
	}

	if (vm.count("preserve-encoding"))
		preserve_encoding = true;

	if (vm.count("reference"))
		reference_file = vm["reference"].as<string>();

	output_file = vm["output"].as<string>();

	if (vm.count("threads"))
		threads_number = vm["threads"].as<unsigned int>();

	auto_errors = true; // default

	if (vm.count("errors")) {
		common_errors_allowed = vm["errors"].as<t_errors>();
		auto_errors = false;
	}

	return true;
}

}
