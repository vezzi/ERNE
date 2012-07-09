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

#include "Options_DCREATE.h"

namespace options {

bool Options_DCREATE::process(int argc, char *argv[]) {

	this->argc = argc;
	this->argv = argv;

	// PROCESS PARAMETERS
	stringstream ss;
	ss << package_description() << endl << endl << "This is ERNE-DCREATE. Allowed options";
	po::options_description desc(ss.str().c_str());
	desc.add_options() ("help", "produce help message")
		("version", "print version and exit")
		("fasta", po::value< vector < string > >(), "input file (can be repeated several time)")
		("k", po::value<int>(), "k (default: autodetected)")
		("bl", po::value<int>(), "word size to be mapped in the structure (default: 15)")
		("nobalancing","no reference balancing while creating (useful for millions of fasta sequences)")
		("reference", po::value<string>(), "output reference file (in our format)")
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
	program_mode = program_dcreate;

	if (vm.count("k")) {
		auto_k = false;
		k = vm["k"].as<int>();
	}
	if (k > 15 or k < 1) {
		ERROR_CHANNEL << "--k parameter must be a positive number less or equal to 15" << endl;
		exit(1);
	}
	if (vm.count("bl"))
		blockLength = vm["bl"].as<int>();

	if (not(vm.count("reference") and vm.count("fasta"))) {
		ERROR_CHANNEL << "Both '--reference', and '--fasta' is required!" << endl;
		ERROR_CHANNEL  << "Try \"--help\" for help" << endl;
		return false;
	}

	balancing = not vm.count("nobalancing");

	input_files = vm["fasta"].as<vector < string > >();
	output_file = vm["reference"].as<string>();

	return true;
}

}
