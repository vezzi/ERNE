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

#ifndef MODULE_MAP_H_
#define MODULE_MAP_H_

#include "modules/Module.h"

#include "io/read_functions.h"
#include "io/Samfile.h"
#include "io/formatting.h"

#include "data_structures/GenericItem.h"

#include "io/Fasta.h"
using namespace fasta;

#include <boost/thread/thread.hpp>
using namespace boost;

using namespace std;

namespace modules {

#define SEQUENCES_FOR_BLOCK 1000
typedef boost::mt19937 RNGType;

class Module_MAP : public Module {
public:

	typedef vector<GenericItem> GenericItems;

	Module_MAP() : Module() { }
	virtual ~Module_MAP() { }

	virtual void execute(const Options & options);

	Hash H;
	Hash CR;
	void free_hash_memory();

	string get_reading_group_id() const { return output_samfile.get_reading_group_id(); }

protected:

	static RNGType rng;
	static boost::uniform_real<> double_range;
	static boost::variate_generator<RNGType&, boost::uniform_real<> > guessRes;

	Samfile output_samfile;

	bool force_fastqformat;
	Fasta::FASTQ_encoding fastqformat;

	//unsigned int blockLength;
	bool verbose;
	int threads_number;
	bool gui_output;

	bool quality_check;
	bool trim;
	Mask::t_min_phred_value_MOTT min_phred_value_MOTT;
	Mask::t_min_phred_value_MOTT min_mean_quality;
	t_pattern_length min_size;
	bool auto_errors;
	t_errors common_errors_allowed;
	t_errors errors_rate;
	bool contamination_check;

	unsigned long int processed;

	bool gap;
	bool methyl_reconstruction;
	t_pattern_length seed_sizes;
	t_errors seed_errors;
	t_length max_gap;
	bool transcriptome;

	bool printAll;
	unsigned int toBePrinted;

	unsigned int delta;
	bool indels;
	int indels_max_value;
	
	string query1;
	string query2;

	bool insert_size_check;
	unsigned int insert_size_min;
	unsigned int insert_size_max;

	void initialize_parameters(const Options & options);

	void process_single_reads();
	static void process_single_reads_thr(Module_MAP * search, int id, Auto_Unzip * input);

	void process_pair_reads();
	static void process_pair_reads_thr(Module_MAP * search, int id, Auto_Unzip * inputFile1, Auto_Unzip * inputFile2);
};

}

#endif /* MODULE_MAP_H_ */
