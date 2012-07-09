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

#ifndef MODULE_BS5_H
#define MODULE_BS5_H

#include "modules/Module.h"
#include <data_structures/Hash.h>
#include <data_structures/GlobalToLocal.h>
#include <options/Options.h>

#include "samtools/bam.h"
#include "samtools/sam.h"
#include "samobjects/BAM1.h"
#include "io/Samfile.h"

#include <vector>
#include <string>

#define CG "CG";
#define CHG "CHG";
#define CHH "CHH";

#include "samobjects/BAM1.h"
using namespace samobjects;

using namespace std;

namespace modules {

class Module_BS5 {

public:
	Module_BS5();
	virtual ~Module_BS5() { }

	void execute(const Options & options);

	struct score{

			double methyl_distance;	//distance between read and reference methylation pattern
			double cov_avg;			//coverage average of the read's positions in the text
			int C_covered;			//number of C covered

			score(double md, double ca, int fc){

				this->methyl_distance = md;
				this->cov_avg = ca;
				this->C_covered = fc;
			}

			bool operator <=(score s){

				if((unsigned int)(methyl_distance)>(unsigned int)(s.methyl_distance))
					return false;
				if((unsigned int)(methyl_distance)==(unsigned int)(s.methyl_distance))
					if(cov_avg>s.cov_avg)
						return false;

				return true;

			}
		};

private:

	Hash* H;
	Hash* CR;
	Samfile samfile_output;//final samfile.bam

	static const int cg = 0;
	static const int chg = 1;
	static const int chh = 2;

	//parameters:
	double C_cov;//min number of C to be covered

	unsigned int uniquely_aligned;

	unsigned char* methylConsensus_T;
	unsigned char* methylConsensus_C;
	unsigned char* mult_coverage;//number of C/T aligned under a C using ONLY multiple-mapping reads

	unsigned int cycle;//cycle counter

	double error_threshold;

	unsigned long int totCoverage;

	static const bool FW = true;
	static const bool RC = false;

	static const int MULTIPLE = -1;
	static const int NO_ALIGNMENTS = -2;

	void align_and_update_profile(bam1_t *b);
	void align_and_update_profile(bam1_t ** ali1, int IH1, bam1_t** ali2, int IH2);
	int check_alignments(bam1_t **ali1, int IH);

	score *methyl_distance(bam1_t *b);

	double alpha(char c,unsigned int shift);

	unsigned int global_coord(unsigned int tid, unsigned int pos);
	int context_in_text(unsigned long i);

	unsigned int covg(unsigned int i);

	unsigned int C_coverage(unsigned int i);
	unsigned int T_coverage(unsigned int i);

	double methyl_level(unsigned int i);

	void increment_C(unsigned int i);
	void increment_T(unsigned int i);
	void increment_mult_coverage(unsigned int i);

	void validation(const Options & options);

};

}
#endif /* MODULE_BS5_H */
