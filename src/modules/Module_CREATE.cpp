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

#include "Module_CREATE.h"

#include "data_structures/Hash.h"

namespace modules {

void Module_CREATE::execute(const Options & options) {
	clock_t start = clock();
	int k;
	if (options.auto_k) {
		VERBOSE_CHANNEL << "Reading file for calculating k ..." << endl;
		k = Hash::calculate_k(options.input_files,options.methyl_hash);
		VERBOSE_CHANNEL << "Detected optimal value for k: " << k << endl;
	}else
		k = options.k;
	Hash * H = new Hash(k, options.blockLength,options.methyl_hash);
	H->createHASH(options.input_files);
	clock_t end = clock();
	VERBOSE_CHANNEL << "building time = " <<((end - start) / (double)CLOCKS_PER_SEC) << "\n";
	H->save(options.output_file.c_str());
	clock_t end2 = clock();
	VERBOSE_CHANNEL << "saving time = " <<((end2 - end) / (double)CLOCKS_PER_SEC) << "\n";
	VERBOSE_CHANNEL << "total time = " <<((end2 - start) / (double)CLOCKS_PER_SEC) << "\n";

}

}
