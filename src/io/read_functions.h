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

#ifndef READ_FUNCTIONS_H_
#define READ_FUNCTIONS_H_

#include "io/Auto_Unzip.h"
using namespace useful;

#include <boost/thread/pthread/mutex.hpp>
using namespace boost;

#include "data_structures/Mask.h"
#include "io/Fasta.h"
using namespace fasta;

int read_sequences(Auto_Unzip & input, int num_seq, Mask sequences[], Fasta::FASTQ_encoding format_type, bool gui_output);
int read_sequences(Auto_Unzip & first, Auto_Unzip & second, int num_seq, Mask sequences[], Fasta::FASTQ_encoding format_type, bool gui_output);

void output_progress(Auto_Unzip & input, bool gui_output);

#endif /* READ_FUNCTIONS_H_ */
