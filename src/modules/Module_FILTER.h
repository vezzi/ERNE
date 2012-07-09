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

#ifndef MODULE_FILTER_H_
#define MODULE_FILTER_H_

#include "modules/Module.h"
using namespace options;

#include "data_structures/Hash.h"
#include "data_structures/Mask.h"
#include "io/read_functions.h"

#include <boost/thread/thread.hpp>
#include <boost/thread/pthread/mutex.hpp>
using namespace boost;


namespace modules {

class Module_FILTER : public Module {
public:
	Module_FILTER() { }
	virtual ~Module_FILTER() { }

	void execute(const Options & options);

private:
	static void filter_reads(int id, Auto_Unzip *reads_fw, Auto_Unzip *reads_rv, const Options * options);
	static int read_sequences_conditional(Auto_Unzip * first, Auto_Unzip * second, int num_seq, Mask sequences[], Fasta::FASTQ_encoding format_type, bool gui_output, bool pe);
};

}

#endif /* MODULE_FILTER_H_ */
