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

#ifndef MODULE_DCREATE_H_
#define MODULE_DCREATE_H_

#include "modules/Module.h"
using namespace options;

#include <mpi.h>

#include <vector>
#include <string>
using namespace std;

#include "io/Fasta.h"

namespace modules {

class Module_DCREATE : public Module {
public:
	Module_DCREATE() { proc_name = NULL;}
	virtual ~Module_DCREATE() { delete proc_name; remove_temporary_files(); }
	void execute(const Options & options);
	void remove_temporary_files();
private:
    int my_rank;
    int nprocs;
    char *proc_name;

    vector<string> temp_files;

    void compute_master(const Options & options);
    void send_sequences_to_slave(int node, int k, int blockLength, const string & input_filename, const string & output_filename) const;

    void receive_from_master() const;
    //void compute_hash(int k, int blockLength, const char * input_file, const char * output_file, bool methyl_hash) const ;
    void compute_hash(int k, int blockLength, const char * input_file, const char * output_file, bool methyl_hash) const;

    static bool sort_reverse_function (Fasta * i, Fasta * j) { return (i->length() > j->length()); }

};

}

#endif /* MODULE_DCREATE_H_ */
