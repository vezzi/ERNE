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

#ifndef MODULE_DMAP_H_
#define MODULE_DMAP_H_

#include "mpi.h"

#include "modules/Module_MAP.h"

#include "io/Hash_To_Samfile.h"

#include <boost/thread/thread.hpp>
using namespace boost;

#define MAX_NODE 16

namespace modules {

class Module_DMAP : public Module_MAP {
public:
	typedef pair<Mask *, int> Transmitting_Result;
	typedef pair<bool, int> Random_Choice_Result;

	Module_DMAP() : Module_MAP() {
		finished = false;
		proc_name = NULL;
	}

	virtual ~Module_DMAP() {
		delete proc_name;
	}

	void execute(const Options & options);

private:

	int my_rank;
	int nprocs;
	int nworkers;
	char *proc_name;

	MPI::Win win;

	//mutex thread_mutex;
    mutex mpi_mutex; // TODO rinominare

	Hash_To_Samfile contig_conversion;

	bool finished;

	Auto_Unzip input_file1;
	Auto_Unzip input_file2;

	void execute_generic_worker(const Options & options);
	void send_to_next(const vector<Mask> &printable_solutions,int id);
	Module_DMAP::Transmitting_Result receive_from_previous(int id);
	Transmitting_Result recv_output(int node, int id);
	void send_output(const vector<Mask> &printable_solutions, int node, int id);
	static Random_Choice_Result random_choice_from_previous(int previous, int actual);
	static void generic_worker_single_thr(Module_DMAP * search, int id);
	static void generic_worker_paired_thr(Module_DMAP * search, int id);


	//static Random_Choice_Result random_choice_from_previous(int previous, int actual);
	//void execute_worker(const Options & options);
	//void execute_master(const Options & options);

	//static void worker_single_thr(Module_DMAP * search, int id, Mask::MaskData* sequences, int numberOfSequences, MPI::Win &win);
	//static void worker_paired_thr(Module_DMAP * search, int id, Mask::MaskData* sequences, int numberOfSequences, MPI::Win &win);

	//void UpdateRMA_mutex_lock(MPI::Win &win);
	//void UpdateRMA_mutex_unlock(MPI::Win &win);


};

}

#endif /* MODULE_DMAP_H_ */
