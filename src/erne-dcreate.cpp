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

using namespace std;

#include "options/Options_DCREATE.h"
using namespace options;

#include "modules/Module_DCREATE.h"
using namespace modules;

int main(int argc, char *argv[]) {
	Options_DCREATE options(argc,argv);

	MPI::Init_thread(argc,argv, MPI::THREAD_MULTIPLE);

    if (MPI::COMM_WORLD.Get_size() <= 1) {
        DEFAULT_CHANNEL << "Number of process must be greater or equal to 2 (do you used mpirun?)" << endl;
        MPI::Finalize();
        exit(2);
    }

	Module_DCREATE create;
	create.execute(options);

	return EXIT_SUCCESS;
}

