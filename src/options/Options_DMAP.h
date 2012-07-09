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

#ifndef OPTIONS_DMAP_H_
#define OPTIONS_DMAP_H_

#include "options/Options.h"

namespace options {

class Options_DMAP : public Options {
public:
	Options_DMAP() : Options() { balancing = true; }
	Options_DMAP(int argc, char *argv[]) : Options() {
		balancing = true;
		if (not process(argc,argv))
			exit(2);
	}
	virtual ~Options_DMAP() { }

	bool process(int argc, char *argv[]);

protected:
};

}

#endif /* OPTIONS_DMAP_H_ */
