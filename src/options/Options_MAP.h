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

#ifndef OPTIONS_MAP_H_
#define OPTIONS_MAP_H_

#include "options/Options.h"

namespace options {

class Options_MAP : public Options {
public:
	Options_MAP() : Options() { }
	Options_MAP(const Options & orig) : Options(orig) { }
	Options_MAP(int argc, char *argv[]) : Options() {
		if (not process(argc,argv))
			exit(2);
	}
	virtual ~Options_MAP() { }

	bool process(int argc, char *argv[]);
};

}

#endif /* OPTIONS_MAP_H_ */
