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

#include "Temporary_File.h"

#include <cstdlib>
#include <vector>

namespace samobjects {

Temporary_File::Temporary_File(const char * path) {
	if (path == NULL)
		filename = open_temp(string("/tmp"),file);
	else
		filename = open_temp(string(path),file);
	open = true;
}

Temporary_File::~Temporary_File() {
	remove(filename.c_str());
}

void Temporary_File::close_file() {
	if (open) {
		file.close();
		open = false;
	}
}


std::string Temporary_File::open_temp(std::string path, std::ofstream& f) {
    path += "/XXXXXX";
    std::vector<char> dst_path(path.begin(), path.end());
    dst_path.push_back('\0');

    int fd = mkstemp(&dst_path[0]);
    if(fd != -1) {
        path.assign(dst_path.begin(), dst_path.end() - 1);
        f.open(path.c_str(),
               std::ios_base::trunc | std::ios_base::out);
        close(fd);
    }
    return path;
}

}
