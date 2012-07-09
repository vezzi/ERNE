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
#ifndef TEMPORARY_FILE_H_
#define TEMPORARY_FILE_H_

#include <fstream>
using namespace std;

namespace samobjects {

class Temporary_File {
public:
	Temporary_File(const char * path = NULL);
	virtual ~Temporary_File();
	ofstream & get_stream() { return file; }
	string get_filename() const { return filename; }
	void close_file();
protected:
	ofstream file;
	string filename;
	bool open;

	string open_temp(string path, ofstream& f);

};

}

#endif /* TEMPORARY_FILE_H_ */
