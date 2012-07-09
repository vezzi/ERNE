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

#ifndef AUTO_UNZIP_H_
#define AUTO_UNZIP_H_

#include <fstream>
#include <iostream>
using namespace std;

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
using namespace boost;

#include "errors/File_Not_Found.h"
using namespace errors;

namespace useful {

class Auto_Unzip {
public:
	Auto_Unzip();
	Auto_Unzip(const char * filename) throw (File_Not_Found);
	Auto_Unzip(const string & filename) throw (File_Not_Found);
	virtual ~Auto_Unzip();

	void open(const char * filename) throw (File_Not_Found);
	void open(const string & filename) throw (File_Not_Found) ;
	istream & filtered() {
		if (filtered_stream == NULL)
			return *input_stream;
		else
			return *filtered_stream;
	}
	ifstream & file() { return *input_stream; }
	streampos tellg();
	streampos length();
	streampos fixed_length() { return dimension; }
	void seekg(streampos pos);

	bool eof() const { return (filtered_stream == NULL) ? input_stream->eof() : filtered_stream->eof(); }

protected:
	ifstream * input_stream;
	iostreams::filtering_istream * filtered_stream;

	streampos dimension;

};

}

#endif /* AUTO_UNZIP_H_ */
