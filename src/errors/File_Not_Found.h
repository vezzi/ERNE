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

#ifndef FILE_NOT_FOUND_H_
#define FILE_NOT_FOUND_H_

#include <exception>
#include <sstream>
using namespace std;

namespace errors
{

class File_Not_Found: public std::exception
{
public:
	virtual ~File_Not_Found() throw();
	File_Not_Found(const char* comment);
	File_Not_Found(const string &comment);
	File_Not_Found(const string* comment);
	virtual const char* what() const throw();

private:
	string _comment;
};

}

#endif /*FILE_NOT_FOUND_H_*/
