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

#ifndef DATA_EXCEPTION_H_
#define DATA_EXCEPTION_H_

#include <exception>
#include <sstream>
using namespace std;

namespace errors {

class Data_Exception: public std::exception
{
public:
	virtual ~Data_Exception() throw();
	Data_Exception(const string &comment);
	Data_Exception(const char *comment);
	Data_Exception(long int min, long int max, long int value);
	Data_Exception(long int min, long int max, long int value, const string &comment);
	Data_Exception(long int min, long int max, long int value, const char *comment);
	virtual const char* what() const throw();
	void add_comment(const string &str);
	void add_comment(const char *str);
private:
	long int _min;
	long int _max;
	long int _value;
	string _comment;
	string output;

	void set_output();
};

}

#endif /*DATA_EXCEPTION_H_*/
