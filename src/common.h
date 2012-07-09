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

#ifndef TYPES_H_
#define TYPES_H_

#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <climits>
#include <cstdlib>

#ifdef INLINE_DISABLED
#define INLINE
#else
#define INLINE inline
#endif

#define DEFAULT_CHANNEL std::cout
#define VERBOSE_CHANNEL std::cerr
#define ERROR_CHANNEL std::cerr

static inline std::string package_description() {
	std::string line(PACKAGE_NAME);
	line.append(" version ");
	line.append(PACKAGE_VERSION);
	return line;
}

#define text_delimitator '$'
#define masked_base 'N'

//#define LONG_LENGTH

#ifdef LONG_LENGTH
	typedef long int t_length; ///< Type for TEXT
	#define MAX_T_LENGTH LONG_MAX ///< Max value for TEXT
	#define SA_NOT_FOUND -1
#else
	typedef int t_length; ///< Type for TEXT
	#define MAX_T_LENGTH INT_MAX ///< Max value for TEXT
	#define SA_NOT_FOUND -1
#endif

typedef unsigned short int t_char; ///< Type for char conversion
typedef unsigned short int t_errors; ///< Type for ERRORS
typedef t_errors t_errors_delta ;
typedef unsigned int t_pattern_length; ///< Type for PATTERN
typedef double t_quality; ///< Type for QUALITY VALUES

typedef std::vector<short unsigned int> t_quality_vector; ///< Vector with quality values

typedef std::string t_edit_string; ///< For future expansion

enum t_strand { forward_strand, reverse_strand, unknown_strand }; ///< Just an enumeration of possible strands

//enum t_alignment { unique_alignment, multiple_alignments, alignments_not_found, quality_discarded, low_complexity, contamination, unknown_alignment };

enum t_masked { masked, not_masked };

/**
 * Conversion from strand type to char type
 */
static inline char strand_to_char(const t_strand s) {
	if (s == unknown_strand)
		return '.';
	else
		return (s == forward_strand) ? '+' : '-';
}

/**
 * Conversion from char type to strand type
 */
static inline t_strand char_to_strand(const char c) {
	if (c == '+' or c == 'F')
		return forward_strand;
	else if (c == '-' or c == 'R')
		return reverse_strand;
	else
		return unknown_strand;
}

static inline char reverse_complement_standalone_char(const char c) {
	switch (c) {
	case 'A' : return 'T';
	case 'T' : return 'A';
	case 'C' : return 'G';
	case 'G' : return 'C';
	case 'U' : return 'A';
	case 'R' : return 'Y';
	case 'Y' : return 'R';
	case 'M' : return 'K';
	case 'K' : return 'M';
	case 'W' : return 'S';
	case 'S' : return 'W';
	case 'B' : return 'V';
	case 'V' : return 'B';
	case 'D' : return 'H';
	case 'H' : return 'D';
	//case 'N' : return 'N';
	//case 'X' : return 'X';
	default : return c;
	}
}

static inline std::string reverse_complement_standalone_str_length(const char *str, size_t length) {
	std::string reverse;
	size_t i = length;
	while (i > 0)
		reverse.push_back(reverse_complement_standalone_char(str[--i]));
	return reverse;
}

static inline std::string reverse_complement_standalone_str(const char *str) {
	return reverse_complement_standalone_str_length(str,strlen(str));
}

static inline std::string reverse_complement_standalone_str(const std::string & str) {
	return reverse_complement_standalone_str_length(str.c_str(),str.length());
}


static inline std::string reverse_standalone_str_length(const char *str, size_t length) {
	std::string reverse;
	size_t i = length;
	while (i > 0)
		reverse.push_back(str[--i]);
	return reverse;
}

static inline std::string reverse_standalone_str(const char *str) {
	return reverse_standalone_str_length(str,strlen(str));
}

static inline std::string tolower(const std::string & old_string) {
	std::string lower_string(old_string);
	for (std::string::iterator iter = lower_string.begin(); iter != lower_string.end(); iter++) {
		*iter = tolower(*iter);
	}
	return lower_string;
}

static inline std::string toupper(const std::string & old_string){
	std::string upper_string(old_string);
	for (std::string::iterator iter = upper_string.begin(); iter != upper_string.end(); iter++) {
		*iter = toupper(*iter);
	}
	return upper_string;
}

static inline std::string bam_integer(int value) {
	// da ottimizzare
	std::string r;
	if (-128 <= value and value <= 127) {
		r.push_back('c');
		r.push_back((char)value);
	} else if (-128*256 <= value and value <= 128*256-1) {
		r.push_back('s');
		r.push_back((char)(value%256));
		r.push_back((char)(value/256));
	} else {
		r.push_back('i');
		r.push_back((char)(value % 256));
		value = value / 256;
		r.push_back((char)(value % 256));
		value = value / 256;
		r.push_back((char)(value % 256));
		value = value / 256;
		r.push_back((char)(value));
	}
	return r;
}

static inline std::string bam_uinteger(unsigned int value) {
	// da ottimizzare
	std::string r;
	if (value <= 255) {
		r.push_back('C');
		r.push_back((char)value);
	} else if (value <= 256*256-1) {
		r.push_back('S');
		r.push_back((char)(value % 256));
		r.push_back((char)(value/256));
	} else {
		r.push_back('I');
		r.push_back((char)(value % 256));
		value = value / 256;
		r.push_back((char)(value % 256));
		value = value / 256;
		r.push_back((char)(value % 256));
		value = value / 256;
		r.push_back((char)(value));
	}
	return r;
}

#endif /*TYPES_H_*/
