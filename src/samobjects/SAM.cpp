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

#include "SAM.h"

#include "Temporary_File.h"

namespace samobjects {

SAM::SAM() {
	unset();
}

SAM::~SAM() throw (SAM_IO_Error) {
	clear();
}

void SAM::unset() {
	sam_file = NULL;
	filename.clear();
	samfile_is_open = false;
}

void SAM::clear() throw (SAM_IO_Error){
	close();
	filename.clear();
}

void SAM::close() {
	if (sam_file != NULL) {
		samfile_is_open = false;
		samclose(sam_file);
		sam_file = NULL;
	}
}

void SAM::revert() throw (SAM_IO_Error) {
	close();
	open();
}

SAM::sam_format_t SAM::detect_file_type_by_name(const char * filename) {
	const regex bam_re("^.*\\.bam$");
	if (regex_match(filename, bam_re))
		return bam_format;

	const regex tam_re("^.*\\.[st]am$");
	if (regex_match(filename, tam_re))
		return tam_format;

	if (strlen(filename) == 1 and filename[0] == '-') {
		return tam_format;
	}

	const regex tam_compressed_re("^.*\\.[st]am\\.(gz|bz2)$");
	if (regex_match(filename, tam_compressed_re))
		return tam_compressed_format;

	return unknown_sam_format;
}

void SAM::popolate_name_list(bam_header_t * header) {
	if (header == NULL)
		return;
	names_list.clear();
	names_to_tid.clear();

	for (int32_t i = 0; i < header->n_targets; i++) {
		names_list.push_back(string(header->target_name[i]));
		names_to_tid[string(header->target_name[i])] = i;
	}
}

int32_t SAM::name_to_tid(const string * name) const {
	name_to_tid_t::const_iterator iter = names_to_tid.find(*name);
	if (iter == names_to_tid.end())
		return -1;
	else
		return iter->second;
}

/*static*/ bam_header_t * SAM::update_header_from_list(bam_header_t *header, names_list_t & list) {
	Temporary_File samfile;
	samfile.close_file();
	samfile_t * sf = samopen(samfile.get_filename().c_str(),"wh",header);
	samclose(sf);

	Temporary_File tempfile;
	ofstream &output = tempfile.get_stream();

	ifstream input(samfile.get_filename().c_str());
	string temp;
	while (not input.eof()) {
		getline(input,temp);
		if ((temp.size() >= 3) and (temp[0] != '@' or temp[1] != 'S' or temp[2] != 'Q'))
			output << temp << '\n';
	}

	for (names_list_t::iterator iter = list.begin(); iter != list.end(); iter++)
		output << "@SQ\tSN:" << iter->first << "\tLN:" << iter->second << '\n';
	tempfile.close_file();

	tamFile fp = sam_open(tempfile.get_filename().c_str());

	bam_header_t * newheader = sam_header_read(fp);
	sam_close(fp);

	return newheader;
}

/*static*/ void SAM::regex_header(bam_header_t * header, const char * pattern, const char * replace) {
	string old_s(header->text, header->l_text);

	regex re(pattern);
	if (regex_search(old_s, re)) {
		string rep(replace);
		string new_s = regex_replace (old_s, re, rep);

		free(header->text);
		header->l_text = new_s.size();
		header->text = (char *)malloc(header->l_text);
		memcpy(header->text,new_s.c_str(),header->l_text);
	}
}

bool SAM::check_ordered() const {
	if (is_open()) {
		string text (sam_file->header->text, sam_file->header->l_text);
		regex re("SO:coordinate");
		if (regex_search(text, re)) {
			return true;
		}
	}
	return false;
}

void SAM::change_to_unsorted() {
	regex_header(sam_file->header,"SO:coordinate","SO:unsorted");
}

const string * SAM::tid_to_name(int32_t i) const {
	if (i < 0 or (size_t)i >= names_list.size())
		return NULL;
	else
		return &(names_list.at(i));
}

}
