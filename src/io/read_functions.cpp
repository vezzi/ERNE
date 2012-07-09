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

#include "io/read_functions.h"

mutex read_mutex;

#define GUI_OUTPUT_MAX_VALUE 10000
int output_prev_value = -1;

void output_progress(Auto_Unzip & input, bool gui_output) {
	if (gui_output) {
		int now = (1.0 * input.tellg() / input.fixed_length()) * GUI_OUTPUT_MAX_VALUE;
		if (now != output_prev_value) {
			output_prev_value = now;
			ERROR_CHANNEL << ">P\t" << output_prev_value << '\n';
		}
	} else {
		int now = (1.0 * input.tellg() / input.fixed_length())*1000;
		if (now != output_prev_value) {
			output_prev_value = now;
			VERBOSE_CHANNEL << (output_prev_value / 10.0) << "% done." << endl;
		}
	}
}

int read_sequences(Auto_Unzip & input, int num_seq, Mask sequences[], Fasta::FASTQ_encoding format_type, bool gui_output) {
	if (&input == NULL)
		return 0;

	Fasta read;
	read.set_FASTQ_type(format_type);
	int n_seq = 0;
	{
		mutex::scoped_lock lock(read_mutex);
		istream & in = input.filtered();
		while (not input.eof() and n_seq < num_seq) {
			in >> read;
			if (read.length() > MAX_READ_LENGTH) {
				cerr << "Read " << read.get_id() << " too long. Max allowed read size is " << MAX_READ_LENGTH << endl;
				exit(5);
			}
			//Reset and take a ref
			Mask & r = sequences[n_seq] = Mask();
			r.set_id(read.get_id());
			r.set_sequence(read.get_sequence());
			r.set_quality(read.get_quality());
			n_seq++;
		}
		output_progress(input, gui_output);
	}
	return n_seq;
}

int read_sequences(Auto_Unzip & first, Auto_Unzip & second, int num_seq, Mask sequences[], Fasta::FASTQ_encoding format_type, bool gui_output) {
	if (&first == NULL or &second == NULL)
		return 0;

	Fasta read;
	read.set_FASTQ_type(format_type);
	int n_seq = 0;
	{
		mutex::scoped_lock lock(read_mutex);
		istream & first_in = first.filtered();
		istream & second_in = second.filtered();
		while (not first.eof() and not second.eof() and (n_seq + 1) < num_seq) {
			first_in >> read;
			if (read.length() > MAX_READ_LENGTH) {
				cerr << "Read " << read.get_id() << " too long. Max allowed read size is " << MAX_READ_LENGTH << endl;
				exit(5);
			}
			Mask & rf = sequences[n_seq];
			rf.set_id(read.get_id());
			rf.set_sequence(read.get_sequence());
			rf.set_quality(read.get_quality());
			n_seq++;
			second_in >> read;
			if (read.length() > MAX_READ_LENGTH) {
				cerr << "Read " << read.get_id() << " too long. Max allowed read size is " << MAX_READ_LENGTH << endl;
				exit(5);
			}
			Mask & rs = sequences[n_seq];
			rs.set_id(read.get_id());
			rs.set_sequence(read.get_sequence());
			rs.set_quality(read.get_quality());
			n_seq++;
			//CHECK!!
			if (rf.id.compare(0, rf.id.size() - 1, rs.id, 0, rs.id.size() - 1) != 0) {
				ERROR_CHANNEL << "wrong paired reads IDs: '" << rf.id << "' and '" << rs.id << '\'' << endl;
				exit(2);
			}
		}
		output_progress(first, gui_output);
	}
	return n_seq;
}

