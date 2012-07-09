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

#include "Module_DCREATE.h"

#include "io/Auto_Unzip.h"
using namespace useful;

#include "data_structures/Hash.h"

#include <fstream>
#include <sstream>

#define COMMUNICATION_CHANNEL 0
#define DATA_CHANNEL 10

namespace modules {

///////////////////// old code //////////////////

void Module_DCREATE::execute(const Options & options) {
    my_rank = MPI::COMM_WORLD.Get_rank();
    nprocs = MPI::COMM_WORLD.Get_size();
    proc_name = new char[MPI_MAX_PROCESSOR_NAME];
    int resultlen;
    MPI::Get_processor_name(proc_name, resultlen);


    if (my_rank == 0) {
        DEFAULT_CHANNEL << "Master process started on machine " << proc_name << endl;
        compute_master(options);
        DEFAULT_CHANNEL << "Master process waiting for finalize on " << proc_name << endl;
    } else {
        DEFAULT_CHANNEL << "Slave process " << my_rank << '/' << (nprocs-1) << " started on machine " << proc_name << endl;
        receive_from_master();
        DEFAULT_CHANNEL << "Slave process " << my_rank << '/' << (nprocs-1) << " before finalize on " << proc_name << endl;
    }
    MPI::Finalize();
}

void Module_DCREATE::remove_temporary_files() {
	for (vector<string>::iterator iter = temp_files.begin(); iter != temp_files.end(); iter++)
		remove(iter->c_str());
	temp_files.clear();
}

void Module_DCREATE::compute_master(const Options & options) {

	string prefix_temp = options.output_file + string("_temp");

	DEFAULT_CHANNEL << '[' << my_rank << "] reading input" << endl;

	// Read all the Fasta files and check for duplicate names
	vector<Fasta *> multi_fasta;
	set<string> names;
	pair<set<string>::iterator,bool> ret;
	bool all_ok = true;
	size_t sum = 0;
	for (vector<string>::const_iterator iter = options.input_files.begin(); iter != options.input_files.end(); iter++) {
		Auto_Unzip input(iter->c_str());
		while (not input.eof()) {
			Fasta * temp = new Fasta();
			input.filtered() >> *temp;
			sum += temp->length();
			multi_fasta.push_back(temp);
			ret = names.insert(temp->get_id());
			if (ret.second == false) {
				ERROR_CHANNEL << "Error: name \"" << temp->get_id() << "\" already exists!" << endl;
				all_ok = false;
			}
		}
	}
	if (not all_ok) {
		for (int node = 1; node < nprocs; node++)
		send_sequences_to_slave(node, 0, 0, string(), string());
		return;
	}

	DEFAULT_CHANNEL << '[' << my_rank << "] sorting" << endl;
	// sort by length
	if (options.balancing)
		sort(multi_fasta.begin(), multi_fasta.end(), sort_reverse_function);

	DEFAULT_CHANNEL << '[' << my_rank << "] preparing header" << endl;

	// prepare file for header
	stringstream header_name;
	header_name << options.output_file << "_h.dht";
	ofstream o(header_name.str().c_str());
	for (vector<Fasta *>::iterator iter = multi_fasta.begin(); iter != multi_fasta.end(); iter++)
		o << (*iter)->get_id() << '\t' << (*iter)->get_sequence().size() << endl;
	o.close();

	DEFAULT_CHANNEL << '[' << my_rank << "] preparing temporary files" << endl;

	// prepare sets and create temp files
	size_t bins = nprocs;
	size_t bin_length[bins];
	ofstream outputs[bins];
	for (size_t i = 0; i < bins; i++) {
		bin_length[i] = 0;
		stringstream filename;
		filename << prefix_temp << '_' << (i+1) << ".fasta";
		temp_files.push_back(filename.str());
		outputs[i].open(filename.str().c_str());
	}

	DEFAULT_CHANNEL << '[' << my_rank << "] writing to files" << endl;

			// write to files
	for (vector<Fasta *>::iterator iter = multi_fasta.begin(); iter != multi_fasta.end(); iter++) {
		size_t min_pos = 0;
		size_t t_min = sum;
		for (size_t i = 0; i < bins; i++)
			if (bin_length[i] < t_min) {
				min_pos = i;
				t_min = bin_length[min_pos];
			}
		bin_length[min_pos] += (*iter)->length();
		outputs[min_pos] << **iter;
	}

	for (size_t i = 0; i < bins; i++) {
		outputs[i].close();
	}

	DEFAULT_CHANNEL << '[' << my_rank << "] sending to slaves" << endl;

	// send to slaves
	for (size_t i = 1; i < bins; i++) {
		stringstream filename_input;
		stringstream filename_output;
		filename_input << prefix_temp << '_' << (i+1) << ".fasta";
		filename_output << options.output_file << '_' << (i+1) << ".eht";
		send_sequences_to_slave(i, options.k, options.blockLength, filename_input.str(), filename_output.str());
	}

	DEFAULT_CHANNEL << '[' << my_rank << "] clearing memory" << endl;
	// clear memory
	for (vector<Fasta *>::iterator iter = multi_fasta.begin(); iter != multi_fasta.end(); iter++)
		delete *iter;

	DEFAULT_CHANNEL << '[' << my_rank << "] computing" << endl;

	// compute by master process
	stringstream filename_input;
	stringstream filename_output;
	filename_input << prefix_temp << "_1.fasta";
	filename_output << options.output_file << "_1.eht";
	compute_hash(options.k,options.blockLength,filename_input.str().c_str(), filename_output.str().c_str(),false); //TODO handle methyl_hash

	DEFAULT_CHANNEL << '[' << my_rank << "] finishing" << endl;

	stringstream filename_numberfile;
	filename_numberfile << options.output_file << "_n.dht";
	ofstream nf(filename_numberfile.str().c_str());
	if (!nf) {
		ERROR_CHANNEL << "I cannot open file " << filename_numberfile.str() << " for writing!" << endl;
		exit(6);
	}
	nf << nprocs << endl;

}

void Module_DCREATE::send_sequences_to_slave(int node, int k, int blockLength, const string & input_filename, const string & output_filename) const {
	DEFAULT_CHANNEL << "Sending informations to slaves " << node << endl;
	MPI::COMM_WORLD.Send(&k,1,MPI::INT,node,COMMUNICATION_CHANNEL);
	MPI::COMM_WORLD.Send(&blockLength,1,MPI::INT,node,COMMUNICATION_CHANNEL);
	MPI::COMM_WORLD.Send(input_filename.c_str(),input_filename.length()+1,MPI::CHAR,node,COMMUNICATION_CHANNEL);
	MPI::COMM_WORLD.Send(output_filename.c_str(),output_filename.length()+1,MPI::CHAR,node,COMMUNICATION_CHANNEL);
}

void Module_DCREATE::receive_from_master() const {
    size_t length;
    MPI::Status status;
    int k;
    int blockLength;

    MPI::COMM_WORLD.Recv(&k,1,MPI::INT,0,COMMUNICATION_CHANNEL);

    MPI::COMM_WORLD.Recv(&blockLength,1,MPI::INT,0,COMMUNICATION_CHANNEL);

    MPI::COMM_WORLD.Probe(0,COMMUNICATION_CHANNEL,status);
    length = status.Get_count(MPI::CHAR);
    char input[length];
    MPI::COMM_WORLD.Recv(input,length,MPI::CHAR,0,COMMUNICATION_CHANNEL);

    MPI::COMM_WORLD.Probe(0,COMMUNICATION_CHANNEL,status);
    length = status.Get_count(MPI::CHAR);
    char output[length];
    MPI::COMM_WORLD.Recv(output,length,MPI::CHAR,0,COMMUNICATION_CHANNEL);

    DEFAULT_CHANNEL << "Informations from master received by node " << my_rank << endl;

    if (strlen(input) != 0 and strlen(output) != 0)
    	compute_hash(k,blockLength,input,output,false); //TODO handle methyl_hash

}

void Module_DCREATE::compute_hash(int k, int blockLength, const char * input_file, const char * output_file, bool methyl_hash) const {
	clock_t start = clock();
	Hash * H = new Hash(k, blockLength, methyl_hash);
	vector<string> temp;
	temp.push_back(string(input_file));
	H->createHASH(temp);
	clock_t end = clock();
	DEFAULT_CHANNEL << "building time = " <<((end - start) / (double)CLOCKS_PER_SEC) << endl;
	H->save(output_file);
	clock_t end2 = clock();
	DEFAULT_CHANNEL << "saving time = " <<((end2 - end) / (double)CLOCKS_PER_SEC) << endl;
	DEFAULT_CHANNEL << "total time = " <<((end2 - start) / (double)CLOCKS_PER_SEC) << endl;
}




///////////////////// new (broken) code //////////////////


/*
void Module_DCREATE::execute(const Options & options) {
    my_rank = MPI::COMM_WORLD.Get_rank();
    nprocs = MPI::COMM_WORLD.Get_size();
    proc_name = new char[MPI_MAX_PROCESSOR_NAME];
    int resultlen;
    MPI::Get_processor_name(proc_name, resultlen);


    if (my_rank == 0) {
        DEFAULT_CHANNEL << "Master process started on machine " << proc_name << endl;
        compute_master(options);
        DEFAULT_CHANNEL << "Master process waiting for finalize on " << proc_name << endl;
    } else {
        DEFAULT_CHANNEL << "Slave process " << my_rank << '/' << (nprocs-1) << " started on machine " << proc_name << endl;
        receive_from_master();
        DEFAULT_CHANNEL << "Slave process " << my_rank << '/' << (nprocs-1) << " before finalize on " << proc_name << endl;
    }
    MPI::Finalize();
    if (my_rank == 0)
    	DEFAULT_CHANNEL << "Creation finished!" << endl;
}

void Module_DCREATE::remove_temporary_files() {
	for (vector<string>::iterator iter = temp_files.begin(); iter != temp_files.end(); iter++)
		remove(iter->c_str());
	temp_files.clear();
}

void Module_DCREATE::compute_master(const Options & options) {

	string prefix_temp = options.output_file + string("_temp");

	DEFAULT_CHANNEL << '[' << my_rank << "] reading input" << endl;

	// Read all the Fasta files and check for duplicate names
	vector<Fasta *> multi_fasta;
	set<string> names;
	pair<set<string>::iterator,bool> ret;
	bool all_ok = true;
	size_t sum = 0;
	for (vector<string>::const_iterator iter = options.input_files.begin(); iter != options.input_files.end(); iter++) {
		Auto_Unzip input(iter->c_str());
		while (not input.eof()) {
			Fasta * temp = new Fasta();
			input.filtered() >> *temp;
			sum += temp->length();
			multi_fasta.push_back(temp);
			ret = names.insert(temp->get_id());
			if (ret.second == false) {
				ERROR_CHANNEL << "Error: name \"" << temp->get_id() << "\" already exists!" << endl;
				all_ok = false;
			}
		}
	}
	if (not all_ok) {
		for (int node = 1; node < nprocs; node++)
		send_sequences_to_slave(node, 0, 0, string(), string());
		return;
	}

	// sort by length
	if (options.balancing) {
		DEFAULT_CHANNEL << '[' << my_rank << "] sorting" << endl;
		sort(multi_fasta.begin(), multi_fasta.end(), sort_reverse_function);
	} else {
		DEFAULT_CHANNEL << '[' << my_rank << "] skipping sorting" << endl;
	}

	DEFAULT_CHANNEL << '[' << my_rank << "] preparing header" << endl;
	// prepare file for header
	stringstream header_name;
	header_name << options.output_file << "_header.erne";
	ofstream o(header_name.str().c_str());
	for (vector<Fasta *>::iterator iter = multi_fasta.begin(); iter != multi_fasta.end(); iter++)
		o << (*iter)->get_id() << '\t' << (*iter)->get_sequence().size() << endl;
	o.close();

	DEFAULT_CHANNEL << '[' << my_rank << "] preparing temporary files" << endl;
	// prepare sets and create temp files
	size_t bins = nprocs - 1;
	size_t bin_length[bins];
	ofstream outputs[bins];
	for (size_t i = 0; i < bins; i++) {
		bin_length[i] = 0;
		stringstream filename;
		filename << prefix_temp << '_' << (i+1) << ".fasta";
		temp_files.push_back(filename.str());
		outputs[i].open(filename.str().c_str());
	}

	DEFAULT_CHANNEL << '[' << my_rank << "] writing to files" << endl;
	// write to files
	for (vector<Fasta *>::iterator iter = multi_fasta.begin(); iter != multi_fasta.end(); iter++) {
		size_t min_pos = 0;
		size_t t_min = sum;
		for (size_t i = 0; i < bins; i++)
			if (bin_length[i] < t_min) {
				min_pos = i;
				t_min = bin_length[min_pos];
			}
		bin_length[min_pos] += (*iter)->length();
		outputs[min_pos] << **iter;
	}

	for (size_t i = 0; i < bins; i++) {
		outputs[i].close();
	}

	DEFAULT_CHANNEL << '[' << my_rank << "] clearing memory" << endl;
	// clear memory
	for (vector<Fasta *>::iterator iter = multi_fasta.begin(); iter != multi_fasta.end(); iter++)
		delete *iter;

	DEFAULT_CHANNEL << '[' << my_rank << "] sending to slaves" << endl;
	// send to slaves
	for (size_t i = 1; i <= bins; i++) {
		stringstream filename_input;
		stringstream filename_output;
		filename_input << prefix_temp << '_' << (i) << ".fasta";
		filename_output << options.output_file << '_' << (i) << ".erne";
		int k;
		if (options.auto_k)
			k = Hash::calculate_k(filename_input.str(),options.methyl_hash);
		else
			k = options.k;
		send_sequences_to_slave(i, k, options.blockLength, filename_input.str(), filename_output.str());
	}

	DEFAULT_CHANNEL << '[' << my_rank << "] computing" << endl;
	// compute by master process
	stringstream filename_input;
	stringstream filename_output;
	filename_input << prefix_temp << "_1.fasta";
	filename_output << options.output_file << "_1.erne";

	int k;
	if (options.auto_k)
		k = Hash::calculate_k(filename_input.str(),options.methyl_hash);
	else
		k = options.k;

	compute_hash(k,options.blockLength,filename_input.str().c_str(), filename_output.str().c_str(), options.methyl_hash);

	DEFAULT_CHANNEL << '[' << my_rank << "] finishing" << endl;
	stringstream filename_numberfile;
	filename_numberfile << options.output_file << "_n.inf";
	ofstream nf(filename_numberfile.str().c_str());
	if (!nf) {
		ERROR_CHANNEL << "I cannot open file " << filename_numberfile.str() << " for writing!" << endl;
		exit(6);
	}
	nf << nprocs << endl;

}

void Module_DCREATE::send_sequences_to_slave(int node, int k, int blockLength, const string & input_filename, const string & output_filename) const {
	DEFAULT_CHANNEL << "Sending informations to slave " << node << endl;
	MPI::COMM_WORLD.Send(&k,1,MPI::INT,node,COMMUNICATION_CHANNEL);
	MPI::COMM_WORLD.Send(&blockLength,1,MPI::INT,node,COMMUNICATION_CHANNEL);
	MPI::COMM_WORLD.Send(input_filename.c_str(),input_filename.length()+1,MPI::CHAR,node,COMMUNICATION_CHANNEL);
	MPI::COMM_WORLD.Send(output_filename.c_str(),output_filename.length()+1,MPI::CHAR,node,COMMUNICATION_CHANNEL);
}

void Module_DCREATE::receive_from_master() const {
    size_t length;
    MPI::Status status;
    int k;
    int blockLength;

    MPI::COMM_WORLD.Recv(&k,1,MPI::INT,0,COMMUNICATION_CHANNEL);

    MPI::COMM_WORLD.Recv(&blockLength,1,MPI::INT,0,COMMUNICATION_CHANNEL);

    MPI::COMM_WORLD.Probe(0,COMMUNICATION_CHANNEL,status);
    length = status.Get_count(MPI::CHAR);
    char input[length];
    MPI::COMM_WORLD.Recv(input,length,MPI::CHAR,0,COMMUNICATION_CHANNEL);

    MPI::COMM_WORLD.Probe(0,COMMUNICATION_CHANNEL,status);
    length = status.Get_count(MPI::CHAR);
    char output[length];
    MPI::COMM_WORLD.Recv(output,length,MPI::CHAR,0,COMMUNICATION_CHANNEL);

    DEFAULT_CHANNEL << "Informations from master received by node " << my_rank << endl;

    if (strlen(input) != 0 and strlen(output) != 0)
    	compute_hash(k,blockLength,input,output,false); // TODO: controllare methylation

}

void Module_DCREATE::compute_hash(int k, int blockLength, const char * input_file, const char * output_file, bool methyl_hash) const {
	clock_t start = clock();
	Hash * H = new Hash(k, blockLength, methyl_hash);
	vector<string> temp;
	temp.push_back(string(input_file));
	H->createHASH(temp);
	clock_t end = clock();
	DEFAULT_CHANNEL << "building time = " <<((end - start) / (double)CLOCKS_PER_SEC) << endl;
	H->save(output_file);
	clock_t end2 = clock();
	DEFAULT_CHANNEL << "saving time = " <<((end2 - end) / (double)CLOCKS_PER_SEC) << endl;
	DEFAULT_CHANNEL << "total time = " <<((end2 - start) / (double)CLOCKS_PER_SEC) << endl;
}
*/

}
