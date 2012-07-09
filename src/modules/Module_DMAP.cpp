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

#include "Module_DMAP.h"

#include <mpi.h>

#define COMMUNICATION_CHANNEL 0
#define DATA_CHANNEL 10
#define NUM_OF_SEQ 11
#define LOCK_FREE 12

namespace modules {

///////////// OLD VERSION ///////////////


void Module_DMAP::execute(const Options & options) {
	clock_t started = clock();
	processed = 0;

	initialize_parameters(options);

	my_rank = MPI::COMM_WORLD.Get_rank();
	nprocs = MPI::COMM_WORLD.Get_size();
	proc_name = new char[MPI_MAX_PROCESSOR_NAME];
	int resultlen;
	MPI::Get_processor_name(proc_name, resultlen);

	// Check for the correct number of files
	if (my_rank == 0) {
		stringstream filename_numberfile;
		filename_numberfile << options.reference_file << "_n.dht";
		ifstream nf(filename_numberfile.str().c_str());
		if (nf.fail()) {
			ERROR_CHANNEL << "File " << filename_numberfile.str() << " not found! Check the reference name!" << endl;
			exit(5);
		}
		string s;
		getline(nf,s);
		int n_files = -1;
		n_files = atoi(s.c_str());
		DEFAULT_CHANNEL << "Number of rNA files in the set: " << s << endl;
		if (n_files != nprocs) {
			ERROR_CHANNEL << "Wrong number of files: expected " << nprocs << " files but the data structure was built for " << n_files << " processes!" << endl;
			exit(5);
		} else {
			DEFAULT_CHANNEL << '[' << my_rank << "] Found " << n_files << '/' << nprocs << " files" << endl;
		}
	}

	if (my_rank == 0)
		DEFAULT_CHANNEL << '[' << my_rank << "] Reading process started on machine " << proc_name << endl;
	else if (my_rank == nprocs-1)
		DEFAULT_CHANNEL << '[' << my_rank << "] Writing process started on machine " << proc_name << endl;
	else
		DEFAULT_CHANNEL << '[' << my_rank << "] Slave process " << my_rank << "/" << (nprocs-1) << " started on machine " << proc_name << endl;

	execute_generic_worker(options);

	if (my_rank == 0)
		DEFAULT_CHANNEL << '[' << my_rank << "] Reading process waiting for finalize on " << proc_name << endl;
	else if (my_rank == nprocs-1)
		DEFAULT_CHANNEL << '[' << my_rank << "] Writing process waiting for finalize on " << proc_name << endl;
	else
		DEFAULT_CHANNEL << '[' << my_rank << "] Slave process " << my_rank << "/" << (nprocs-1) << " waiting for finalize on " << proc_name << endl;

	MPI::Finalize();

	if (my_rank == 0) {
		clock_t finished = clock();
		double dif = (finished - started) / (double)CLOCKS_PER_SEC;

		DEFAULT_CHANNEL << "Time spent: " << dif << "s (";
		print_formatted_time(DEFAULT_CHANNEL, dif);
		DEFAULT_CHANNEL << ')' << endl;
		//DEFAULT_CHANNEL << "Queries per second: " << (processed/dif)	<<endl;
	}

}


void Module_DMAP::execute_generic_worker(const Options & options) {
	/* TODO: handle contamination check
	if (options.contamination_check) {
		if (options.verbose)
			DEFAULT_CHANNEL << "Reading contamination SA in file " << options.contamination_file << endl;
		clock_t start = clock();
		CR.load(options.contamination_file.c_str());
		clock_t end = clock();

		DEFAULT_CHANNEL << "Hash (contamination) loading time = ";
		print_formatted_time(DEFAULT_CHANNEL, (end - start) / (double)CLOCKS_PER_SEC);
		DEFAULT_CHANNEL << endl;
	}
	 */
	{
		DEFAULT_CHANNEL << '[' << my_rank << "] Loading hash: " << options.reference_file << '_' << (my_rank + 1) << ".eht" << endl;

		clock_t start = clock();
		stringstream filename;
		filename << options.reference_file << '_' << (my_rank + 1) << ".eht";
		H.load(filename.str().c_str());
		clock_t end = clock();

		DEFAULT_CHANNEL << '[' << my_rank << "] Hash (reference) loading time = ";
		print_formatted_time(DEFAULT_CHANNEL, (end - start) / (double)CLOCKS_PER_SEC);
		DEFAULT_CHANNEL << endl;
	}

	if (my_rank == 0) {
		if (not force_fastqformat) {
			fastqformat = Fasta::check_FASTQ_type_file(query1);
			if (fastqformat == Fasta::unknown_fastq_encoding) {
				ERROR_CHANNEL << "The program can not autodetect input format file, please use --force-illumina or --force-standard options" << endl;
				exit(4);
			}
			DEFAULT_CHANNEL << "[0] FASTQ format: " << (fastqformat == Fasta::illumina_fastq_encoding ? "ILLUMINA" : "STANDARD") << " (autodetected)" << endl;
		}
		// open file for input
		if (options.paired_ends) {
			DEFAULT_CHANNEL << '[' << my_rank << "] Opening " << query1 << " and " << query2 << " files" << endl;
			input_file1.open(query1);
			input_file2.open(query1);
		} else {
			DEFAULT_CHANNEL << '[' << my_rank << "] Opening " << query1 << " file" << endl;
			input_file1.open(query1);
		}
	} else if (my_rank == (nprocs -1)) {
		// open file for output
		DEFAULT_CHANNEL << '[' << my_rank << "] Opening output " << options.output_file << " file" << endl;
		output_samfile.open_file_mpi(options);
	} // else do nothig: generic worker

	// retrieve information for contig conversion
	{
		stringstream header_name;
		header_name << options.reference_file << "_h.dht";

		contig_conversion.link(H.globalToLocal,header_name.str());
	}

	thread_group threads;
	if (options.paired_ends) {
		// PROCESS PAIR ENDS FILES
		if (my_rank == 0)
			DEFAULT_CHANNEL << "[0] Processing files " << options.query1 << " and " << options.query2 << endl;
		for (int id=0; id < threads_number; id++)
			threads.create_thread(boost::bind(&this->generic_worker_paired_thr, this, id));
		threads.join_all();
	} else {
		// PROCESS SINGLE READ FILES
		if (my_rank == 0)
			DEFAULT_CHANNEL << "[0] Processing file " << options.query1 << endl;
		for (int id=0; id < threads_number; id++)
			threads.create_thread(boost::bind(&this->generic_worker_single_thr, this, id));
		threads.join_all();
		// Send finalization to next node
		if (my_rank != (nprocs -1)) {
			int number = 0;
			MPI::COMM_WORLD.Send(&number,1,MPI::INT,my_rank+1,COMMUNICATION_CHANNEL);
		}
	}

	DEFAULT_CHANNEL << '[' << my_rank << "] Finished!" << endl;
}

/* static */
Module_DMAP::Random_Choice_Result Module_DMAP::random_choice_from_previous(int previous, int actual) {
	int random_value = (int)(guessRes() * (previous + actual));
	if (random_value < previous)
		return Random_Choice_Result(true,0);
	else
		return Random_Choice_Result(false, random_value - previous);
}

/* static */
void Module_DMAP::generic_worker_single_thr(Module_DMAP * search, int id) {
	Mask * sequences;
	int read_seq;
	vector<Mask> printable_solutions;
	Transmitting_Result received;

	if (search->my_rank == 0) {
		sequences = new Mask[SEQUENCES_FOR_BLOCK];
		read_seq = read_sequences(search->input_file1, SEQUENCES_FOR_BLOCK, sequences, search->fastqformat, search->gui_output);
	} else {
		received = search->receive_from_previous(id);
		sequences = received.first;
		read_seq = received.second;
	}

	while (read_seq != 0) {
		Items solutions;
		//ItemsGapped solutions_gapped;
		t_errors calculated_errors;

		for (int i = 0; i < read_seq; i++) {
			Mask & s = sequences[i];
			if (search->my_rank == 0 and search->trim) {
					//	check_routine(sequences[i], 0);
					s.quality_trimming_MOTT(search->min_phred_value_MOTT,search->min_mean_quality,search->min_size);
			}
			if (s.status_discarded()) {
				if (s.status_low_complexity())
					s.low_complexity  =true; //s.set_type(low_complexity);
				else
					s.low_quality = true; //s.set_type(quality_discarded);
				printable_solutions.push_back(s);
				continue;
			}
			if (search->auto_errors)
				calculated_errors = round((double)s.get_good_length() / search->errors_rate);
			else
				calculated_errors = search->common_errors_allowed;
			if (search->my_rank != 0 and s.algn > 0 and s.NM < calculated_errors)
				calculated_errors = s.NM;
			t_errors count = 0;
			for (t_pattern_length i = s.get_good_region_start()-1; (i < s.get_good_region_stop()) and (count <= calculated_errors); i++)
				if (s.sequence[i] == 'N' or s.sequence[i] == 'n')
					count++;
			if (count > calculated_errors) {
				//s.set_type(alignments_not_found);
				printable_solutions.push_back(s);
				continue;
			}

			/** ALIGNMENT **/
			solutions.clear();

			if (search->my_rank == 0 and search->contamination_check) {
				search->CR.search(s.get_good_sequence(),solutions,calculated_errors);
				if (solutions.size() > 0)
					s.contaminated = true;
			}

			if (not s.contaminated)
				search->H.search(s.get_good_sequence(),solutions,calculated_errors);
			if (solutions.size() == 0) {
				/** Try gapped **/
				/*
				solutions_gapped.clear();
				if (search->gap)
					search->H.search_gapped(s.get_good_sequence(),solutions_gapped,search->seed_sizes,search->seed_errors,calculated_errors,search->max_gap);
				*/
				/*
				if (solutions_gapped.size() == 0) {// size 0 means no alignment found
				*/
					//s.set_type(alignments_not_found);
					printable_solutions.push_back(s);
					continue;
				/*
				} else {
					if (not search->printAll) {
						Random_Choice_Result r;
						bool improved = (s.NM + s.NM_gap) > (solutions_gapped.at(0).errors1 + solutions_gapped.at(0).errors2);

						if (improved)
							r = Search_MPI::random_choice_from_previous(0,solutions_gapped.size());
						else
							r = Search_MPI::random_choice_from_previous(s.algn,solutions_gapped.size());
						if (not improved and r.first) {
							// take the previous solution
							s.algn += solutions_gapped.size();
						} else {
							// update solution
							const ResultItemGapped & HM = solutions_gapped.at(r.second);
							s.globalPosition = HM.GlobalPosition1;
							if (improved)
								s.algn = solutions_gapped.size();
							else
								s.algn += solutions_gapped.size();
							s.HI = 1;
							s.IH = 1;
							s.primary = true;
							s.strand = HM.strand;
							s.NM = HM.errors1;
							s.NM_gap = HM.errors2;
							s.contig = HM.contig;
							s.position = HM.GlobalPosition1 - search->H.globaltolocal.startPositions[HM.contig] + 1 ;
							s.position_gap = HM.GlobalPosition2 - search->H.globaltolocal.startPositions[HM.contig] + 1 ;
							s.length1_gap = HM.length1;
							s.length2_gap = HM.length2;
							s.contig = search->contig_conversion.convert(s.contig);

						}
						printable_solutions.push_back(s);


					} else { // printALL
				*/
						/*
						unsigned int processed=0;
						unsigned int alignments;
						(search->toBePrinted < solutions_gapped.size()) ? alignments = search->toBePrinted : alignments = solutions_gapped.size() ;
						while(processed < alignments) {
							const ResultItemGapped & HM = solutions_gapped.at(processed);
							s.globalPosition = HM.GlobalPosition1;
							s.algn = alignments;
							s.HI =1;
							s.IH =1;
							(processed == 0 ) ? s.primary = true : s.primary = false;
							s.strand = HM.strand;
							s.NM = HM.errors1;
							s.NM_gap = HM.errors2;
							s.contig = HM.contig ;
							s.position = HM.GlobalPosition1 - search->H.globaltolocal.startPositions[HM.contig] + 1 ;
							s.position_gap = HM.GlobalPosition2 - search->H.globaltolocal.startPositions[HM.contig] + 1 ;
							s.length1_gap = HM.length1;
							s.length2_gap = HM.length2;
							s.contaminated = contaminated;
							printable_solutions.push_back(s);
							processed++;
						}
						 */
/*
						ERROR_CHANNEL << "--print-all option not implemented yet!" << endl;
						exit(3);
					}
				}
				*/
			} else if (not search->printAll) {
				sort(solutions.begin(), solutions.end(), ResultItem::less()); // sort solutions
				solutions.erase(unique(solutions.begin(), solutions.end(), ResultItem::equal()), solutions.end());

				Random_Choice_Result r;
				bool improved = (s.NM + s.NM_gap) > (solutions.at(0).errors);

				if (improved)
					r = Module_DMAP::random_choice_from_previous(0,solutions.size());
				else {
					r = Module_DMAP::random_choice_from_previous(s.algn,solutions.size());
					s.algn += solutions.size();
				}
				if (not r.first) {
					const ResultItem & HM = solutions.at(r.second);
					s.HI = 1;
					s.IH = 1;
					s.primary = true;
					s.globalPosition = HM.globalPosition;
					s.strand = HM.strand;
					s.NM = HM.errors;
					s.NM_gap = 0;
					if ((search->my_rank == 0) and s.contaminated) {
						s.contig = search->CR.globalToLocal.searchContig(HM.globalPosition); // find the contig/scaffold
						s.position = HM.globalPosition - search->CR.globalToLocal.startPositions[s.contig] + 1;
						s.contig = search->contig_conversion.convert(s.contig);
					} else {
						s.contig = search->H.globalToLocal.searchContig(HM.globalPosition); // find the contig/scaffold
						s.position = HM.globalPosition - search->H.globalToLocal.startPositions[s.contig] + 1;
						s.contig = search->contig_conversion.convert(s.contig);
					}
				}
				printable_solutions.push_back(s);
				continue;
			} else { // printAll
				/*
				// memorize all printable solutions
				sort(solutions.begin(), solutions.end(), ResultItem::less()); // sort solutions
				solutions.erase(unique(solutions.begin(), solutions.end(), ResultItem::equal()), solutions.end());
				unsigned int processed=0;
				unsigned int alignments;
				(search->toBePrinted < solutions.size()) ? alignments = search->toBePrinted : alignments = solutions.size() ;
				while(processed < alignments) {
					// while I print enough solution or there are no more solutions
					const ResultItem & HM = solutions.at(processed);
					s.globalPosition = HM.GlobalPosition;
					s.algn = solutions.size();
					s.IH = alignments;
					s.HI = processed+1;
					(processed == 0 ) ? s.primary = true : s.primary = false;
					s.strand = HM.strand;
					s.NM = HM.errors;
					if (contaminated) {
						s.contig = search->CR.globaltolocal.searchContig(HM.GlobalPosition); // find the contig/scaffold
						s.position = HM.GlobalPosition - search->CR.globaltolocal.startPositions[s.contig] + 1;
					} else {
						s.contig = search->H.globaltolocal.searchContig(HM.GlobalPosition); // find the contig/scaffold
						s.position = HM.GlobalPosition - search->H.globaltolocal.startPositions[s.contig] + 1;
					}
					s.contaminated = contaminated;
					printable_solutions.push_back(s);
				}
				 */

				ERROR_CHANNEL << "--print-all option not implemented yet!" << endl;
				exit(3);
			}
		}

		if (search->my_rank == (search->nprocs-1)) {
			// now print all
			for(unsigned int i=0; i < printable_solutions.size(); i++)
				search->output_samfile.print_output(printable_solutions.at(i));
			search->processed += read_seq;
		} else // send data to next node
			search->send_to_next(printable_solutions,id);

		delete [] sequences;
		printable_solutions.clear();

		if (search->my_rank == 0) {
			sequences = new Mask[SEQUENCES_FOR_BLOCK];
			read_seq = read_sequences(search->input_file1, SEQUENCES_FOR_BLOCK, sequences, search->fastqformat, search->gui_output);
			if (read_seq == 0)
				delete [] sequences;
		} else {
			received = search->receive_from_previous(id);
			sequences = received.first;
			read_seq = received.second;
		}
	}

}

/* static */
void Module_DMAP::generic_worker_paired_thr(Module_DMAP * search, int id) {

}

void Module_DMAP::send_to_next(const vector<Mask> &printable_solutions,int id) {
	if (printable_solutions.size() > 0) {
		send_output(printable_solutions,my_rank+1,id);
	}
}

Module_DMAP::Transmitting_Result Module_DMAP::receive_from_previous(int id) {
	Transmitting_Result result;
	result = recv_output(my_rank-1,id);
	return result;
}


void Module_DMAP::send_output(const vector<Mask> &printable_solutions, int node, int id) {
	int count = printable_solutions.size();
	if (count == 0)
		return;
	//DEFAULT_CHANNEL << '[' << my_rank << ',' << id << "] Send " << count << " OUTPUTs from node " << my_rank << " to node " << node << endl;
	size_t sum = 0;

	for (vector<Mask>::const_iterator iter = printable_solutions.begin(); iter != printable_solutions.end(); iter++)
		sum += iter->id.size() + 2*iter->sequence.size() + 3;

	char * informations = new char[sum];
	char * h = informations;
	size_t temp;
	for (vector<Mask>::const_iterator iter = printable_solutions.begin(); iter != printable_solutions.end(); iter++) {
		memcpy(h, iter->get_id().c_str(), iter->get_id().size() + 1);
		h += iter->id.size() + 1;
		temp = iter->sequence.size();
		memcpy(h, iter->sequence.c_str(), temp + 1);
		h += temp + 1;
		memcpy(h, iter->quality.c_str(), temp + 1);
		h += temp + 1;
	}

	unsigned long int * positions = new unsigned long int[count*2];
	unsigned long int * global_positions = new unsigned long int[count*2];
	int * contigs = new int[count];
	int * NMs = new int[count*2];
	int * lengths = new int[count*2];
	int * algn = new int[count];
	//t_alignment * types = new t_alignment[count];
	unsigned short int * bools = new unsigned short int[count];
	unsigned int * trim_info = new unsigned int[count*2];
	unsigned short int bool_temp;
	for (int i = 0; i < count; i++) {
		const Mask & r = printable_solutions.at(i);
		positions[i*2] = r.position;
		positions[i*2+1] = r.position_gap;
		global_positions[i*2] = r.globalPosition;
		global_positions[i*2+1] = r.globalPosition_gap;
		lengths[i*2] = r.length1_gap;
		lengths[i*2+1] = r.length2_gap;
		contigs[i] = r.contig;
		//types[i] = r.type;
		NMs[i*2] = r.NM;
		NMs[i*2+1] = r.NM_gap;
		algn[i] = r.algn;
		trim_info[i*2] = r.good_region_start;
		trim_info[i*2+1] = r.good_region_stop;
		bool_temp = 0;
		if (r.strand) bool_temp |= 0x01;
		if (r.masked) bool_temp |= 0x02;
		if (r.low_quality) bool_temp |= 0x04;
		if (r.trimmed) bool_temp |= 0x08;
		if (r.discarded) bool_temp |= 0x10;
		if (r.low_complexity) bool_temp |= 0x20;
		if (r.contaminated) bool_temp |= 0x40;
		if (r.gapped) bool_temp |= 0x80;
		bools[i] = bool_temp;
	}

	// sending !!!
	{
		mutex::scoped_lock lock(mpi_mutex);
		MPI::COMM_WORLD.Send(&count,1,MPI::INT,my_rank+1,COMMUNICATION_CHANNEL);
		MPI::COMM_WORLD.Send(informations,sum,MPI::CHAR,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Send(positions,count*2,MPI::UNSIGNED_LONG,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Send(global_positions,count*2,MPI::UNSIGNED_LONG,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Send(contigs,count,MPI::INT,node,DATA_CHANNEL);
		//MPI::COMM_WORLD.Send(types,count * sizeof(t_alignment),MPI::CHAR,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Send(NMs,count*2,MPI::INT,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Send(lengths,count*2,MPI::INT,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Send(algn,count,MPI::INT,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Send(bools,count,MPI::UNSIGNED_SHORT,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Send(trim_info,count*2,MPI::UNSIGNED,node,DATA_CHANNEL);
	}
	delete [] positions;
	delete [] contigs;
	delete [] NMs;
	delete [] algn;
	delete [] bools;
	delete [] trim_info;
	//delete [] types;
	delete [] lengths;
	delete [] global_positions;
	delete [] informations;

	//DEFAULT_CHANNEL << '[' << my_rank << ',' << id << "] Sent " << count << " OUTPUTs from node " << my_rank << " to node " << node << endl;
}

Module_DMAP::Transmitting_Result Module_DMAP::recv_output(int node, int id) {
	int count;
	unsigned long int * positions;
	unsigned long int * global_positions;
	int * contigs;
	//t_alignment * types;
	int * NMs;
	int * lengths;
	int * algn;
	unsigned short int * bools;
	unsigned int * trim_info;
	char * informations;
	unsigned short int bool_temp;
	Mask * reads;

	{
		mutex::scoped_lock lock(mpi_mutex);
		//DEFAULT_CHANNEL << '[' << my_rank << ',' << id << "] Waiting info from node " << node << " to node " << my_rank << endl;
		if (finished)
			return Transmitting_Result(NULL,0);
		MPI::COMM_WORLD.Recv(&count,1,MPI::INT,my_rank-1,COMMUNICATION_CHANNEL);
		//DEFAULT_CHANNEL << '[' << my_rank << ',' << id << "] Receive " << count << " OUTPUTs from node " << node << " to node " << my_rank << endl;
		if (count == 0) {
			finished = true;
			return Transmitting_Result(NULL,0);
		}

		positions = new unsigned long int[count*2];
		global_positions = new unsigned long int[count*2];
		contigs = new int[count];
		//types = new t_alignment[count];
		NMs = new int[count*2];
		lengths = new int[count*2];
		algn = new int[count];
		bools = new unsigned short int[count];
		trim_info = new unsigned int[count*2];

		size_t sum;
		MPI::Status status;
		MPI::COMM_WORLD.Probe(node,DATA_CHANNEL,status);
		sum = status.Get_count(MPI::CHAR);
		informations = new char[sum];
		MPI::COMM_WORLD.Recv(informations,sum,MPI::CHAR,node,DATA_CHANNEL);

		MPI::COMM_WORLD.Recv(positions,count*2,MPI::UNSIGNED_LONG,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Recv(global_positions,count*2,MPI::UNSIGNED_LONG,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Recv(contigs,count,MPI::INT,node,DATA_CHANNEL);
		//MPI::COMM_WORLD.Recv(types,count*sizeof(t_alignment),MPI::CHAR,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Recv(NMs,count*2,MPI::INT,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Recv(lengths,count*2,MPI::INT,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Recv(algn,count,MPI::INT,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Recv(bools,count,MPI::UNSIGNED_SHORT,node,DATA_CHANNEL);
		MPI::COMM_WORLD.Recv(trim_info,count*2,MPI::UNSIGNED,node,DATA_CHANNEL);
	}
	reads = new Mask[count];
	char * h = informations;
	for (int i = 0; i < count; i++) {
		Mask & r = reads[i];
		r.id = string(h);
		h += r.id.size() + 1;
		r.sequence = string(h);
		h += r.sequence.size() + 1;
		r.quality = string(h);
		h += r.sequence.size() + 1;

		r.position = positions[i*2];
		r.position_gap = positions[i*2+1];
		r.globalPosition = global_positions[i*2];
		r.globalPosition_gap = global_positions[i*2+1];
		r.contig = contigs[i];
		r.length1_gap = lengths[i*2];
		r.length2_gap = lengths[i*2+1];
		//r.type = types[i];
		r.NM = NMs[i*2];
		r.NM_gap = NMs[i*2+1];
		r.algn = algn[i];
		r.good_region_start = trim_info[i*2];
		r.good_region_stop  = trim_info[i*2+1];

		bool_temp = bools[i];
		r.strand = bool_temp & 0x01;
		r.masked = bool_temp & 0x02;
		r.low_quality = bool_temp & 0x04;
		r.trimmed = bool_temp & 0x08;
		r.discarded = bool_temp & 0x10;
		r.low_complexity = bool_temp & 0x20;
		r.contaminated = bool_temp & 0x40;
		r.gapped = bool_temp & 0x80;
	}
	delete [] positions;
	delete [] contigs;
	//delete [] types;
	delete [] NMs;
	delete [] algn;
	delete [] bools;
	delete [] trim_info;
	delete [] lengths;
	delete [] global_positions;
	delete [] informations;

	//DEFAULT_CHANNEL << '[' << my_rank << ',' << id << "] Received " << count << " OUTPUTs from node " << node << " to node " << my_rank << endl;
	return Transmitting_Result(reads,count);
}



///////////// NEW (BROKEN) VERSION ///////////////


/*
void Module_DMAP::execute(const Options & options) {
	clock_t started = clock();
	processed = 0;

	initialize_parameters(options);

	my_rank = MPI::COMM_WORLD.Get_rank();
	nprocs = MPI::COMM_WORLD.Get_size();
	nworkers = nprocs - 1;
	proc_name = new char[MPI_MAX_PROCESSOR_NAME];
	int resultlen;
	MPI::Get_processor_name(proc_name, resultlen);


	if (my_rank == 0) {
		DEFAULT_CHANNEL << "Master process started on machine " << proc_name << endl;
		execute_master(options);
	} else {
		DEFAULT_CHANNEL << "Worker process " << my_rank << '/' << nworkers << " started on machine " << proc_name << endl;
		execute_worker(options);
	}

	MPI::Finalize();

	if (my_rank == 0) {
		DEFAULT_CHANNEL << "Master process waiting for finalize on " << proc_name << endl;

		clock_t finished = clock();
		double dif = (finished - started) / (double) CLOCKS_PER_SEC;

		DEFAULT_CHANNEL << "Time spent: " << dif << "s (";
		print_formatted_time(DEFAULT_CHANNEL, dif);
		DEFAULT_CHANNEL << ')' << endl;
		DEFAULT_CHANNEL << "Queries per second: " << (processed / dif) << endl;
	}

}
*/

/*
void Module_DMAP::execute_master(const Options & options) {

	//Load input FASTQ file(s)
	//Check format
	if (not force_fastqformat) {
		fastqformat = Fasta::check_FASTQ_type_file(query1);
		if (fastqformat == Fasta::unknown_fastq_encoding) {
			ERROR_CHANNEL << "The program can not autodetect input format file, please use --force-illumina or --force-standard options" << endl;
			exit(4);
		}
		DEFAULT_CHANNEL << "FASTQ format: " << (fastqformat == Fasta::illumina_fastq_encoding ? "ILLUMINA" : "STANDARD") << " (autodetected)" << endl;
	}
	//Open file
	if (options.paired_ends) {
		DEFAULT_CHANNEL << "Opening " << query1 << " and " << query2 << " files...";
		input_file1.open(query1);
		input_file2.open(query1);
		DEFAULT_CHANNEL << "\t[OK]" << endl;
	} else {
		DEFAULT_CHANNEL << "Opening " << query1 << " file...";
		input_file1.open(query1);
		DEFAULT_CHANNEL << "\t[OK]" << endl;
	}

	//Open output SAM/BAM file
	DEFAULT_CHANNEL << "Opening output " << options.output_file << " file...";
	output_samfile.open_file_mpi(options);
	DEFAULT_CHANNEL << "\t[OK]" << endl;

	//Read sequenced
	Mask sequences[SEQUENCES_FOR_BLOCK];

	//Allocate share sequences data and lock structure
	void* mem = MPI::Alloc_mem((sizeof (Mask::MaskData) * SEQUENCES_FOR_BLOCK) + (sizeof (bool) * nprocs), MPI::INFO_NULL);
	Mask::MaskData* sharedSequencesData = (Mask::MaskData*) mem;
	bool* lockStruct = (bool*) ((char *)mem + sizeof (Mask::MaskData) * SEQUENCES_FOR_BLOCK);

	//Initialize lock structure
	for (int i = 0; i < nprocs; i++) {

		//Unset lock
		lockStruct[i] = false;
	}

	//Enable shared memory
	DEFAULT_CHANNEL << "Initializing shared memory and synchronizing Workers...";
	win = MPI::Win::Create(sharedSequencesData, (sizeof (Mask::MaskData) * SEQUENCES_FOR_BLOCK) + (sizeof (bool) * nprocs), 1, MPI::INFO_NULL, MPI::COMM_WORLD);
	DEFAULT_CHANNEL << "\t[OK]" << endl;

	//Auxiliary data
	int readSequences;
	int goodSequences;
	t_errors calculated_errors;

	while ((readSequences = read_sequences(input_file1, SEQUENCES_FOR_BLOCK, sequences, fastqformat, gui_output)) != 0) {

		goodSequences = 0;

		//For every sequence...
		for (int i = 0; i < readSequences; i++) {

			//Take a reference
			Mask &s = sequences[i];

			//Trim sequence
			if (trim) {
				s.quality_trimming_MOTT(min_phred_value_MOTT, min_mean_quality, min_size);
			}

			//Check status
			if (s.status_discarded()) {
				output_samfile.print_output(s);
				continue;
			}

			//Determinate error
			if (auto_errors)
				calculated_errors = round((double) s.get_good_length() / errors_rate);
			else
				calculated_errors = common_errors_allowed;

			//Check indeterminate nucleotide
			t_errors count = 0;
			for (t_pattern_length i = s.get_good_region_start() - 1; (i < s.get_good_region_stop()) and (count <= calculated_errors); i++)
				if (s.sequence[i] == 'N' or s.sequence[i] == 'n')
					count++;
			if (count > calculated_errors) {
				output_samfile.print_output(s);
				continue;
			}

			//Add sequence to share data
			s.toData(sharedSequencesData[goodSequences++]);
		}

		//If there are some sequences to align...
		if (goodSequences > 0) {

			//(Re)Start alignment on all Workers
			for (int rank = 1; rank < nprocs; rank++) {
				//Send number of sequences
				MPI::COMM_WORLD.Send(&goodSequences, 1, MPI::INT, rank, NUM_OF_SEQ);
			}

			//Global Synchronization
			MPI::COMM_WORLD.Barrier();
		}

		//Print solution
		Mask s;
		for (int i = 0; i < goodSequences; i++) {

			//Save solution to output file
			s.updateFromData(sharedSequencesData[i]);
			output_samfile.print_output(s);
		}

		//Update processed read counter
		processed += readSequences;
	}

	//Close SAM/BAM output file
	DEFAULT_CHANNEL << "Closing output " << options.output_file << " file...";
	output_samfile.close_file();
	DEFAULT_CHANNEL << "\t[OK]" << endl;

	//Stop all Workers
	goodSequences = 0;
	for (int rank = 1; rank < nprocs; rank++) {

		DEFAULT_CHANNEL << "Stopping  Worker " << rank << '/' << nworkers << " ...";

		//Send stop message
		MPI::COMM_WORLD.Send(&goodSequences, 1, MPI::INT, rank, NUM_OF_SEQ);

		DEFAULT_CHANNEL << "\t[OK]" << endl;

	}

	//Close shared memory
	DEFAULT_CHANNEL << "Removing shared memory and shutting down Workers...";
	win.Free();
	DEFAULT_CHANNEL << "\t[OK]" << endl;

	return;
}
*/

/*
void Module_DMAP::execute_worker(const Options & options) {

	//Loading reference file
	stringstream filename;
	filename << options.reference_file << '_' << my_rank << ".eht";
	H.load(filename.str().c_str());
	H.set_indels_max_value(indels_max_value);


	//Retrieve information for contig conversion
	stringstream header_name;
	header_name << options.reference_file << "_header.eht";
	contig_conversion.link(H.globaltolocal, header_name.str());

	/* TODO: handle contamination check
        if (options.contamination_check) {
            if (options.verbose)
                    DEFAULT_CHANNEL << "Reading contamination SA in file " << options.contamination_file << endl;
            clock_t start = clock();
            CR.load(options.contamination_file.c_str());
            clock_t end = clock();

            DEFAULT_CHANNEL << "Hash (contamination) loading time = ";
            print_formatted_time(DEFAULT_CHANNEL, (end - start) / (double)CLOCKS_PER_SEC);
            DEFAULT_CHANNEL << endl;
        }
	 */
/*
	//Create thread group
	thread_group threads;

	//Enable shared memory
	win = MPI::Win::Create(MPI::BOTTOM, 0, 1, MPI::INFO_NULL, MPI::COMM_WORLD);

	//Share data
	int numberOfSequences;
	Mask::MaskData localSequencesData[SEQUENCES_FOR_BLOCK];

	//Receive start message
	MPI::COMM_WORLD.Recv(&numberOfSequences, 1, MPI::INT, 0, NUM_OF_SEQ);


	while (numberOfSequences != 0) {

		//Lock memory in shared mode
		win.Lock(MPI::LOCK_SHARED, 0, 0);

		//Read all share data
		win.Get(&localSequencesData, sizeof (Mask::MaskData) * numberOfSequences, MPI::BYTE, 0, 0, sizeof (Mask::MaskData) * numberOfSequences, MPI::BYTE);

		//Unlock shared memory
		win.Unlock(0);

		//Sequence elaboration
		//Paired ends reads
		if (options.paired_ends) {

			//Create and launch thread
			for (int id = 0; id < threads_number; id++)
				threads.create_thread(boost::bind(&worker_paired_thr, this, id, localSequencesData, numberOfSequences, win));

			//Wait all threads
			threads.join_all();

		} else {//Single reads

			//Create and launch thread
			for (int id = 0; id < threads_number; id++)

				threads.create_thread(boost::bind(&worker_single_thr, this, id, localSequencesData, numberOfSequences, win));


			//Wait all threads
			threads.join_all();
		}

		//Global Synchronization
		MPI::COMM_WORLD.Barrier();

		//Receive start message
		MPI::COMM_WORLD.Recv(&numberOfSequences, 1, MPI::INT, 0, NUM_OF_SEQ);

	}

	//Close shared memory
	win.Free();

	return;
}
*/




/* static */
/*
Module_DMAP::Random_Choice_Result Module_DMAP::random_choice_from_previous(int previous, int actual) {
	int random_value = (int) (guessRes() * (previous + actual));
	if (random_value < previous)
		return Random_Choice_Result(true, 0);
	else
		return Random_Choice_Result(false, random_value - previous);
}
*/

/* static */
/*
void Module_DMAP::worker_single_thr(Module_DMAP * search, int id, Mask::MaskData* sequences, int numberOfSequences, MPI::Win &win) {

	//Auxiliary data
	Items solutions;
	t_errors calculated_errors;
	bool improved;
	Mask::MaskData updatedSeq;

	//For every sequence associated to thread
	for (int currentSequence = id; currentSequence < numberOfSequences; currentSequence += search->threads_number) {


		//Retrieve sequence
		Mask s;
		s.updateFromData(sequences[currentSequence]);


		//Determinate error
		if (search->auto_errors)
			calculated_errors = round((double) s.get_good_length() / search->errors_rate);
		else
			calculated_errors = search->common_errors_allowed;

		//Alignment
		solutions.clear();
		improved = false;

		//Check if sequence is contaminated
		if (search->contamination_check) {
			search->CR.search(s.get_good_sequence(), solutions, calculated_errors);
			if (solutions.size() > 0)
				s.contaminated = true;
		}

		//Alignment
		if (not s.contaminated)
			search->H.search_MPI(s.get_good_sequence(), solutions, calculated_errors, currentSequence, win, search->thread_mutex);

		//If no alignment found...
		/*
		if (solutions.size() == 0) {

			//Only one thread by time can read/write to Master
			mutex::scoped_lock lock(search->thread_mutex);

			//Mutual exclusion in update operations
			search->UpdateRMA_mutex_lock(win);

			//Lock memory
			win.Lock(MPI::LOCK_SHARED, 0, 0);

			//Get update
			win.Get(&updatedSeq, sizeof (Mask::MaskData), MPI::BYTE, 0, sizeof (Mask::MaskData) * currentSequence, sizeof (Mask::MaskData), MPI::BYTE);

			//Unlock shared memory
			win.Unlock(0);

			if (updatedSeq.algn == 0) {

				updatedSeq.type = alignments_not_found;

				//Lock memory
				win.Lock(MPI::LOCK_EXCLUSIVE, 0, 0);

				//Update Master
				win.Put(&updatedSeq, sizeof (Mask::MaskData), MPI::BYTE, 0, currentSequence * sizeof (Mask::MaskData), sizeof (Mask::MaskData), MPI::BYTE);

				//Unlock shared memory
				win.Unlock(0);

			}

			search->UpdateRMA_mutex_unlock(win);

		} else */ /*
		if (not search->printAll) { //Alignment(s) found ands save only one random alignment position

			//Only one thread by time can read/write to Master
			mutex::scoped_lock lock(search->thread_mutex);

			//Mutual exclusion in update operations
			search->UpdateRMA_mutex_lock(win);

			//Lock memory
			win.Lock(MPI::LOCK_SHARED, 0, 0);

			//Get update
			win.Get(&updatedSeq, sizeof (Mask::MaskData), MPI::BYTE, 0, sizeof (Mask::MaskData) * currentSequence, sizeof (Mask::MaskData), MPI::BYTE);

			//Unlock shared memory
			win.Unlock(0);

			//Sort solutions
			sort(solutions.begin(), solutions.end(), ResultItem::less());
			//Remove duplicates
			solutions.erase(unique(solutions.begin(), solutions.end(), ResultItem::equal()), solutions.end());

			//If no better solution found...
			if ((updatedSeq.algn > 0) and (solutions.size() > 0) and (solutions.at(0).errors > (updatedSeq.NM + updatedSeq.NM_gap))) {

				search->UpdateRMA_mutex_unlock(win);

				//Skip sequence
				continue;
			}

			Random_Choice_Result r;

			//Determinate if was found a best occurrence of sequence
			if (solutions.size() > 0)
				improved = (solutions.at(0).errors) < (updatedSeq.NM + updatedSeq.NM_gap);
			else
				improved = false;

			//If best occurrence found...
			if (improved) {
				//Choose a random alignment among all improved alignments
				r = Module_DMAP::random_choice_from_previous(0, solutions.size());
				//Update number of alignments
				updatedSeq.algn = solutions.size();
			} else {
				//Keep previous alignment
				r = Module_DMAP::random_choice_from_previous(updatedSeq.algn, solutions.size());
				//Add found alignments to previous alignments
				updatedSeq.algn += solutions.size();
			}

			if (solutions.size() > 0 and not r.first) {//New align position
				//Extract a alignment

				const ResultItem & HM = solutions.at(r.second);

				//Save alignment position
				updatedSeq.globalPosition = HM.globalPosition;

				//Update other data
				updatedSeq.primary = true;
				updatedSeq.IH = 1;
				updatedSeq.HI = 1;
				updatedSeq.strand = HM.strand;
				updatedSeq.NM = HM.errors;
				updatedSeq.NM_gap = 0;
				if (updatedSeq.contaminated) {
					updatedSeq.contig = search->contig_conversion.convert(search->CR.globaltolocal.searchContig(HM.globalPosition));
					updatedSeq.position = HM.globalPosition - search->CR.globaltolocal.startPositions[updatedSeq.contig] + 1;
					updatedSeq.contig = search->contig_conversion.convert(updatedSeq.contig);
				} else {
					updatedSeq.contig = search->H.globaltolocal.searchContig(HM.globalPosition);
					updatedSeq.position = HM.globalPosition - search->H.globaltolocal.startPositions[updatedSeq.contig] + 1;
					updatedSeq.contig = search->contig_conversion.convert(updatedSeq.contig);
				}
			}

			//Lock memory
			win.Lock(MPI::LOCK_EXCLUSIVE, 0, 0);

			//Update Master
			win.Put(&updatedSeq, sizeof (Mask::MaskData), MPI::BYTE, 0, sizeof (Mask::MaskData) * currentSequence, sizeof (Mask::MaskData), MPI::BYTE);

			//Unlock shared memory
			win.Unlock(0);

			search->UpdateRMA_mutex_unlock(win);

		} else { //Save all alignments positions

			ERROR_CHANNEL << "--print-all option not implemented yet!" << endl;
			exit(3);
		}

	}
}
*/

/* static */
/*
void Module_DMAP::worker_paired_thr(Module_DMAP * search, int id, Mask::MaskData* sequences, int numberOfSequences, MPI::Win & win) {

}
*/

/*
void Module_DMAP::UpdateRMA_mutex_lock(MPI::Win &win) {

	//Auxiliary data
	bool activeLock = true;
	bool shareStruct[nprocs];

	//Lock shared memory
	win.Lock(MPI::LOCK_EXCLUSIVE, 0, 0);

	//Add request lock
	win.Put(&activeLock, 1, MPI::BOOL, 0, (sizeof (Mask::MaskData) * SEQUENCES_FOR_BLOCK) + (sizeof (bool) * my_rank), 1, MPI::BOOL);


	//Get lock status
	win.Get(&shareStruct, nprocs, MPI::BOOL, 0, sizeof (Mask::MaskData) * SEQUENCES_FOR_BLOCK, nprocs, MPI::BOOL);

	//Unlock shared memory
	win.Unlock(0);

	bool lockSet = false;

	//For all locks...
	for (int i = 1; i < nprocs; i++) {

		//Except mine...
		if (i == my_rank)
			continue;

		//Check if it is locked
		if (shareStruct[i]) {
			lockSet = true;
			break;
		}
	}

	if (not lockSet) {
		return;
	} else {

		int dummyBuffer;

		//Wait granted lock message
		MPI::COMM_WORLD.Recv(&dummyBuffer, 0, MPI::INT, MPI::ANY_SOURCE, LOCK_FREE);

		return;
	}
}

void Module_DMAP::UpdateRMA_mutex_unlock(MPI::Win &win) {

	//Auxiliary data
	bool unactiveLock = false;
	bool shareStruct[nprocs];

	//Lock shared memory
	win.Lock(MPI::LOCK_EXCLUSIVE, 0, 0);

	//Add request lock
	win.Put(&unactiveLock, 1, MPI::BOOL, 0, (sizeof (Mask::MaskData) * SEQUENCES_FOR_BLOCK) + (sizeof (bool) * my_rank), 1, MPI::BOOL);

	//Get lock structure
	win.Get(shareStruct, nprocs, MPI::BOOL, 0, sizeof (Mask::MaskData) * SEQUENCES_FOR_BLOCK, nprocs, MPI::BOOL);

	//Unlock shared memory
	win.Unlock(0);

	int waitingLock = -1;

	//For all locks...
	for (int i = 1; i < nprocs; i++) {

		//Except mine...
		if (i == my_rank)
			continue;

		//Search a waiting lock
		if (shareStruct[i] == true) {
			waitingLock = i;
			break;
		}
	}

	//If there is a waiting lock
	if (waitingLock != -1) {

		int dummyBuffer;

		//Send granted lock message
		MPI::COMM_WORLD.Isend(&dummyBuffer, 0, MPI::INT, waitingLock, LOCK_FREE);
	}

	return;
}
*/

}

