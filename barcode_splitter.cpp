#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib> 
#include <cstring>
#include <set>

#include "BKTree.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/regex.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

namespace po = boost::program_options;

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

class bc_splitter {

	public:
	//bc_splitter();
	bool parse_args(int argc, char* argv[]);
	int distance(std::string source, std::string target);	
	unsigned long updateMaps(std::string& barcode_str,
    	std::string& lword1, std::string& lword2,
    	std::string& lword3, std::string& lword4,
    	std::string& rword1, std::string& rword2,
    	std::string& rword3, std::string rword4,
    	unsigned long totalcap);	
	void writeMapsToFile();
	void split_engine();
	void write_log();
	BKTree<std::string>& getTree();
	void initialize();
	void print_help();

	private:
	int cutoff;
	std::string ltype;
	std::string dict_file;
	std::string file1_str;
	std::string file2_str;
	std::string outdirpath;
	int barcode_start;
	int barcode_size;
	int umi_start;
	int umi_size;
	int allowed_MB;
	std::map<std::string, std::vector<std::string>> lQueueMap;
    std::map<std::string, std::vector<std::string>> rQueueMap;
	std::set<std::string> barcode_set;
	std::set<std::string> outfile_set;
	std::map<std::string, unsigned long> zero_dist_map;
	std::map<std::string, unsigned long> one_dist_map;
	std::map<std::string, unsigned long> higher_dist_map;
    BKTree<std::string> tree;
	po::options_description desc;
	std::map<int, int> distmap;

	unsigned long totalcap = 0;
	unsigned long  match_total = 0;
	unsigned long ambiguous_total = 0;
	unsigned long no_match_total = 0;

	bool validUmi = false;

};

void bc_splitter::print_help() {
    std::cout << desc << "\n";
	std::cout << "Usage: bc_splitter -t allseq -d <dict_file> --file1 <file1> --file2 <file2> -o outdir\n\n";
}


BKTree<std::string>& bc_splitter::getTree() {
	return tree;
}


void bc_splitter::initialize() {
	std::ifstream iff(dict_file);
    boost::archive::text_iarchive iar(iff);
        
    iar >> tree;

	struct stat st = {0};

	if (stat(outdirpath.c_str(), &st) == -1) {
		mkdir(outdirpath.c_str(), 0755);
	}
}

bool bc_splitter::parse_args(int argc, char* argv[]) {
	bool all_set = true;
	desc.add_options()
		("help,h", "produce help message")
		("type,t", po::value(&ltype), "allseq/rnatagseq")
		("dict_file,d", po::value<std::string>(&dict_file), "Dictionary file")
		("bc-start", po::value(&barcode_start)->default_value(6), 
			"Optional/Barcode start position")
		("bc-size", po::value(&barcode_size)->default_value(6), 
			"Optional/Barcode size")
		("file1", po::value<std::string>(&file1_str), "First file")
		("file2", po::value<std::string>(&file2_str), "Second file")
		("outdir,o", po::value<std::string>(&outdirpath), "Output directory")
		("mismatch,m", po::value(&cutoff)->default_value(1), 
			"Optional/Maximum allowed mismatches.")
		("umi-start", po::value(&umi_start)->default_value(0), 
			"Optional/Umi start position")
		("umi-size", po::value(&umi_size)->default_value(6), 
			"Optional/Umi size")
		("allowed-mb", po::value(&allowed_MB)->default_value(2048),
			"Optional/Estimated memory requirement in MB.")
	;

	po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

	
	if (vm.count("help")) {
		print_help();
        return 0;
    }

	if (vm.count("type")) {
		boost::to_lower(ltype);
		boost::trim(ltype);
		// check the ltype
		std::cout << "Type is set to: " << ltype << "\n";
		if (ltype.compare("allseq") == 0) {
			// Let the default values for barcode and umi
			validUmi = true;
		} else if ((ltype.compare("rnatagseq") == 0) ||
			(ltype.compare("rnatag-seq") == 0)) {
			// May not be nextseek; Have to confirm with Jonathan.
			barcode_start = 0;
			barcode_size = 8;
			umi_start = 0;
			umi_size = 0;
		} else if (ltype.compare("customized") == 0) {
			std::cout << "Note! The user need to set the variables for" <<
			"barcode_start, barcode_size, umi_start, umi_size\n";
		} else {
			std::cout << "Error: Invalid type option.\n";
			all_set = false;
		}
	} else {
		std::cout << "Error: Type not set.";
		all_set = false;
	}

	std::cout << "Max mismatch is set to " << cutoff << ".\n";
	std::cout << "Umi-start is set to " << umi_start << ".\n";
	std::cout << "Umi-size is set to " << umi_size << ".\n";
	std::cout << "Barcode-start is set to " << barcode_start << ".\n";
	std::cout << "Barcode-size is set to " << barcode_size << ".\n";

	if (vm.count("file1")) {
		std::cout << "First fastq file is set to: " << file1_str << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: First fastq file is not set.\n";
	}

	if (vm.count("file2")) {
		std::cout << "Second fastq file is set to: " << file2_str << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: Second fastq file is not set.\n";
	}

	if (vm.count("dict_file")) {
		std::cout << "Dict_file is set to " << dict_file << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: Dict_file is not set.\n";
	}

	if (vm.count("outdir")) {
		std::cout << "Outdir is set to: " << outdirpath << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: Outdir is not set.\n";
	}
	
	return all_set;
}

int bc_splitter::distance(std::string source, std::string target) {

    const int n = source.length();
    const int m = target.length();
    if (n == 0) {
        throw std::invalid_argument("Length of source is zero.");
    }
    if (m == 0) {
        throw std::invalid_argument("Length of target is zero.");
    }

    if (m !=n ) {
        throw std::invalid_argument("Source and target have different length");
	}

    int ldist = 0;

    for (int j = 0; j < n; j++) {
        if (source[j] != target[j]) {
            ldist++;
        }
    }

    return ldist;
}	


unsigned long 
bc_splitter::updateMaps(std::string& barcode_str, 
	std::string& lword1, std::string& lword2, 
	std::string& lword3, std::string& lword4, 
	std::string& rword1, std::string& rword2, 
	std::string& rword3, std::string rword4, 
	unsigned long totalcap) {

	// Get the capacity of the the eight strings as we have to insert
	// them into the two maps.
	int lcap1 = lword1.capacity();
	int lcap2 = lword2.capacity();
	int lcap3 = lword3.capacity();
	int lcap4 = lword4.capacity();

	int rcap1 = rword1.capacity();
	int rcap2 = rword2.capacity();
	int rcap3 = rword3.capacity();
	int rcap4 = rword4.capacity();
	
	int allcap = lcap1 + lcap2 + lcap3 + lcap4 + rcap1 + rcap2 + rcap3 + rcap4;

	totalcap += allcap;

	lQueueMap[barcode_str].push_back(lword1);
	lQueueMap[barcode_str].push_back(lword2);
	lQueueMap[barcode_str].push_back(lword3);
	lQueueMap[barcode_str].push_back(lword4);

	rQueueMap[barcode_str].push_back(rword1);
	rQueueMap[barcode_str].push_back(rword2);
	rQueueMap[barcode_str].push_back(rword3);
	rQueueMap[barcode_str].push_back(rword4);

	return totalcap;
}

void bc_splitter::writeMapsToFile() {
	
	for (auto& kv : lQueueMap) {
	    std::string barcode = kv.first;

        const std::string file1 = outdirpath + "/" + barcode + "_1.fastq";
        const std::string file2 = outdirpath + "/" + barcode + "_2.fastq";

        std::ofstream ofs1;
        std::ofstream ofs2;
        if (outfile_set.count(barcode) > 0) {
	        ofs1 = std::ofstream(file1, std::ofstream::out|std::ofstream::app);
            ofs2 = std::ofstream(file2, std::ofstream::out|std::ofstream::app);
        } else {

            ofs1 = std::ofstream(file1, std::ofstream::out|std::ofstream::trunc);
            ofs2 = std::ofstream(file2, std::ofstream::out|std::ofstream::trunc);
            outfile_set.insert(barcode);
        }

        // Dump the content of the two maps to the two files

        std::vector<std::string> & valSet1 = lQueueMap[barcode];
        std::vector<std::string> & valSet2 = rQueueMap[barcode];

        for (auto const& kv1 : valSet1) {
 	        std::string val1 = kv1;
            ofs1 << val1 << '\n';
        }

        ofs1.close();

		for (auto const& kv2 : valSet2) {
        	std::string val2 = kv2;
            ofs2 << val2 << '\n';
        }

        ofs2.close();
	}

	lQueueMap.clear();
	rQueueMap.clear();

} 

void bc_splitter::split_engine() {

	// The words starting with lword is for read 1
	std::string lword1;
	std::string lword2;
	std::string lword3;
	std::string lword4;

	// The words starting with rword is for read 2
	std::string rword1;
	std::string rword2;
	std::string rword3;
	std::string rword4;

	for (int j = 0; j <= cutoff; j++) {
		distmap[j] = 0;
	}   
	const unsigned long MB_SIZE = 1024 * 1024;
	unsigned long total_allowed = allowed_MB * MB_SIZE;

	//const std::string logfile_detailed = outdirpath + "/logfile_detailed.txt";
	//std::ofstream log_detailed(logfile_detailed);

	// We shall start reading the first line. The assumption is that second line 
	// contains the barcode.
	std::ifstream file1(file1_str);
	std::ifstream file2(file2_str);

	while (std::getline(file1, lword1)) {

		if (!std::getline(file1, lword2)) {break;}
		if (!std::getline(file1, lword3)) {break;}
		if (!std::getline(file1, lword4)) {break;}
		
		if (!std::getline(file2, rword1)) {break;}
		if (!std::getline(file2, rword2)) {break;}
		if (!std::getline(file2, rword3)) {break;}
		if (!std::getline(file2, rword4)) {break;}

		// The barcode stays at the second line of each four lines of first
		// read file.
		std::string barcode_str = lword2.substr(barcode_start, barcode_size);	
		std::string umi_str;
		if (validUmi) {
			umi_str = lword2.substr(umi_start, umi_size);
		}
		std::vector<std::string> results;
		
		results = tree.find(barcode_str, cutoff);

		// calculate the minimum distance between the target and references

		int smallest_dist = cutoff + 1;
		std::string smallest_barcode = "JJJJJJ";

		std::vector<int> dist_vec;
		for (auto const& val : results) {
			int ldist = distance(val, barcode_str);
			if (ldist < smallest_dist) {
				smallest_dist = ldist;
				smallest_barcode = val;
			}
			dist_vec.push_back(ldist);
		}

		int smallest_count = 0;
		for (auto const& temp_dist : dist_vec) {
			if (temp_dist == smallest_dist) {
				smallest_count++;
			}
		}

		//log_detailed << "actual_barcode: " << barcode_str << 
		//	", smallest barcode: " <<  smallest_barcode <<  
		//	", sallest dist: " << smallest_dist << 
		//	", smallest_count: " << smallest_count << "\n";

		// So the smallest dist has to be unique, otherwise we shall put 
		// 	them in a fil called unknow.

		std::string write_barcode;

		if (smallest_count == 1) {
			write_barcode = smallest_barcode;
			if (smallest_dist == 0) {
				zero_dist_map[write_barcode]++;
			} else if (smallest_dist == 1) {                        
				one_dist_map[write_barcode]++;
			} else {
				higher_dist_map[write_barcode]++;
			}
			match_total++;
				
		} else if (smallest_count > 0) {
			write_barcode = "ambiguous";
			ambiguous_total++;
		} else {
			write_barcode = "no_match";
			no_match_total++;
		}
		barcode_set.insert(write_barcode);
	
		// Adding umi_string to the output file.	
		std::string rword1A;
		if (validUmi) {
			rword1A = rword1 + ":" + umi_str;
		} else {
			rword1A = rword1;
		}
	
		totalcap = updateMaps(write_barcode, lword1, lword2, lword3, lword4, 
			rword1A, rword2, rword3, rword4, totalcap);
			
		//std::cout << "total cap: " << totalcap << "\n";

		if (totalcap > total_allowed) {
			writeMapsToFile();
			// Write all the data in the respective files sequentially
			totalcap = 0;
		}

		distmap[smallest_dist]++;
	}

	// final writing to the files
	writeMapsToFile();
	totalcap = 0;

	for (const auto& it : distmap) {
		std::cout << it.first
				<< ':'
				<<it.second
				<< std::endl;
	}	

	int lcount2 = 1;
	for (const auto& it2 : barcode_set) {
		std::cout << lcount2++ << " " << it2 << "\n";
	}

	//log_detailed.close();
}


void bc_splitter::write_log() {

	// Writing the logs
	const std::string logfile1 = outdirpath + "/frequency_logfile.txt";
	std::ofstream log_freq(logfile1);

	std::setprecision(2);

	unsigned long total_reads = match_total + ambiguous_total + no_match_total;

	double ambiguous_percent = ((double) ambiguous_total / (double) total_reads) * 100;
	double no_match_percent = ((double) no_match_total / (double) total_reads) * 100;

	log_freq << "Total reads: " << total_reads << "\n..................\n";

	log_freq << "Ambiguous:\n";
	log_freq << ".................." << "\n";
	log_freq << "Total ambiguous reads: " << ambiguous_total << " (" << ambiguous_percent << "%)\n\n";

	log_freq << "No match:\n";
    log_freq << ".................." << "\n";
    log_freq << "Total non-match reads: " << no_match_total << " (" << no_match_percent << "%)\n\n";


	for (const auto& lbarcode : barcode_set) {

		if (lbarcode == "ambiguous" || lbarcode == "no_match") {continue;}
	
		unsigned long zero_dist_count = zero_dist_map[lbarcode];
		unsigned long one_dist_count = one_dist_map[lbarcode];
		unsigned long higher_dist_count = higher_dist_map[lbarcode];

		unsigned long total_correct_count = zero_dist_count + one_dist_count + higher_dist_count;
		double zero_dist_percent = ((double)zero_dist_count / (double)total_correct_count) * 100.0;
		double one_dist_percent = ((double)one_dist_count / (double)total_correct_count) * 100.0;
		double higher_dist_percent = 100 - zero_dist_percent - one_dist_percent;
		double barcode_read_percent = ((double)total_correct_count / (double)total_reads) * 100;  

		log_freq << "Barcode: " << lbarcode << "\n";
		log_freq << ".................." << "\n";
		log_freq << "Zero base mismatch: " << zero_dist_percent << "%\n";
		log_freq << "One base mismatch: " << one_dist_percent << "%\n";
		log_freq << "Total read for this barcode: " << total_correct_count << 
			" (percent of total reads: " << barcode_read_percent << "%)\n";
		log_freq << "\n";

	}

	log_freq.close();
}


int main(int argc, char* argv[]) { 

	bc_splitter lbs;

	// To avoid performance hit but keeping a large number of files open 
	// and writing them in parallel, we have decided to maintain a queue
	// of lines according to barcodes. So, it reality, there would be
	// n number of queues where n is the number of barcodes.
	// Once the total size of the queue exceeds a predefined limit such as
	// 2 GB, we shall write the queu to files. We shall open those files
	// at that time sequentially, update them and close them sequentially. 
	// This is a sort of lazy update, I hope that it would provide good 
	// performence.
	
	bool all_set = true;
	try {
		all_set = lbs.parse_args(argc, argv);	
	} catch(std::exception& e) {
		std::cerr << "error: " << e.what() << "\n";
        //lbs.print_help();
        return 1;

	} catch (...) {
		//lbs.print_help();
		return 0;
	}

 	if (!all_set) {
		lbs.print_help();
		return 0;
	}

	lbs.initialize();
	try {
		lbs.split_engine();
	} catch(std::invalid_argument& e) {
        std::cerr << "error: " << e.what() << "\n";
		//lbs.print_help();
		return 1;
    }

	lbs.write_log();
        
    return 0;
}
