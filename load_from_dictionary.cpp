/*
    load_and_search.cpp
	
    ./bkLoad [from] [to]
    ./bkLoad /usr/share/dict/words wordTree.dat

*/

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib> 
#include <cstring>

#include "BKTree.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/regex.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

class load_from_dict{

	public:
	bool parse_args(int argc, char* argv[]);
	bool build_tree_simple();
	bool build_tree_complex();
	void save_tree();
	std::string& get_type();
	void print_help();
	load_from_dict();
	~load_from_dict();

	private:
	BKTree<std::string> tree;
	std::string infile;
	std::string outfile;
	std::string ltype;
    po::options_description desc;

};


void load_from_dict::print_help() {
	std::cout << desc << "\n";
    std::cout << "usage: bkLoad -i <infile> -o <outfile> -t <type>\n\n";
}

load_from_dict::load_from_dict() {
	tree = BKTree<std::string>();
}

std::string& load_from_dict::get_type() {
	return ltype;
}

load_from_dict::~load_from_dict() {
	// Nothing yet	
}

bool load_from_dict::build_tree_simple() {
	std::ifstream words(infile);
	std::string lstr;

	if (!words.is_open()) {
    	std::cerr << "The infile cannot be open!\n";
		return false;
	} else {
		while (std::getline(words, lstr)) {
			tree.insert(lstr);
		}
	}
	std::cout<< "Loaded " <<tree.size()<< " entries" << std::endl;
	return true;
}

void load_from_dict::save_tree() {
    
	std::ofstream ofs(outfile);
    boost::archive::text_oarchive oa(ofs);
    oa << tree;
	ofs.close();
	return;
}


bool load_from_dict::build_tree_complex() {
 std::ifstream words (infile);
    std::string str;
    
    boost::regex expr ("NNNNNN(\\w{6})");
	boost::smatch what;	
    
    if (!words.is_open()) {
    	std::cerr << "The infile cannot be open!\n";
		return false;
    }
    else {
        while(std::getline(words, str)) {
				
			// Extract the six bases from each of the line after NNNNNN
			bool res = boost::regex_search(str, what, expr);
			if (res) {				
				std::string lstr = what[1];
        		tree.insert(lstr);
			}
		}
    }
    
    std::cout<<"Loaded "<<tree.size()<< " entries"<<std::endl;
	return true;
}

bool load_from_dict::parse_args(int argc, char* argv[]) {
	 bool all_set = true;

    desc.add_options()
        ("help,h", "produce help mesage")
        ("infile,i", po::value<std::string>(&infile), "Input file")
        ("outfile,o", po::value<std::string>(&outfile), "Output file")
		("type,t", po::value<std::string>(&ltype), "Input file type")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        return 0;
    }

    if (!vm.count("infile")) {
        std::cout << "Infile not set.\n";
        all_set = false;
    } else {
        std::cout << "Infile is set to: " << infile << ".\n";
    }

    if (!vm.count("outfile")) {
        std::cout << "Outfile not set.\n";
        all_set = false;
    } else {
        std::cout << "Outfile is set to: " << outfile << ".\n";
    }

	if (!vm.count("type")) {
		std::cout << "Type not set.\n";
	} else {
		std::cout << "Type set to: " << ltype << ".\n";
	}

	return all_set;
}

int main(int argc, char* argv[]) { 
	load_from_dict ldict = load_from_dict();
	bool all_set = true;
	
	po::options_description desc("Allowed options");
	try {

		all_set = ldict.parse_args(argc, argv);
	} catch(std::exception& e) {
       	std::cerr << "error: " << e.what() << "\n";
		ldict.print_help();
       	return 1;
   	}
   	catch(...) {
       	std::cerr << "Exception of unknown type!\n";
		ldict.print_help();
   	}

	if (!all_set) {
		ldict.print_help();
		return 0;
	}

	std::string ltype = ldict.get_type();
	if (ltype.compare("simple") == 0) {
		ldict.build_tree_simple();
	} else if (ltype.compare("complex") == 0) {
		ldict.build_tree_complex();
	} else {
		std::cerr << "Invalid type: " << ltype << ".\n";
	}
	ldict.save_tree();
        
    return 0;
}

