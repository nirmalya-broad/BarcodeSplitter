/*
    load_and_search.cpp
	
	By Nirmalya Bandyopadhyay 
	adapted from 
	---
    By Stephen Holiday 2011
    http://stephenholiday.com
    (Exception, Distance Algorithm by Anders Sewerin Johansen)

    The code is under the [Apache 2.0](http://www.apache.org/licenses/LICENSE-2.0) license.

    This is a utility for loading a BKTree from a serialized file.
    
    ./bkSearch [serialized tree] [word] [max edxit distance]
    ./bkSearch bkTree.dat word 3

*/

#include <iostream>
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


#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

int distance(std::string source, std::string target) {

    const int n = source.length();
    const int m = target.length();
    if (n == 0) {
        throw std::invalid_argument("Length of source is zero.");
    }
    if (m == 0) {
        throw std::invalid_argument("Length of target is zero.");
    }

    if (m !=n )
        throw std::invalid_argument("Source and target have different length.");

    int ldist = 0;

    for (int j = 0; j < n; j++) {
        if (source[j] != target[j]) {
            ldist++;
        }
    }

    return ldist;
}	


unsigned long 
updateMaps(std::string& barcode_str, 
	std::string& lword1, std::string& lword2, 
	std::string& lword3, std::string& lword4, 
	std::string& rword1, std::string& rword2, 
	std::string& rword3, std::string rword4, 
	std::map<std::string, std::vector<std::string>>& lQueueMap,
    std::map<std::string, std::vector<std::string>>& rQueueMap,
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

void writeMapsToFile(std::map<std::string, std::vector<std::string>>& lQueueMap,
		std::map<std::string, std::vector<std::string>>& rQueueMap,
		std::set<std::string>& outfile_set, std::string const& outdirpath) {
	
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

} 

int main(int argc, char* argv[]) { 
    BKTree<std::string> tree;
    {
        std::ifstream iff(argv[1]);
        boost::archive::text_iarchive iar(iff);
        
        iar >> tree;
    }
    
   
	// argv[2] contains read 1 containing all the barcode and
	// argv[3] contains read 2
	// argv[4] contains the cutoff for mismatch
	// argv[5] outdir for output files.
	// any mismatches would be thrown to a file called unknown.
	
	std::string file1_str = argv[2];
	std::string file2_str = argv[3];
	int cutoff = atoi(argv[4]);
	std::string outdirpath = argv[5];


	std::set<std::string> barcode_set;
	std::set<std::string> outfile_set;

	std::map<int, int> distmap;
	for (int j = 0; j <= cutoff; j++) {
		distmap[j] = 0;
	}	


	struct stat st = {0};

	if (stat(outdirpath.c_str(), &st) == -1) {
		mkdir(outdirpath.c_str(), 0700);
	}


	// To avoid performance hit but keeping a large number of files open 
	// and writing them in parallel, we have decide to maintain a queue
	// of lines according to barcode. So, it reality, there would be
	// n number of queues where n is the number of barcodes.
	// Once the total size of the queue reaches a predefined number such as
	// 2 GB, we shall write the 75% or so to files. We shall open those files
	// at that time sequentially. This is a sort of lazy update, I hope that
	// it would provide good performence.

	std::map<std::string, std::vector<std::string>> lQueueMap;
	std::map<std::string, std::vector<std::string>> rQueueMap;

	int barcode_start = 6;
	int barcode_size = 6;

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

	const unsigned long GB_SIZE = 1024*1024*1024;
	int allowed_GB = 2;
	unsigned long total_allowed = allowed_GB * GB_SIZE;

	unsigned long totalcap = 0;

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
	
		std::vector<std::string> results = tree.find(barcode_str, cutoff);

		// calculate the minimum distance between the target and references

		int smallest_dist = cutoff;
		for (auto const& val : results) {
			int ldist = distance(val, barcode_str);
			if (ldist < smallest_dist) {
				smallest_dist = ldist;
			}
		}

		if (smallest_dist == 0) {
			barcode_set.insert(barcode_str);
			totalcap = updateMaps(barcode_str, lword1, lword2, lword3, lword4, rword1, rword2, 
				rword3, rword4, lQueueMap, rQueueMap, totalcap);

		}


		if (totalcap > total_allowed) {
			writeMapsToFile(lQueueMap, rQueueMap, outfile_set, outdirpath);
			// Write all the data in the respective files sequentially
			lQueueMap.clear();
			rQueueMap.clear();
			totalcap = 0;
		}

		distmap[smallest_dist]++;
	}

	// final writing to the files
	if (totalcap > total_allowed) {
		writeMapsToFile(lQueueMap, rQueueMap, outfile_set, outdirpath);
		// Write all the data in the respective files sequentially
		lQueueMap.clear();
		rQueueMap.clear();
		totalcap = 0;
	}


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
        
    return 0;
}
