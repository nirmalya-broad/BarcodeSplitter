#include <stxxl/sorter>
#include <stxxl/stats>
#include <limits>
#include <fstream>
#include <cstdio>
#include <iostream>

#include <regex>

/* Parses from the simple input that Jonathan provides */

struct jp_record {

	long long start_c;
	long long end_c;
	char strand;
	char is_paired;
	std::string umi_str;
	std::string rname;


	jp_record(long long _start_c, long long _end_c, char _strand, 
		char _is_paired, std::string _umi_str, std::string _rname):
		start_c(_start_c), 
		end_c(_end_c), 
		strand(_strand), 
		is_paired(_is_paired), 
		umi_str(_umi_str),
		rname(_rname)
	{}

	jp_record() {}

};

struct jpr_Comparator {
	bool operator () (const jp_record& a, const jp_record& b) const {
		int umi_c = a.umi_str.compare(b.umi_str);
		if (umi_c < 0) {
			return true;
		} else if (umi_c == 0) {
			if (a.strand < b.strand) {
				return true;
			} else if (a.strand == b.strand) {
				if (a.start_c < b.start_c) {
					return true;
				}
			}
		}
		return false;
	}

	jp_record min_value() const {
		jp_record dummy;
		dummy.start_c = std::numeric_limits<long long>::min();
		dummy.strand = '+';
		dummy.umi_str = "AAAAAA";

		return dummy;
	}

	jp_record max_value() const {

		jp_record dummy;
		dummy.start_c = std::numeric_limits<long long>::max();
		dummy.strand = '-';
		dummy.umi_str = "TTTTTT";

		return dummy;

	}
};


class UmiNorm {
	public:
	int parse_jpr_file(std::string infile_str);
	void sort_jp_records();
	void print_jp_records(std::string outfile_str, std::string freq_file_str);
	stxxl::sorter<jp_record, jpr_Comparator, 1*1024*1024> jpr_sorter_vec;
	std::map<std::string, int> umiCountMap;

	UmiNorm(int mem_size); 

	private:
	int lmem_size = 64;

};

UmiNorm::UmiNorm(int mem_size) :
	lmem_size(mem_size),
	jpr_sorter_vec(jpr_Comparator(), mem_size * 1024 * 1024)
{
}


int UmiNorm::parse_jpr_file(std::string infile_str) {
	std::regex expr("(\\d+)\\.\\s+(\\d+)\\.\\s+(\\S)\\s+(\\w)\\s+(\\S+)\\:(\\w{6})");
	std::smatch what;
	
	std::string lstr;
	std::ifstream infile(infile_str);

	while (std::getline(infile, lstr)) {
	
		bool res = std::regex_search(lstr, what, expr);
		//std::cout << lstr << "\n";

		if (res) {
			std::string start_str = what[1];
			std::string end_str = what[2];
			std::string strand_str = what[3];
			std::string is_paired_str = what[4];
			std::string rname = what[5];
			std::string umi_str = what[6];
		
			long long lstart = std::stoll(start_str);
			long long lend = std::stoll(end_str);
			char lstrand;
			if (0 == strand_str.compare("+")) {
				lstrand = '+';
			} else if (0 == strand_str.compare("-")) {
				lstrand = '-';
			}

			char is_paired;

			if (0 == is_paired_str.compare("Y") || 0 == is_paired_str.compare("y")) {
				is_paired = 'Y';
			} else if (0 == is_paired_str.compare("N") || 0 == is_paired_str.compare("n")){
				is_paired = 'N';
			}

			//std::cout << lstart << " " << lend << " " << lstrand << " " << is_paired 
			//	<< " " << umi_str << "\n";

			jp_record ljp_record = {
				lstart,
				lend,
				lstrand,
				is_paired,
				umi_str,
				rname
			};
			
			jpr_sorter_vec.push(ljp_record);
			umiCountMap[umi_str]++;
		}
	}

	return 0;
}

void UmiNorm::sort_jp_records() {
	jpr_sorter_vec.sort();
}


void UmiNorm::print_jp_records(std::string outfile_str, std::string freq_file_str) {

	std::ofstream outfile(outfile_str);
	std::ofstream freq_file(freq_file_str);
	std::string last_umi = "";
	char last_strand = '.';

	while (!jpr_sorter_vec.empty()) {
		if ((*jpr_sorter_vec).umi_str.compare(last_umi) != 0) {
			outfile << "----------------------------------------\n";
		} else if ((*jpr_sorter_vec).strand != last_strand) {
			outfile << "............\n";
		}

		last_umi = (*jpr_sorter_vec).umi_str;
		last_strand = (*jpr_sorter_vec).strand;

        outfile << (*jpr_sorter_vec).umi_str << " " 
			<< (*jpr_sorter_vec).strand << " "
			<< (*jpr_sorter_vec).start_c 
			<< " " << (*jpr_sorter_vec).end_c 
			<< " " << (*jpr_sorter_vec).rname << "\n";

        ++jpr_sorter_vec;
    }
		
	outfile.close();

	std::vector<std::pair<int, std::string>> items;

	for(auto& rec : umiCountMap) {
		//freq_file << rec.first << "\t" << rec.second << "\n";
		std::pair<int, std::string> temp_pair;
		temp_pair.first = rec.second;
		temp_pair.second = rec.first;
		items.push_back(temp_pair);
		
	}

	std::sort(items.begin(), items.end(), [](std::pair<int, std::string> a,
		std::pair<int, std::string>	b) {
		return b.first < a.first ||
			b.first == a.first && a.second < b.second;
	});

	for (auto& rec2 : items) {
		freq_file  << rec2.first << "\t" << rec2.second << "\n";
		std::cout  << rec2.first << "\t" << rec2.second << "\n";
	}

	freq_file.close();
}


int main(int argc, char* argv[]) {

	std::string infile_str = argv[1];
	std::string outfile_str = argv[2];
	
	//std::cout << argv[1] << " " << argv[2] << "\n";
	int mem_size = std::atoi(argv[3]);
	std::string freq_file_str = argv[4];

	UmiNorm umi_obj(mem_size);
	umi_obj.parse_jpr_file(infile_str);
	umi_obj.sort_jp_records();
	umi_obj.print_jp_records(outfile_str, freq_file_str);

	return 0;

}

