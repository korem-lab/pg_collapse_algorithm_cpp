#ifndef PE_EXT_READER_H
#define PE_EXT_READER_H
#include <string>
#include <memory>
#include "aux_functions.h"
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include "global_copan.h"
#include "async_logger.h"
#include "ext_contig.h"
#include <algorithm>

typedef std::unordered_map<int, int> IdMap;
typedef std::shared_ptr<IdMap> IdMapPtr;

class PeExtReader
{
    public:
        unsigned int  Parse(std::string& fasta_file_path, ContigContainerPtr out_contigContainer, IdMapPtr out_idmap);
        unsigned int  Parse(std::vector<std::string>& fasta_path_list, ContigContainerPtr out_contigContainer, IdMapPtr out_idmap);
    private:
        void ParseFile(std::string& fasta_file_path, ContigContainer& out_contigContainer, IdMapPtr out_idmap, int fasta_file_num);
        void SaveContigContainer(ContigContainer& ccon);
};

#endif