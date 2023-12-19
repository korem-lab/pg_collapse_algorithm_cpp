
#include "global_copan.h"
#include "async_logger.h"
#include <string>
#include <iostream>
#include "pe_ext_reader.h"
#include "ext_contig.h"
#include "alignment.h"
#include "gluepoints.h"
#include "copangraph.h"

 Util::ConfigReader config;
 Util::AsyncLogger logger;


int main(int argc, char *argv[])
{
    auto config_name = Util::ParseConfigName(argc, argv);
    //Read from config file and initialize a logger
    Util::InitializeServices(config_name);    
    
    PeExtReader per{};
    ContigContainerPtr ccptr(new ContigContainer());
    IdMapPtr imptr(new IdMap());
    // Reads fasta file and populates contig container and idmap
    auto num_samples = per.Parse(ccptr, imptr);
    
    auto alignments = compute_alignments(ccptr);   

    GluepointMapPtr gmp = nullptr;
    auto node_intervals = cpp_gen_gp(alignments, ccptr);

    auto copan = CoPanGraph(node_intervals, ccptr, imptr, num_samples);
    copan.build_graph();        
    copan.output_graph();

    Util::Close();
}   
   