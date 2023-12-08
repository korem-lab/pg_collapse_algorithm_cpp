
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

void InitializeServices(const char* config_name)
{
    try
    {
        config.Initialize(config_name);

        auto log_file_dir = config.GetValue<std::string>("log_file_dir");
        auto app_name = config.GetValue<std::string>("app_name");

        logger.Start(log_file_dir, app_name);
        logger.Info("Initializing... ");
        logger.Info(*config.GetAllKeyValues());
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        logger.Stop();
    }
}

const char* ParseConfigName(int argc, char *argv[])
{
    if (argc > 1)
        return argv[1];
    return "config.ini";
}

 
int main(int argc, char *argv[])
{
    auto config_name = ParseConfigName(argc, argv);
    InitializeServices(config_name);
    
    auto fasta_input_file = config.GetValue<std::string>("sample_list");
    
    PeExtReader per{};
    ContigContainerPtr ccptr(new ContigContainer());
    IdMapPtr imptr(new IdMap());
    
    auto num_samples = per.Parse(fasta_input_file, ccptr, imptr);
    
    auto rg_pref = config.GetValue<std::string>("graph_pref");

    auto alignments = compute_alignments(ccptr);
    
    auto max_separation = config.GetValue<unsigned int>("max_separation");
    GluepointMapPtr gmp = nullptr;
    auto node_intervals = cpp_gen_gp(alignments, ccptr, max_separation);

    auto copan = CoPanGraph(node_intervals, ccptr, imptr, num_samples);
    copan.build_graph();

    logger.Info("writing to files");
    auto output_dir = config.GetValue<std::string>("output_dir");
    auto file_name = output_dir + "/" + rg_pref;
    auto gfa_file = file_name+ config.GetValue<std::string>("gfa_file_ext");
    auto fasta_file = file_name + config.GetValue<std::string>("fasta_file_ext");
    auto ncolor_file = file_name + config.GetValue<std::string>("node_color_file_ext");
    auto ecolor_file = file_name + config.GetValue<std::string>("edge_color_file_ext");
    
    copan.output_graph(gfa_file, fasta_file ,ncolor_file, ecolor_file);

    logger.Info("Complete!"); 
    logger.Stop();
}  
   