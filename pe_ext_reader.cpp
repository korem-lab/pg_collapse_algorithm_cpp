#include "pe_ext_reader.h"

unsigned int PeExtReader::Parse(std::string& fasta_file_path, ContigContainerPtr out_contigContainer, IdMapPtr out_idmap)
{
    ParseFile(fasta_file_path, *out_contigContainer, out_idmap, 0);   
    SaveContigContainer(*out_contigContainer);
    return 1;
}

unsigned int PeExtReader::Parse(std::vector<std::string>& fasta_path_list, ContigContainerPtr out_contigContainer, IdMapPtr out_idmap)
{
    unsigned int fasta_file_num = 0;
    for(auto fasta_file_path : fasta_path_list)
    {
        ParseFile(fasta_file_path, *out_contigContainer, out_idmap, fasta_file_num++);
    }
    SaveContigContainer(*out_contigContainer);
    return fasta_file_num - 1;
}

void PeExtReader::ParseFile(std::string& fasta_file_path, ContigContainer& out_contigContainer, IdMapPtr out_idmap, int fasta_file_num)
{
    logger.Info("Parsing file " + fasta_file_path);   

    auto length_threshold = config.GetValue<int>("length_threshold");
    auto remove_duplicates = config.GetValue<bool>("remove_duplicates");
    auto min_overlap = config.GetValue<int>("min_overlap");
    auto epsilon = config.GetValue<int>("epsilon");
    auto flank_size = config.GetValue<int>("flank_size");
    auto min_adj_overlap = min_overlap - epsilon;

    std::ifstream configFile(fasta_file_path);
    std::unordered_set<std::string> unique_contigs;
    
    if(configFile.is_open())
    {
        std::string line;
        while(configFile)
        {
            std::getline(configFile, line);
            Util::trim_str(line);
            if (line.size() > 0 && line[0] == '>')
            {
                ExtContig econtig(line, out_contigContainer.size());    
                econtig.m_sample_id = fasta_file_num; 
                
                (*out_idmap)[econtig.m_contig_id] = econtig.m_sample_id;
                
                while(configFile.peek() != '>' && !configFile.eof())                                      
                {
                    std::getline(configFile, line);
                    Util::trim_str(line);
                    econtig.m_seq.append(line);
                }

                econtig.m_unext_len = econtig.m_seq.size() - econtig.m_bwd_ext_sz - econtig.m_fwd_ext_sz;

                if (remove_duplicates && unique_contigs.find(econtig.m_name) != unique_contigs.end() && 
                    econtig.m_unext_len >= length_threshold)
                {
                    if (econtig.m_bwd_ext_sz >= min_adj_overlap)
                    {
                        econtig.m_hdr = econtig.m_sample_name + ":" + econtig.m_name + ":" + Util::convert_to_string(econtig.m_bwd_ext_sz) + ":0";
                        Util::get_first_chars(econtig.m_seq, flank_size + econtig.m_bwd_ext_sz);
                        econtig.m_is_tag = true;
                        out_contigContainer.push_back(econtig);
                    }
                    else if(econtig.m_fwd_ext_sz >= min_adj_overlap)
                    {
                        econtig.m_hdr = econtig.m_sample_name + ":" + econtig.m_name + ":0:" + Util::convert_to_string(econtig.m_fwd_ext_sz);
                        Util::get_last_chars(econtig.m_seq, flank_size + econtig.m_fwd_ext_sz);
                        econtig.m_is_tag = true;
                        out_contigContainer.push_back(econtig);
                    }               
                }
                else
                {
                    if (econtig.m_bwd_ext_sz > 0 && econtig.m_bwd_ext_sz < min_adj_overlap)
                    {
                        Util::remove_first_chars(econtig.m_seq, econtig.m_bwd_ext_sz);
                    }
                    if (econtig.m_fwd_ext_sz > 0 && econtig.m_fwd_ext_sz < min_adj_overlap)
                    {
                        Util::remove_last_chars(econtig.m_seq, econtig.m_fwd_ext_sz);
                    }                    

                    if (remove_duplicates)
                        unique_contigs.insert(econtig.m_name);

                    out_contigContainer.push_back(econtig);  
                }
            }
        }
        logger.Debug("Contigs size: " + Util::convert_to_string(out_contigContainer.size()));
    }
    configFile.close();
}

void PeExtReader::SaveContigContainer(ContigContainer& ccon)
{
    auto prefix = config.GetValue<std::string>("graph_pref");
    auto out_dir = config.GetValue<std::string>("output_dir");
    auto processed_input_fasta_path = out_dir + "/" + prefix + "_processed_input.fasta";

    std::ofstream processed_input_fasta_file(processed_input_fasta_path);
    for(auto& econtig : ccon)
    {
        processed_input_fasta_file << (*econtig.ToString());
    }
    processed_input_fasta_file.close();

    logger.Info("Saved processed fasta file to " + processed_input_fasta_path);
}

