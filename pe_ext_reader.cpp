#include "pe_ext_reader.h"


std::vector<std::string> GetSampleList(std::string& sample_list_file_path)
{
    std::ifstream sample_list_file(sample_list_file_path);
    std::vector<std::string> sample_list;
    std::string str{};
    logger.Info("Retrieving a sample list from " + sample_list_file_path);
    if (sample_list_file.is_open())
    {
        while(std::getline(sample_list_file, str))
        {
            if (Util::file_exists(str))
            {
                sample_list.push_back(str);
                logger.Info("  * sample file: " + str);
            }
            else
                logger.Error("  File " + str + " does not exist.");
        }
    }
    sample_list_file.close();
    return sample_list;
}
int PeExtReader::Parse(ContigContainerPtr out_contigContainer, IdMapPtr out_idmap)
{
    auto fasta_input_file = config.GetValue<std::string>("sample_list");
    return Parse(fasta_input_file, out_contigContainer, out_idmap);
}
int PeExtReader::Parse(std::string& fasta_file_path, ContigContainerPtr out_contigContainer, IdMapPtr out_idmap)
{
    if (Util::get_file_ext(fasta_file_path) != "fasta")
    {
        auto sample_list = GetSampleList(fasta_file_path);
        return Parse(sample_list, out_contigContainer, out_idmap);
    }
    else
    {
        if (Util::file_exists(fasta_file_path))
        {
            ParseFile(fasta_file_path, *out_contigContainer, out_idmap, 0);   
            SaveContigContainer(*out_contigContainer);    
            return 1;
        }
        else
        {
            logger.Error("File '" + fasta_file_path + "' is not found");
            return -1;
        }
    }
}

int PeExtReader::Parse(std::vector<std::string>& fasta_path_list, ContigContainerPtr out_contigContainer, IdMapPtr out_idmap)
{
    unsigned int fasta_file_num = 0;
    for(auto fasta_file_path : fasta_path_list)
    {
        ParseFile(fasta_file_path, *out_contigContainer, out_idmap, fasta_file_num++);
    }
    SaveContigContainer(*out_contigContainer);
    return fasta_file_num;
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

    std::ifstream contig_file(fasta_file_path);
    std::unordered_set<std::string> unique_contigs;
    
    if(contig_file.is_open())
    {
        std::string line;
        while(contig_file)
        {
            std::getline(contig_file, line);
            Util::trim_str(line);
            if (line.size() > 0 && line[0] == '>')
            {
                ExtContig econtig(line);    
                econtig.m_sample_id = fasta_file_num; 
                econtig.m_contig_id = out_contigContainer.size(); 

                (*out_idmap)[econtig.m_contig_id] = econtig.m_sample_id;
                
                while(contig_file.peek() != '>' && !contig_file.eof())                                      
                {
                    std::getline(contig_file, line);
                    Util::trim_str(line);
                    econtig.Seq.append(line);
                }                
                
                auto unext_len = econtig.Seq.size() - econtig.m_bwd_ext_sz - econtig.m_fwd_ext_sz;   
                if (remove_duplicates && unique_contigs.find(econtig.Name) != unique_contigs.end() && 
                    unext_len >= length_threshold)
                {
                    if (econtig.m_bwd_ext_sz >= min_adj_overlap)
                    {
                        econtig.Hdr = econtig.SampleName + ":" + econtig.Name + ":" + Util::to_str(econtig.m_bwd_ext_sz) + ":0";
                        Util::get_first_chars(econtig.Seq, flank_size + econtig.m_bwd_ext_sz);
                        econtig.m_is_tag = true;
                        econtig.m_unext_len = econtig.Seq.size() - econtig.m_bwd_ext_sz - econtig.m_fwd_ext_sz;                              
                        out_contigContainer.push_back(econtig.Clone());
                    }
                    if(econtig.m_fwd_ext_sz >= min_adj_overlap)
                    {
                        econtig.Hdr = econtig.SampleName + ":" + econtig.Name + ":0:" + Util::to_str(econtig.m_fwd_ext_sz);
                        Util::get_last_chars(econtig.Seq, flank_size + econtig.m_fwd_ext_sz);
                        econtig.m_unext_len = econtig.Seq.size() - econtig.m_bwd_ext_sz - econtig.m_fwd_ext_sz;       
                        econtig.m_is_tag = true;
                        out_contigContainer.push_back(econtig.Clone());
                    }               
                }
                else
                {
                    if (econtig.m_bwd_ext_sz > 0 && econtig.m_bwd_ext_sz < min_adj_overlap)
                    {
                        Util::remove_first_chars(econtig.Seq, econtig.m_bwd_ext_sz);
                    }
                    if (econtig.m_fwd_ext_sz > 0 && econtig.m_fwd_ext_sz < min_adj_overlap)
                    {
                        Util::remove_last_chars(econtig.Seq, econtig.m_fwd_ext_sz);
                    }            
                    econtig.m_unext_len = unext_len;        
                    if (remove_duplicates)
                        unique_contigs.insert(econtig.Name);

                    out_contigContainer.push_back(econtig);  
                }                        
            }
        }
        logger.Debug("Contigs size: " + Util::to_str(out_contigContainer.size()));
    }
    else
    {
        logger.Error("File " + fasta_file_path + " does not exist.");
    }
    contig_file.close();
}
/* std::string ExtContig::ToString()
{
    return GetHeader() + "\n" + this.Seq + "\n";
}

stringptr ExtContig::GetHeader()
{
    return std::make_unique<std::string>(">" + Util::to_str<int>(-1) + ":" + SampleName + ":" + 
        Util::to_str<int>(-1) + ":" + Name + ":" + Util::to_str<unsigned int>(m_bwd_ext_sz) + ":" +
        Util::to_str<unsigned int>(m_fwd_ext_sz) + ":" + Util::to_str(m_is_tag));
} */
void PeExtReader::SaveContigContainer(ContigContainer& ccon)
{
    auto prefix = config.GetValue<std::string>("graph_pref");
    auto out_dir = config.GetValue<std::string>("output_dir");
    auto processed_input_fasta_path = out_dir + "/" + prefix + "_processed_input.fasta";

    std::ofstream processed_input_fasta_file(processed_input_fasta_path);
    for(auto& econtig : ccon)
    {
        processed_input_fasta_file << ">" + Util::to_str<int>(-1) + ":" + econtig.SampleName + ":" + 
        Util::to_str<int>(-1) + ":" + econtig.Name + ":" + Util::to_str<unsigned int>(econtig.m_bwd_ext_sz) + ":" +
        Util::to_str<unsigned int>(econtig.m_fwd_ext_sz) + ":" + Util::to_str(econtig.m_is_tag) + "\n";
        unsigned int i = 0;
        const unsigned int fold = 80;
        
        while(i < econtig.Seq.size())
        {
            unsigned int rest = econtig.Seq.size() - i;
            rest = std::min(rest, fold);
            processed_input_fasta_file << econtig.Seq.substr(i, rest) << "\n";
            i += fold;
        }
    }
    processed_input_fasta_file.close();

    logger.Info("Saved processed fasta file to " + processed_input_fasta_path);
}

