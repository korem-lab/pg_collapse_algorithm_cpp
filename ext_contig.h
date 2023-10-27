#ifndef EXT_CONTIG_H
#define EXT_CONTIG_H

#include <memory>
#include <string>
#include <vector>
#include "aux_functions.h"

class ExtContig
{
    public:
        ExtContig(){}
        ExtContig(const std::string& hdr, unsigned int contig_id)
        {
            parseHdr(hdr);
            m_contig_id = contig_id;
        }
        stringptr GetName() {return std::make_unique<std::string>(m_name);}
        int GetId() {return m_contig_id;}
        stringptr GetSampleName() {return std::make_unique<std::string>(m_sample_name);}
        stringptr GetSeq() {return std::make_unique<std::string>(m_seq);}
        unsigned int GetFwdExt() {return m_bwd_ext_sz; }
        unsigned int GetBwdExt() {return m_bwd_ext_sz; }  
        unsigned int GetUnExtLen() {return m_unext_len;}
        unsigned int GetLen() {return m_seq.size();}
        bool GetIsTag() {return m_is_tag;}
        stringptr GetHeader();
        stringptr ToString();
        friend class PeExtReader;

    private:
        void parseHdr(const std::string& hdr);
        std::string m_hdr;
        std::string m_name;
        std::string m_sample_name;
        std::string m_seq;
        unsigned int m_fwd_ext_sz;
        unsigned int m_bwd_ext_sz; 
        int m_sample_id;
        int m_contig_id = -1;
        unsigned int m_unext_len;
        bool m_is_tag = false;
        static unsigned int contig_id_generator;      
};


typedef std::vector<ExtContig> ContigContainer;
typedef std::shared_ptr<ContigContainer> ContigContainerPtr;

#endif