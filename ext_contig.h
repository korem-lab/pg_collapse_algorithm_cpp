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
        ExtContig(const std::string& hdr)
        {
            parseHdr(hdr);
        }
        std::string Hdr;
        std::string Name;
        std::string SampleName;
        std::string Seq;
        
        int GetId() const {return m_contig_id;}        
        unsigned int GetFwdExt() const {return m_fwd_ext_sz; }
        unsigned int GetBwdExt() const {return m_bwd_ext_sz; }  
        unsigned int GetUnExtLen() {return m_unext_len;}
        unsigned int GetLen() const {return Seq.size();}
        bool GetIsTag() {return m_is_tag;}
        friend class PeExtReader;

    private:
        void parseHdr(const std::string& hdr);

        unsigned int m_fwd_ext_sz;
        unsigned int m_bwd_ext_sz; 
        int m_sample_id;
        int m_contig_id = -1;
        unsigned int m_unext_len;
        bool m_is_tag = false;
};


typedef std::vector<ExtContig> ContigContainer;
typedef std::shared_ptr<ContigContainer> ContigContainerPtr;

#endif