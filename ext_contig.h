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
        uint32_t GetFwdExt() const {return m_fwd_ext_sz; }
        uint32_t GetBwdExt() const {return m_bwd_ext_sz; }  
        uint32_t GetUnExtLen() {return m_unext_len;}
        uint32_t GetLen() const {return Seq.size();}
        bool GetIsTag() {return m_is_tag;}
        friend class PeExtReader;

    private:
        void parseHdr(const std::string& hdr);

        uint32_t m_fwd_ext_sz;
        uint32_t m_bwd_ext_sz; 
        int m_sample_id;
        int m_contig_id = -1;
        uint32_t m_unext_len;
        bool m_is_tag = false;
};


typedef std::vector<ExtContig> ContigContainer;
typedef std::shared_ptr<ContigContainer> ContigContainerPtr;

#endif