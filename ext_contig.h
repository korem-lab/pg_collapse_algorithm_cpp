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
        
        int GetId() {return m_contig_id;}        
        unsigned int GetFwdExt() {return m_bwd_ext_sz; }
        unsigned int GetBwdExt() {return m_bwd_ext_sz; }  
        unsigned int GetUnExtLen() {return m_unext_len;}
        unsigned int GetLen() {return Seq.size();}
        bool GetIsTag() {return m_is_tag;}
        ExtContig Clone()
        {
            ExtContig c;
            c.Name = Name;
            c.SampleName = SampleName;
            c.Seq = Seq;
            c.Hdr = Hdr;
            c.m_bwd_ext_sz = m_bwd_ext_sz;
            c.m_fwd_ext_sz = m_fwd_ext_sz;
            c.m_sample_id = m_sample_id;
            c.m_contig_id = m_contig_id;
            c.m_unext_len = m_unext_len;
            c.m_is_tag = m_is_tag;            
            return c;
        }
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