#include "aux_functions.h"
#include "ext_contig.h"


void ExtContig::parseHdr(const std::string& hdr)
{
    auto hdr_parts = Util::split(hdr, ':');
    if (hdr_parts->size() >= 4)
    {
        auto sname = hdr_parts->at(0);
        SampleName = sname.size() > 0 && sname[0] == '>' ? sname.substr(1, sname.size() - 1) : sname;
        Name = hdr_parts->at(1);
        m_bwd_ext_sz = Util::convert_to_int(hdr_parts->at(2));
        m_fwd_ext_sz = Util::convert_to_int(hdr_parts->at(3));     
        m_unext_len = 0;  
        
    }
}

