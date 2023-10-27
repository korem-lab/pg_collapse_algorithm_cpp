
#include "ext_contig.h"


void ExtContig::parseHdr(const std::string& hdr)
{
    auto hdr_parts = Util::split(hdr, ':');
    if (hdr_parts->size() >= 4)
    {
        auto sname = hdr_parts->at(0);
        m_sample_name = sname.size() > 0 && sname[0] == '>' ? sname.substr(1, sname.size() - 1) : sname;
        m_name = hdr_parts->at(1);
        m_bwd_ext_sz = Util::convert_to_int(hdr_parts->at(2));
        m_fwd_ext_sz = Util::convert_to_int(hdr_parts->at(3));     
        m_unext_len = 0;  
        if (m_contig_id > 84460)
        {
            bool stop = true;
        }
    }
}
stringptr ExtContig::ToString()
{
    m_hdr = *GetHeader();
    return std::make_unique<std::string>(m_hdr + "\n" + m_seq + "\n");
}

stringptr ExtContig::GetHeader()
{
    return std::make_unique<std::string>(">" + Util::convert_to_string(-1) + ":" + m_sample_name + ":" + 
        Util::convert_to_string(-1) + ":" + m_name + ":" + Util::convert_to_string(m_bwd_ext_sz) + ":" +
        Util::convert_to_string(m_fwd_ext_sz) + ":" + Util::convert_to_string(m_is_tag));
}
