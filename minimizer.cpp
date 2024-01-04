#include "minimizer.h"
#include <fstream>

using namespace mzr;

uint8_t mzr::char_to_num(char c)
{
    switch (c)
    {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    default:
        return 255;
    }
}
uint8_t mzr::char_to_num_complement(char c)
{
    switch (c)
    {
    case 'A':
        return 3;
    case 'C':
        return 2;
    case 'G':
        return 1;
    case 'T':
        return 0;
    default:
        return 255;
    }
}

char mzr::num_to_char(uint8_t b)
{
    switch (b)
    {
    case 0:
        return 'A';
    case 1:
        return 'C';
    case 2:
        return 'G';
    case 3:
        return 'T';
    default:
        return 'N';
    }
}
    
std::string mzr::unpack_kmer(uint32_t packed_kmer, uint32_t kmer_len)
{
    uint32_t mask = 0x00000003;
    int8_t shift = kmer_len * 2 - 2;
    std::string kmer;
    while (shift >= 0)
    {
        uint32_t base = mask & (packed_kmer >> shift);
        kmer = kmer + mzr::num_to_char(base);
        shift -= 2;
    }
    return (kmer);
       
}

void mzr::count_kmers(const std::string &s, std::unordered_map<uint32_t, uint32_t> &kc, uint8_t k)
{
    KmerGenerator kmer_gen(s, k, 0x3FFFFFFF, true);
    while (!kmer_gen.empty())
    {
        Kmer kmer = kmer_gen.get_kmer();
        kc[kmer.seq]++;    
    }
}

std::unordered_set<uint32_t> mzr::get_high_frequency_kmers(std::unordered_map<uint32_t, uint32_t> kc)
{
    std::unordered_set<uint32_t> high_freq_kmers;   
    auto percentile = 1 - config.GetValue<double>("high_freq_kmer_filter");

    if (percentile < 1) // if percentile == 1 then don't filter
    {
        std::vector<uint64_t> freqs;
        freqs.reserve(kc.size());
        for (auto k : kc)
        {
            freqs.push_back(k.second);
        } 
        std::sort(freqs.begin(), freqs.end());
        size_t p_idx = percentile * freqs.size() + 1;
        high_freq_kmers.reserve(freqs.size() - p_idx + 1);
        if (p_idx >= freqs.size())
        {
            p_idx = freqs.size() - 1;
        }
        uint64_t cut_off = freqs[p_idx];
        for (auto k : kc)
        {
            if (k.second >= cut_off)
            {
                high_freq_kmers.insert(k.first);
            }
        }
    }
    return high_freq_kmers;
}

std::vector<Kmer> mzr::sketch_string(std::string const &s, uint8_t w, uint8_t k, std::unordered_set<uint32_t> const &hfk)
{
    // setup
    std::vector<Kmer> sketch;
    sketch.reserve(s.length() / w);
    KmerGenerator kmer_gen(s, k, 0x3FFFFFFF, true);

    Kmer min_k = kmer_gen.get_kmer();

    // Initialize min_k to first minimizer of the window w
    uint32_t window_pos = 1; // the global window position (window is s[window_pos-w: window_pos)

    while (window_pos < w && !kmer_gen.empty())
    {
        Kmer candidate = kmer_gen.get_kmer();
        if (candidate.seq < min_k.seq)
        {
            min_k = candidate;
        }
        window_pos++;
    }
    // store first minimizer in sketch
    sketch.emplace_back(min_k);

    // construct remaining sketch

    while (!kmer_gen.empty())
    {
        Kmer candidate = kmer_gen.get_kmer(); // move window forward
        window_pos++;
        // upon getting the next kmer, the current minimum kmer will be outside the window.
        // we therefore need to compute the new minimum from scratch. This may be the next kmer
        // or, it may be any of the other w-1 kmers in the window.
        bool lfk = hfk.find(min_k.seq) == hfk.end();
        if (lfk && min_k.pos < window_pos - w)
        {
            min_k = kmer_gen.min_kmer_in_window(w);;                
            sketch.emplace_back(min_k);
        }
        else if (lfk && candidate.seq < min_k.seq)
        {
            min_k = candidate;
            sketch.emplace_back(min_k);
        }
    }

    // DONE: removed unnecessary copying
    // std::vector<Kmer> ret;
    // for (auto const &k : sketch)
    // {
    //     ret.push_back(k);
    // }
    return sketch;
}


std::vector<std::vector<Kmer>> mzr::sketch_contigs(ContigContainerPtr contigs, uint8_t w, uint8_t k, std::vector<std::vector<Kmer>>& sketches)
{

    // count the kmers in the input
    logger.Info("counting kmers...");
    auto kc_size = (app_stats.TotalContigSeqSize / contigs->size() - k + 1) * contigs->size();
    std::unordered_map<uint32_t, uint32_t> kc;
    kc.reserve(kc_size);
    
    for (const ExtContig &c : *contigs)
    {
        count_kmers(c.Seq, kc, k);
    }

    logger.Info("total kmers: " + Util::to_str(kc.size())); 

    // extract high-frequency kmers
    std::unordered_set<uint32_t> hfk = get_high_frequency_kmers(kc);

    logger.Info("total avoided high-frequency kmers " + Util::to_str(hfk.size()));

    logger.Info("Sketching kmers...");

    //std::vector<std::vector<Kmer>> sketches;
    sketches.reserve(contigs->size());
    long total_kmers = 0;
    for (const ExtContig &c : *contigs)
    {
        auto contig_sketches = sketch_string(c.Seq, w, k, hfk);
        total_kmers += contig_sketches.size();
        sketches.push_back(contig_sketches);
    }
    logger.Debug("Sketched " + Util::to_str(total_kmers) + " kmers");

    return (sketches);
}