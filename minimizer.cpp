#include "minimizer.h"

using namespace mzr;

uint8_t mzr::char_to_num(char c, bool complement)
{
    switch (c)
    {
    case 'A':
        return !complement ? 0 : 3;
    case 'C':
        return !complement ? 1 : 2;
    case 'G':
        return !complement ? 2 : 1;
    case 'T':
        return !complement ? 3 : 0;
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

void mzr::count_kmers(std::string const &s, std::unordered_map<unsigned int, unsigned int> &kc, uint8_t k)
{
    KmerGenerator kmer_gen(s, k, 0x3FFFFFFF, true);
    while (!kmer_gen.empty())
    {
        Kmer kmer = kmer_gen.get_kmer();
        kc[kmer.seq]++;
    }
}

std::unordered_set<unsigned int> mzr::get_high_frequency_kmers(std::unordered_map<unsigned int, unsigned int> kc, double percentile)
{
    std::vector<uint64_t> freqs;
    freqs.reserve(kc.size());
    std::unordered_set<unsigned int> high_freq_kmers;
    for (auto k : kc)
    {
        freqs.push_back(k.second);
    }
    std::sort(freqs.begin(), freqs.end());
    size_t p_idx = percentile * freqs.size() + 1;
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
    return high_freq_kmers;
}

std::shared_ptr<std::vector<Kmer>> mzr::sketch_string(std::string const &s, uint8_t w, uint8_t k, 
    std::unordered_set<unsigned int> const &hfk)
{
    // setup
    std::vector<Kmer> sketch;
    sketch.reserve(s.length() / w);
    KmerGenerator kmer_gen(s, k, 0x3FFFFFFF, true);

    Kmer min_k = kmer_gen.get_kmer();

    // Initialize min_k to first minimizer of the window w
    unsigned int window_pos = 1; // the global window position (window is s[window_pos-w: window_pos)

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
        if (lfk && min_k.pos < (window_pos - w))
        {
            auto temp_min_i = kmer_gen.min_kmer_in_window(w);
            min_k = temp_min_i;
            
            sketch.emplace_back(min_k);
        }
        else if (lfk && candidate.seq < min_k.seq)
        {
            min_k = candidate;
            if (min_k.seq == 50051280)
                bool stop = false;
            sketch.emplace_back(min_k);
        }
    }

    // DONE: removed unnecessary copying
    // std::vector<Kmer> ret;
    // for (auto const &k : sketch)
    // {
    //     ret.push_back(k);
    // }

    return std::make_shared<std::vector<Kmer>>(sketch);
}

std::shared_ptr<std::vector<std::vector<Kmer>>> mzr::sketch_contigs(ContigContainerPtr contigs, uint8_t w, uint8_t k, double percentile)
{

    // count the kmers in the input
    logger.Info("counting kmers...");
    std::unordered_map<unsigned int, unsigned int> kc(contigs->size());

    for (ExtContig &c : *contigs)
    {
        count_kmers(*(c.GetSeq()), kc, k);
    }

    logger.Info("total kmers: " + Util::convert_to_string(kc.size())); 

    // extract high-frequency kmers
    std::unordered_set<unsigned int> hfk = get_high_frequency_kmers(kc, percentile);

    logger.Info("total avoided high-frequency kmers " + Util::convert_to_string(hfk.size()));

    logger.Info("Sketching kmers...");

    std::vector<std::vector<Kmer>> sketches;
    for (ExtContig &c : *contigs)
    {
        sketches.push_back(*sketch_string(*c.GetSeq(), w, k, hfk));
    }
    logger.Debug("Sketched " + Util::convert_to_string(sketches.size()) + " kmers");

    return std::make_shared<std::vector<std::vector<Kmer>>>(sketches);
}