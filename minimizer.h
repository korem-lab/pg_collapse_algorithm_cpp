#ifndef MINIMIZER_H
#define MINIMIZER_H

#include "ext_contig.h"
#include "aux_functions.h"
#include "global_copan.h"
#include <assert.h>
#include <string>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <iomanip>
#include <bitset>
#include <iostream>
#include <fstream>

namespace mzr
{
const uint8_t POSITIVE = 0;
const uint8_t NEGATIVE = 1;

uint8_t char_to_num(char c);
uint8_t char_to_num_complement(char c);

char num_to_char(uint8_t b);

/* void my_assert(bool b, std::string m = "")
{
    if (!b)
    {
        logger.Error("condition failed with message: " + m);
        exit(-1);
    }
} */

const char int_to_char_map[] = {'A','C','G','T'};



std::string unpack_kmer(unsigned int packed_kmer, unsigned int kmer_len);


/* inline char num_to_char(uint8_t b)
{
    if (b < 4)
        return int_to_char_map[b];
    return 'N';
} */

struct Kmer
{
    unsigned int seq;
    unsigned int pos;
    uint8_t sign;
    Kmer(unsigned int s, unsigned int p, uint8_t sg) : seq(s), pos(p), sign(sg) {}
    Kmer(const Kmer& k)
    {
        seq = k.seq;
        pos = k.pos;
        sign = k.sign;
    }
    Kmer() : seq(0), pos(0), sign(POSITIVE){};

    std::string as_string(int8_t k) const
    {        
        return unpack_kmer(seq, k);
    }    
};

struct KmerGenerator
{
    KmerGenerator(std::string const &s, uint8_t k, unsigned int mask, bool lex_low = false) : seq(s), k(k), mask(mask), lxl(lex_low)
    {
        seq_it = seq.begin() + k;
        kmer = Kmer(pack(seq.begin(), seq_it), 0, POSITIVE);

        // replace the kmer with the lexicographically lowest version between it
        // and it's reverse complement
        if (lex_low)
        {
            // std::string rev_seq(s);
            // std::reverse(rev_seq.begin(), rev_seq.end());
            // auto rev_seq_it = rev_seq.end() - k -  1;
            // rev_kmer = Kmer(pack(rev_seq_it + 1, rev_seq.end(), true), 0, NEGATIVE);

            rev_kmer = Kmer(pack_reverse(seq.begin(), seq_it, true), 0, NEGATIVE);
        }
        init = true;
    }

    bool empty() const { return seq_it >= seq.end(); }

/*     std::string kmer_as_string()
    {
        return seq.substr(kmer.pos, k);
    } */

    Kmer get_kmer()
    {
        
        if (!init) // it's an intermediate call, and we need to pack a new base
        {            
            unsigned int pos = seq_it - k - seq.begin() + 1;
            uint32_t n = char_to_num(*seq_it);

            kmer = Kmer(mask & (kmer.seq << 2 | n), pos, POSITIVE);            

            if (lxl)
            {
                uint32_t nc = char_to_num_complement(*seq_it);
                rev_kmer = Kmer(mask & (rev_kmer.seq >> 2 | (nc << (k * 2 - 2))), pos, NEGATIVE);
            }
            seq_it++;
            return (rev_kmer.seq < kmer.seq && lxl) ? rev_kmer : kmer;
        }
        else // first call
        {
            init = false;
            return (rev_kmer.seq < kmer.seq && lxl) ? rev_kmer : kmer;
        }
    }

    Kmer min_kmer_in_window(uint8_t w) const
    {
        std::string::const_iterator f = seq_it - w + 1; // move f to the start of the forward window
        unsigned int _kmer = pack(f - k, f);
        unsigned int _rev_kmer = pack_reverse(f - k, f, true);

        // select lowest kmer
        unsigned int min_k = (_kmer < _rev_kmer) ? _kmer : _rev_kmer;
        unsigned int min_pos = f - seq.begin() - k;
        uint8_t sign = (_kmer < _rev_kmer) ? POSITIVE : NEGATIVE;

        // point to the next unprocessed characters
       
        while (f <= seq_it)
        {
            if (_kmer < min_k)
            {
                min_k = _kmer;
                min_pos = f - seq.begin() - k;
                sign = POSITIVE;
            }
            if (_rev_kmer < min_k)
            {
                min_k = _rev_kmer;
                min_pos = f - seq.begin() - k;
                sign = NEGATIVE;
            }            
            
            _kmer = mask & (_kmer << 2 | char_to_num(*f));
            _rev_kmer = mask & (_rev_kmer >> 2 | char_to_num_complement(*f) << (k * 2 - 2));
            f++;
        }
        return {min_k, min_pos, sign};
    }

    unsigned int pack(std::string::const_iterator begin, std::string::const_iterator end, bool complement = false) const
    {
        unsigned int _kmer = 0;
        unsigned int base;
        while (begin < end)
        {
            base = char_to_num(*begin);
            _kmer = (_kmer << 2) | base;
            begin++;
        }
        return mask & (!complement ? _kmer : ~_kmer);
    }
    unsigned int pack_reverse(std::string::const_iterator begin, std::string::const_iterator end, bool complement = false) const
    {
        unsigned int _kmer = 0;
        unsigned int base;
        int i = 0;
        while (begin < end)
        {
            base = char_to_num(*begin);
            _kmer = _kmer | base << 2 * i;
            i++;
            begin++;
        }
        return mask & (!complement ? _kmer : ~_kmer);
    }

    std::string seq;
    std::string::iterator seq_it;
    uint8_t k;
    unsigned int mask;
    bool lxl, init;
    Kmer kmer, rev_kmer;
    
};

void count_kmers(std::string const &s, std::unordered_map<unsigned int, unsigned int> &kc, uint8_t k);

std::unique_ptr<std::unordered_set<unsigned int>> get_high_frequency_kmers(std::unordered_map<unsigned int, unsigned int> kc);

std::unique_ptr<std::vector<Kmer>> sketch_string(std::string const &s, uint8_t w, uint8_t k, std::unordered_set<unsigned int> const &hfk);

std::vector<std::vector<Kmer>> sketch_contigs(ContigContainerPtr contigs, uint8_t w, uint8_t k, std::vector<std::vector<Kmer>>& sketches);

}

#endif