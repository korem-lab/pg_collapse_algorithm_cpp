#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <omp.h>
#include <iostream>
#include <cmath>
#include "ext_contig.h"
#include "global_copan.h"
#include "aux_functions.h"
#include "minimizer.h"
#include <tuple>


const std::string FWD = "forward";
const std::string REF = "reverse";
const uint8_t COMPATIBLE_TYPES= 3;

struct Anchor
{
    unsigned int q_pos;  // query position
    unsigned int s_pos;  // subject position
    unsigned int s_cid;  // subject contig id
    uint8_t type;
};

struct KmerMatch
{
    unsigned int q_pos;
    unsigned int s_pos;
};
struct Alignment
{
    unsigned int q_cid;   // query contig id
    unsigned int s_cid;   // subject contig id
    unsigned int q_begin; // query alignment start position
    unsigned int s_begin; // subject alignment start position
    unsigned int q_end;   // query alignment end position
    unsigned int s_end;   // subject alignment end position
    //string ori;    // orientation
    //float seq_div; // sequence divergence
    std::vector<KmerMatch> km;  // The q_pos ordered kmer coordinates
};
struct PyAlignment
{
    PyAlignment(){}
    PyAlignment(const Alignment& a): q_cid{a.q_cid}, s_cid{a.s_cid}, q_begin{a.q_begin}, 
        s_begin{a.s_begin}, q_end{a.q_end}, s_end{a.s_end}{}
    unsigned int q_cid;
    unsigned int s_cid;
    unsigned int q_begin;
    unsigned int s_begin;
    unsigned int q_end;
    unsigned int s_end;
};

struct Kmer
{
    Kmer(){};
    Kmer(std::string& seq_, unsigned int pos_, uint8_t sign_):kmer{seq_}, pos{pos_}, sign{sign_}{}
    std::string kmer;
    unsigned int pos;
    uint8_t sign;
};

struct Contig
{
    Contig(){}
    Contig(unsigned int cid_, unsigned int len_, std::vector<Kmer>& sketch_, std::string& hdr_):
        cid{cid_}, len{len_}, hdr{hdr_} {sketch = sketch_;}
    unsigned int cid;  // contig id
    unsigned int len;  // length of the contig
    std::vector<Kmer> sketch;
    std::string hdr;
};

struct Coordinate
{
    unsigned int s_pos;
    unsigned int s_cid;
    uint8_t sign;
};

void build_kmer_map(std::unordered_map<std::string, std::vector<Coordinate>> & kmer_map, const std::vector<Contig> & contigs, unsigned int k)
{
    logger.Info("Building kmer map");

    std::vector<Contig>::const_iterator contig_it = contigs.begin();
    std::unordered_map<std::string, std::vector<Coordinate>>::iterator kmap_it;
    unsigned int i, j;
    
    logger.Info("Indexing " + Util::convert_to_string(contigs.size()) +  " contigs");

    while(contig_it != contigs.end())
    {
        for(i = 0; i < contig_it->sketch.size(); i++)
        {
            Coordinate coord;
            coord.s_pos = contig_it->sketch[i].pos;
            coord.s_cid = contig_it->cid;
            coord.sign = contig_it->sketch[i].sign;

            auto it = kmer_map.find(contig_it->sketch[i].kmer);
            if (it == kmer_map.end())
            {
                kmer_map[contig_it->sketch[i].kmer].push_back(coord);
            }
            else
            {
                it->second.push_back(coord);
            }
        }
        contig_it++;
    }

    int count = 0;
    for(auto [kmer, coord] : kmer_map)
    {
        for(auto c : coord)
        {
            count++;
        }
        count ++;
    }
    //kmer_map.size()
    logger.Info("Built kmer map: " + Util::convert_to_string(count) + " entries");

}

bool sort_by_contig_then_pos(const Anchor& one, const Anchor& two)
{
    if (one.s_cid < two.s_cid)
        return true;
    else if (one.s_cid == two.s_cid and one.q_pos < two.q_pos)
        return true;
    else
    {
        // The other options:
        // one.s_cid >= two.s_cid, or one.q_pos >= two.q_pos
        // both of which we say one is bigger than two
        return false;
    }
}

unsigned int abs_uint32_t(unsigned int a, unsigned int b)
{
    return (a>b)*(a-b) + (a<=b)*(b-a);
}

double compute_seq_divergence(const std::vector<Kmer> & sketch, unsigned int start_pos, unsigned int end_pos, unsigned int chain_len, uint32_t k)
{
    /*
    computes maximum likelihood estimate of seq divergence according to minimap2 paper section 2.1.4
    :param first_q_pos: q_pos of first anchor in chain
    :param last_q_pos: q_pos of last anchor in chain
    :param chain_len: number of anchors in chain
    :param k: k-mer size
    :return: sequence divergence
    */
    unsigned int tmp;
    if (start_pos <= end_pos)
        end_pos = end_pos + k;
    else{
        tmp = end_pos;
        end_pos = start_pos + k;
        start_pos = tmp;}

    size_t i = 0;
    double n_seeds = 0;

    while (i < sketch.size())
    {
        if (start_pos <= sketch[i].pos && sketch[i].pos <= end_pos)
            n_seeds += 1;
        i++;
    }
    double match_rate = n_seeds / chain_len;
    //match_rate = min_double(match_rate, 1)
    double seq_divergence = log(match_rate) / k;

    return seq_divergence;
}

int64_t get_max_chain_end(std::vector<double> & scrt)
{
    size_t max_idx = 0;
    size_t it = 1;
    double max_scr = scrt[0];

    while (it < scrt.size())
    {
        if (max_scr < scrt[it])
        {
            max_scr = scrt[it];
            max_idx = it;
        }
        it++;
    }
    return max_idx;
}

std::vector<Alignment> chain_and_backtrack(std::vector<Contig>& container_buf, const Contig& contig, const std::unordered_map<std::string, std::vector<Coordinate>> & kmer_map,
    const std::string& orientation, double divergence_threshold, unsigned int kmer_len, bool asymmentric, double large_gap = 2.0, 
    double small_gap = 0.5, unsigned int max_jump = 200, unsigned int min_overlap = 100)
{
    unsigned int q_len = contig.len;
    std::unordered_map<std::string, std::vector<Coordinate>>::const_iterator kmap_it;
    std::vector<Coordinate>::const_iterator coord_it;
    std::vector<Alignment> alignments;
    std::vector<Anchor> anchors;
    std::vector<Kmer>::const_iterator sketch_it;
    unsigned int q;
    int64_t r;
    int sign_flip;
    unsigned int print_it = 0;

    if(orientation == FWD)
    {
        sketch_it = contig.sketch.cbegin();
        sign_flip = 1;
    }
    else
    {
        sketch_it = contig.sketch.cend() - 1;
        sign_flip = -1;
    }
    
    for(int r = 0; r < contig.sketch.size(); r++)
    {
        auto flip_sketch_it = sketch_it + sign_flip*r;
        kmap_it = kmer_map.find(flip_sketch_it->kmer);
        if (kmap_it == kmer_map.cend())
            logger.Error("Failed to find a k-mer that should definitely be in the hash table");
        coord_it = kmap_it->second.begin();
        while(coord_it != kmap_it->second.cend())
        {
            if (contig.cid == coord_it->s_cid && flip_sketch_it->pos >= coord_it->s_pos)
            {
                coord_it++;
                continue;
            }
            if (asymmentric && (contig.cid > coord_it->s_cid))
            {
                coord_it++;
                continue;
            }
            Anchor anchor;            
            anchor.q_pos = flip_sketch_it->pos;
            anchor.s_pos = coord_it->s_pos;
            anchor.s_cid = coord_it->s_cid;
            anchor.type = 2 * flip_sketch_it->sign + coord_it->sign;
            anchors.push_back(anchor);
            coord_it++;        
        }
    }

    std::sort(anchors.begin(), anchors.end(), sort_by_contig_then_pos);
    std::vector<Anchor>::const_iterator s_left_ptr, s_right_ptr;
    
    std::vector<KmerMatch> kmer_matches;

    unsigned int min_q, max_q, min_s, max_s;
    unsigned int q_next, s_next, q_prev, s_prev;
    unsigned int jump_div, s_jump, q_jump;
    unsigned int chain_start, chain_len;
    int64_t pos, new_pos, mcs, mce, tmp;
    int64_t EXHAUST = -2;
    int64_t NEW_ALN = -1;
    unsigned int last_match, first_match;
    unsigned int max_id;
    uint8_t type_next, type_prev;
    double alpha;
    double gap_cost;
    double next_score;
    double chain_max_score;
    double max_score;
    float seq_div_q, seq_div_s, seq_div;

    s_right_ptr = anchors.cbegin();
    while(s_right_ptr < anchors.cend())
    {
        //Set left and right pointing to the same entry
        s_left_ptr = s_right_ptr;
        // Looping through until we reach the first anchor of a different contig
        while(s_right_ptr < anchors.cend() && s_left_ptr->s_cid == s_right_ptr->s_cid)
        {
            s_right_ptr++;
        }
        std::vector<Anchor> s_anchors;
        // Copy over anchors of subject contig
        s_anchors.insert(s_anchors.cbegin(), s_left_ptr, s_right_ptr);

        // Optimization:
        // This if statement is true when the distance between the most extreme anchors
        // is less than min_overlap. In this case, its not possible for an alignment to be >= min_overlap
        // so continue.

        min_q = s_anchors[0].q_pos;
        max_q = s_anchors[s_anchors.size() - 1].q_pos;
        min_s = s_anchors[0].s_pos;
        max_s = min_s;

        for(int t = 0; t < s_anchors.size(); t++)
        {
            min_s = std::min(min_s, s_anchors[t].s_pos);
            max_s = std::max(max_s, s_anchors[t].s_pos);
        }
        if (((max_q - min_q) < min_overlap) || ((max_s - min_s) < min_overlap))
            continue;

        // ================================================
        // ------ Compute score and backtrack table -------
        // ================================================
        std::vector<double> score_table(s_anchors.size(), 0.0);
        std::vector<int64_t> backtrack_table(s_anchors.size(), -1);

        for(int i = 1; i < s_anchors.size(); i++)
        {
            max_score = 0;
            max_id = 0;
            q_next = s_anchors[i].q_pos;
            s_next = s_anchors[i].s_pos;
            type_next = s_anchors[i].type;
            // Iterate backwards from i to 0
            for(int j = i - 1; j > -1; j--)
            {
                q_prev = s_anchors[j].q_pos;
                s_prev = s_anchors[j].s_pos;
                type_prev = s_anchors[j].type;
                
                // If true, the distances are small enough to be considered part of the same alignment
                if ((0 < (q_next - q_prev) && (q_next - q_prev) < max_jump) &&
                    (0 < abs_uint32_t(s_next, s_prev) && abs_uint32_t(s_next, s_prev) < max_jump) &&
                    (type_prev == type_next || type_prev + type_next == COMPATIBLE_TYPES))
                {
                    alpha = std::min(std::min(q_next - q_prev, abs_uint32_t(s_next, s_prev)), kmer_len);
                    auto q_jump = q_next - q_prev;
                    auto s_jump = abs_uint32_t(s_next, s_prev);
                    auto jump_div = abs_uint32_t(s_jump, q_jump);
                    gap_cost = jump_div > 100 ? (large_gap * jump_div) : (small_gap * jump_div);
                    next_score = score_table[j] + alpha - gap_cost;
                    
                    if (next_score > max_score)
                    {
                        max_score = next_score;
                        max_id = j;
                        // This is a heuristic that reduced the extent we iterate back in the inner forloop.
                        // in the worst case, we iterate back to 0, in which case O(N^2). However,
                        // The assumption is made that:
                        // 1) if the alignments are perfectly aligned: jump_div == 0
                        // 2) The alignments are identical across the k previous characters: (q_next - q_prev < k),
                        //   means that there was no sequence difference within the previous k characters. If there
                        //   was, their could be no exact anchor match, the closest possible next anchor would be
                        //   > k away.
                        // Then, this previous anchor is probably the optimal choice. Therefore, don't bother
                        // Iterating back to guarantee the optimal choice. This heuristic takes Heng Li's
                        // heuristic (look back h places) to the extreme, h = 1. So, alignments are computed in
                        // O(N) time with high probability.
                        if (jump_div == 0 && q_next - q_prev < kmer_len)
                            break;

                    }
                }
                // If the next anchor is so far away it's considered part of a separate alignment,
                // (q_next - q_prev > max_jump) then don't consider it.
                // TODO PROPOSAL, should we check both query and subject?
                if (q_next - q_prev > max_jump)// or (q_next != q_prev and abs_uint32_t(s_next, s_prev) > max_jump):
                    // if we break the outer loop, we're done aligning, so set the alignment to neutral
                    break;
            }
            // UPDATE SCORE AND BACKTRACK

            score_table[i] = std::max(max_score, static_cast<double>(kmer_len));
            if (max_score > kmer_len)
                backtrack_table[i] = max_id;
        }   

        // ================================================
        // -------- Backtrack to get alignments -----------
        // ================================================

        mce = get_max_chain_end(score_table);

        while (score_table[mce] > min_overlap)
        {
            mcs = mce;
            chain_max_score = score_table[mcs];
            chain_len = 1;
            while (backtrack_table[mcs] >= 0)
            {
                tmp = mcs;
                if (score_table[mcs] > chain_max_score)
                {
                    logger.Info("Encountered edge case!\n");
                }

                mcs = backtrack_table[mcs];
                backtrack_table[tmp] = EXHAUST;
                score_table[tmp] = EXHAUST;
                chain_len += 1;

                // We dont need all the kmer matches, so, a reasonable minimzation is to take the kmers
                // a distance of at least k away
                //match_pos.q_pos = s_anchors[tmp].q_pos
                //match_pos.s_pos = s_anchors[tmp].s_pos
                //if kmer_matches.empty() or kmer_matches.back().q_pos - match_pos.q_pos > k*2:
                //    kmer_matches.push_back(match_pos)                
            }
            if (backtrack_table[mcs] == EXHAUST)
            {
                // then we must be backtracking along a chain that's already part of a chain with a higher score
                mce = get_max_chain_end(score_table);
                continue;
            }
            // else, we are at the start of a new alignment, i.e backtrack_table[msc] == NEW_ALN
            // so process the alignment

            // first, exhaust this anchor
            backtrack_table[mcs] = EXHAUST;
            score_table[mcs] = EXHAUST;

            seq_div_q = compute_seq_divergence(contig.sketch, s_anchors[mcs].q_pos, s_anchors[mce].q_pos, chain_len, kmer_len);
            seq_div_s = compute_seq_divergence(container_buf[s_anchors[0].s_cid].sketch, s_anchors[mcs].s_pos, s_anchors[mce].s_pos, chain_len, kmer_len);
            // next, store the alignment
            seq_div = (seq_div_q + seq_div_s) / 2;
  
            //printf(
            //    "{\"q_cid\":%d, \"s_cid\":%d, \"q_beg\":%d, \"q_end\":%d, \"s_beg\":%d, \"s_end\":%d, \"sd\":%f, \"ori\":\"%s\", \"jd\":%d}\n",
            //    contig.cid, s_anchors[0].s_cid, s_anchors[mcs].q_pos, s_anchors[mce].q_pos + kmer_len,
            //    s_anchors[mcs].s_pos, s_anchors[mce].s_pos + kmer_len, seq_div, orientation.c_str(),
            //    abs_uint32_t(abs_uint32_t(s_anchors[mcs].s_pos, s_anchors[mce].s_pos),
            //                 abs_uint32_t(s_anchors[mcs].q_pos, s_anchors[mce].q_pos))
            //)
            //jump_div = abs_uint32_t(abs_uint32_t(s_anchors[mcs].s_pos, s_anchors[mce].s_pos),
            //            abs_uint32_t(s_anchors[mcs].q_pos, s_anchors[mce].q_pos))

            if (abs_uint32_t(s_anchors[mce].q_pos + kmer_len, s_anchors[mcs].q_pos) >= min_overlap && seq_div <= divergence_threshold)
            {

                Alignment alignment;
                alignment.q_cid = contig.cid;
                alignment.s_cid = s_anchors[0].s_cid;
                alignment.q_begin = s_anchors[mcs].q_pos;
                alignment.q_end = s_anchors[mce].q_pos + kmer_len;
                alignment.s_begin = s_anchors[mcs].s_pos;
                alignment.s_end = s_anchors[mce].s_pos + kmer_len;
                    //if s_anchors[mcs].s_pos < s_anchors[mce].s_pos + kmer_len:
                    //    alignment.s_begin = s_anchors[mcs].s_pos
                    //    alignment.s_end = s_anchors[mce].s_pos + kmer_len
                    //else:
                    //    alignment.s_begin = s_anchors[mce].s_pos
                    //    alignment.s_end = s_anchors[mcs].s_pos + kmer_len


                    //alignment.ori = orientation
                    //alignment.seq_div = seq_div
                    //alignment.km = kmer_matches
                alignments.push_back(alignment);

            }
                // find the next best alignment
            mce = get_max_chain_end(score_table);
        }
    }

    return alignments;
}


std::vector<PyAlignment> align_contigs(ContigContainerPtr container, int k, double divergence_threshold, int window_sz=1, int num_threads=10, bool asymmetric = false,
    double lg = 2.0, double sg=0.5, int mj = 200, int mo=100, double hfkf = 1e-05)
{
    size_t N = container->size();
    std::vector<Contig> container_buf;
    size_t idx, i;
    size_t nthreads =  num_threads;
    unsigned int kmer_len = k;
    std::unordered_map<std::string, std::vector<Coordinate>> kmer_map;
    size_t print_chunk = std::ceil(N/10);
    double high_freq_kmer_filt = hfkf;

    logger.Info("Sketching " + Util::convert_to_string(N) + " contigs.");

    auto sketches = * mzr::sketch_contigs(container, window_sz, kmer_len, 1-high_freq_kmer_filt);
    logger.Info("Completed sketch_contigs: " + Util::convert_to_string(sketches.size()));

    for(int i = 0; i < N; i++)
    {
        if (i % 10000 == 0)
            logger.Debug("Transferring " + Util::convert_to_string(i) + " of sketches... ");

        std::vector<Kmer> _s;
    
        for (auto _k : sketches[i])
        {
            _s.push_back(Kmer(*_k.as_string(kmer_len), _k.pos, _k.sign));
        }

        ExtContig ex = container->at(i);
        Contig c(ex.GetId(), ex.GetSeq()->size(), _s, *ex.GetHeader());
        container_buf.push_back(c);
    }     

    build_kmer_map(kmer_map, container_buf, kmer_len);

    logger.Info("Launching chain_and_backtrack. Multithreaded? " + Util::convert_to_string(num_threads > 0 ? "true" : "false"));

    std::vector<Kmer> kmer_iter_fwd, kmer_iter_rev;
    std::vector<std::vector<Alignment>> thrd_aln(N);
    bool asym = asymmetric;
    double large_gap = lg;
    double small_gap = sg;
    unsigned int max_jump = mj;
    unsigned int min_overlap = mo;

    #pragma omp parallel shared(thrd_aln, container_buf, kmer_map)
    {
        #pragma omp for schedule(dynamic)
        for(int i = 0; i < N; i++)
        {
            #pragma omg task
            {
                if (i % 5000 == 0)
                    logger.Debug("analyzed contig " + Util::convert_to_string(i));
            }
            #pragma omp task
            {                
                thrd_aln[i] = (chain_and_backtrack(container_buf, container_buf[i], kmer_map, FWD, divergence_threshold, kmer_len, asym, large_gap, 
                          small_gap, max_jump, min_overlap));
            }
        }    
    }
    kmer_map.clear();

    logger.Info("Merging per-thread alignments");
    
    std::vector<PyAlignment> pyalignments;
    
    std::vector<std::vector<Alignment>>::const_iterator aln_contig_it = thrd_aln.begin();
    std::vector<Alignment>::const_iterator aln_it;
    std::vector<KmerMatch>::const_iterator km_it;

    while(aln_contig_it != thrd_aln.end())
    {
        aln_it = aln_contig_it->begin();
        while(aln_it != aln_contig_it->end())
        {
            pyalignments.push_back(PyAlignment(*aln_it));
            aln_it++;
        }
        aln_contig_it++;
    }
   
    logger.Info("Returning from Align");
    logger.Info("Num py_alignments: " + Util::convert_to_string(pyalignments.size()));   
    return pyalignments;
}


std::vector<PyAlignment> compute_alignments(ContigContainerPtr container)
{
    auto k = config.GetValue<int>("kmer_size");
    auto divergence_threshold = config.GetValue<double>("divergence_threshold");
    auto sensitive_mode = config.GetValue<bool>("sensitive_mode");
    divergence_threshold -= sensitive_mode ? 0.02 : 0;
    divergence_threshold = std::max(0.05, divergence_threshold);
    auto window_sz = config.GetValue<int>("window_size");
    auto num_of_threads = config.GetValue<int>("num_threads");
    return align_contigs(container, k, divergence_threshold, window_sz, num_of_threads);
}

#endif