#ifndef COPAN_GRAPH_H
#define COPAN_GRAPH_H
#include <string>
#include <vector>
#include <unordered_map>
#include <exception>
#include <set>
#include <iostream>
#include <ostream>
#include <functional>
#include <memory>
#include "ext_contig.h"
#include "gluepoints.h"
#include "pe_ext_reader.h"
#include "ext_contig.h"

template <>
struct std::hash<std::pair<unsigned int, unsigned int>>
{
    std::size_t operator()(const std::pair<const unsigned int, const unsigned int> &v) const
    {
        size_t h = (size_t(v.first) << 32) + size_t(v.second);
        h *= 1231231557ull; // "random" uneven integer.
        h ^= (h >> 32);
        return h;
    }
};

std::pair<unsigned int, unsigned int> get_key(unsigned int x, unsigned int y)
{
    if (x < y)
        return std::make_pair(x, y);
    return std::make_pair(y, x);
}

struct GFAType
{
    static const char S = 'S';
    static const char L = 'L';
    static const char H = 'H';
};

class Element
{
public:
    virtual std::string ToString() { return "Not implemented"; }
    virtual std::string ToToken() { return "Not implemented"; }
};

typedef std::shared_ptr<Element> pElement;

class GFALine : public Element
{
public:
    GFALine(const std::string &l_nid, char l_ori, const std::string &r_nid, char r_ori, const std::string &cigar) : m_l_nid{l_nid}, m_l_ori{l_ori}, m_r_nid{r_nid}, m_r_ori{r_ori}, m_cigar{cigar}
    {
    }
    std::string ToString() override
    {
        return "L\t" + m_l_nid + "\t" + m_l_ori + "\t" + m_r_nid + "\t" + m_r_ori + "\t" + m_cigar;
    }
    std::string ToToken() override
    {
        return ToString();
    }

private:
    const std::string m_l_nid;
    char m_l_ori;
    const std::string m_r_nid;
    char m_r_ori;
    const std::string m_cigar;
    static const char gfa_type = GFAType::L;
};

class Rest : public Element
{
public:
    static const char gfa_type = GFAType::S;
    std::string nid;
    std::string full_name;
    int sid;
    std::string contig_name;
    std::string start;
    std::string end;
    char ori;
    bool is_tag = false;
    std::string ToString() override { return "\tINFO:" + nid + ":" + full_name + ":" + Util::convert_to_string(sid) + ":" + contig_name + 
        ":" + start + ":" + end + ":" + ori + ":" + Util::convert_to_string(is_tag); }
    std::string ToToken() override { return nid + ":" + Util::convert_to_string(sid) + ":" + contig_name + ":" + 
        start + ":" + end + ":" + ori + ":" +  Util::convert_to_string(is_tag); }
};
class GFASegment : public Element
{
public:
    GFASegment(const char gfa_type, std::string &nid, std::string &seq, Rest &rest) : m_gfa_type{gfa_type}, m_nid{nid}, m_seq{seq}, m_rest{rest} {}

    std::string ToString() override
    {
        return "S\t" + m_nid + "\t" + m_seq + m_rest.ToString();
    }
    //INFO:{nid}:{full_name}:{sid}:{extended_contig.contig_name}:{start}:{end}:{ori}:{is_tag}
    //   _, nid, fn, sid, cid, lp, rp, ori, tag
    //   f'{nid}:{sid}:{cid}:{lp}:{rp}:{ori}'
    std::string ToToken() override
    {
        return m_rest.nid + ":" + Util::convert_to_string(m_rest.sid) + ":" + m_rest.contig_name + ":" + (m_rest.start) + ":" +
               (m_rest.end) + ":" + m_rest.ori;
    }
    std::string GetLine()
    {
        return "GFASegment(nid=" + m_nid + ", seq=" + m_seq.substr(0, 10) + ")";
    }
    Rest GetRest()
    {
        return m_rest;
    }

private:
    const char m_gfa_type;
    std::string m_nid;
    std::string m_seq;
    Rest m_rest;
};

class CopanNodeSeq : public Element
{
public:
    CopanNodeSeq(std::string &nid, const std::string &full_name, unsigned int sid, unsigned int cid, const std::string &contig_name,
                 unsigned int lp, unsigned int rp, char ori, bool is_tag, const std::string &seq, int fold = 80) : m_nid{nid}, m_full_name{full_name},
                                                                                                                   m_cid{cid}, m_sid{sid}, m_contig_name{contig_name}, m_lp{lp}, m_rp{rp}, m_ori{ori}, m_is_tag{is_tag}, m_seq{seq}, m_fold{fold} {}

    std::string ToString() override
    {
        auto hdr = ">" + m_nid + ":" + m_full_name + ":" + Util::convert_to_string(m_sid) + ":" + Util::convert_to_string(m_cid) +
                   ":" + m_contig_name + ":" + Util::convert_to_string(m_lp) + ":" + Util::convert_to_string(m_rp) + ":" + m_ori + ":" + Util::convert_to_string(m_is_tag);
        hdr += "\n" + m_seq;
        return hdr;
    }
    std::string ToToken() override
    {
        return m_nid + ":" + Util::convert_to_string(m_sid) + ":" + Util::convert_to_string(m_cid) +
               Util::convert_to_string(m_lp) + ":" + Util::convert_to_string(m_rp) + ":" + m_ori;
    }

private:
    std::string m_nid;
    std::string m_full_name;
    unsigned int m_sid;
    unsigned int m_cid;
    std::string m_contig_name;
    unsigned int m_lp;
    unsigned int m_rp;
    char m_ori;
    std::string m_seq;
    bool m_is_tag;
    int m_fold;
};

class Node
{
public:
    Node() {}
    Node(int i, unsigned int alpha, unsigned int beta, unsigned int gamma, unsigned int delta) : uid{++uid_counter}, intervals{i}, init_alpha{alpha}, init_beta{beta}, init_gamma{gamma}, init_delta{delta}, oris{i, '+'}
    {
        name = Util::convert_to_string(uid);
        common_numerals.insert(get_key(init_alpha, init_beta));
        common_numerals.insert(get_key(init_gamma, init_delta));
    }
    Node(const Node& n)
    {
        common_numerals = n.common_numerals;
        name = n.name;
        intervals = n.intervals;
        init_alpha = n.init_alpha;
        init_beta = n.init_beta;
        init_delta = n.init_delta;
        init_gamma = n.init_gamma;
        uid = n.uid;
        bool sorted = n.sorted;
        oris = n.oris;
    }
    std::string GetFullName()
    {
        std::string fullName = "[";
        for (auto p : common_numerals)
        {
            fullName.append("(");
            fullName.append(Util::convert_to_string(p.first));
            fullName.append(",");
            fullName.append(Util::convert_to_string(p.second));
            fullName.append(")");
        }
        fullName.append("]");
        return fullName;
    }

    std::string GetName()
    {
        return name;
    }

    std::vector<int> GetIntervalIndices()
    {
        if (!sorted)
        {
            std::sort(intervals.begin(), intervals.end());
            sorted = true;
        }
        return intervals;
    }
    std::string ToString()
    {
        return "Node(a=" + Util::convert_to_string(init_alpha) + ",b=" + Util::convert_to_string(init_beta) +
               ",g=" + Util::convert_to_string(init_gamma) + ",d=" + Util::convert_to_string(init_delta) + ", name=" +
               GetName() + ",ivl=" + Util::convert_to_csv(intervals) + ",oris=" + std::string{oris.second} + ")";
    }
    void add(int i, unsigned int alpha, unsigned int beta, unsigned int gamma, unsigned int delta,
             std::unordered_map<std::pair<unsigned int, unsigned int>, Node *> &lookup)
    {
        intervals.push_back(i);

        auto key_a_b = get_key(alpha, beta);
        auto key_g_d = get_key(gamma, delta);

        std::set<std::pair<unsigned int, unsigned int>> intersection{};
        std::set<std::pair<unsigned int, unsigned int>> set2{};
        set2.insert(key_a_b);
        set2.insert(key_g_d);

        std::set_intersection(common_numerals.begin(), common_numerals.end(), set2.begin(), set2.end(),
                              std::inserter(intersection, intersection.begin()));

        common_numerals = intersection;

        auto key_a_b_init = get_key(init_alpha, init_beta);
        if (common_numerals.find(key_a_b_init) == common_numerals.end())
            lookup[key_a_b_init] = nullptr;

        auto key_g_d_init = get_key(init_gamma, init_delta);
        if (common_numerals.find(key_g_d_init) == common_numerals.end())
            lookup[key_g_d_init] = nullptr;
        
        if ((alpha == init_alpha && beta == init_beta) || (alpha == init_gamma && beta == init_delta) ||
            (gamma == init_alpha && delta == init_beta) || (gamma == init_gamma && delta == init_delta))
            oris = {i, '+'};
        else if ((alpha == init_beta && beta == init_alpha) || (alpha == init_delta && beta == init_gamma) ||
                 (gamma == init_beta && delta == init_alpha) || (gamma == init_delta && delta == init_gamma))
            oris = {i, '-'};
        else
            throw "No relation between interval numerals and node numerals.";

        if (name == "9560")
            std::cout << name << " -> (" << oris.first << "," << oris.second << ") " << alpha << " " << beta << " " << gamma << " " << delta << "\n";
    }

private:
    friend class CoPanGraph;
    std::string name;
    std::vector<int> intervals;
    std::pair<int, char> oris;
    unsigned int init_alpha, init_beta, init_gamma, init_delta;
    int uid;
    std::set<std::pair<unsigned int, unsigned int>> common_numerals;
    bool sorted = false;
    static int uid_counter;
};

int Node::uid_counter = 0;

class CoPanGraph
{
private:
    std::vector<SequenceInterval> intervals;
    ContigContainerPtr contig_container;
    IdMapPtr cid_to_sid;
    const unsigned num_samples;
    std::vector<Node*> nodes{};
    std::vector<Node*> ivl_to_node{};
    std::unordered_map<std::pair<unsigned int, unsigned int>, Node *> node_lookup{};

public:
    CoPanGraph(std::vector<SequenceInterval> &seq_intervals, ContigContainerPtr contig_container_ptr, IdMapPtr contig_to_sample_map, unsigned int num_of_samples) : intervals{seq_intervals}, contig_container{contig_container_ptr}, cid_to_sid{contig_to_sample_map}, num_samples{num_of_samples}
    {
    }
    unsigned int get_start_of(int i)
    {
        return intervals[i].start;
    }
    unsigned int get_end_of(int i)
    {
        return intervals[i].end;
    }
    char get_orientation_of(int i)
    {
        return ivl_to_node[i]->oris.second;
    }
    int get_sid_of(int i)
    {
        return cid_to_sid->at(intervals[i].cid);
    }
    int get_cid_of(int i)
    {
        return intervals[i].cid;
    }

    void build_graph()
    {
        for (int i = 0; i < intervals.size(); i++)
        {
            auto ivl = intervals[i];
            auto key_a_b = get_key(ivl.alpha, ivl.beta);
            auto key_g_d = get_key(ivl.gamma, ivl.delta);

            Node *node = nullptr;
            if (node_lookup.find(key_a_b) != node_lookup.end())
            {
                node = node_lookup[key_a_b];
            }
            if (!node && node_lookup.find(key_g_d) != node_lookup.end())
            {
                node = node_lookup[key_g_d];
            }

            if (node)
            {
                node->add(i, ivl.alpha, ivl.beta, ivl.gamma, ivl.delta, node_lookup);
                ivl_to_node.push_back(new Node(*node));
            }
            else
            {
                node = new Node(i, ivl.alpha, ivl.beta, ivl.gamma, ivl.delta);
                ivl_to_node.push_back(new Node(*node));

                nodes.push_back(node);

                node_lookup[key_a_b] = node;
                node_lookup[key_g_d] = node;
            }
        }
        logger.Debug("Node_lookup size: " + Util::convert_to_string(node_lookup.size()));
        logger.Debug("ivl_to_node size: " + Util::convert_to_string(ivl_to_node.size()));
    }

    void output_graph(const std::string &gfa_name, const std::string &fasta_name, const std::string &node_table_name, const std::string &edge_table_name)
    {
        write_node_table(node_table_name);
        write_edge_table(edge_table_name);
        write_gfa(gfa_name);
        write_fasta(fasta_name);
    }

    void write_node_table(const std::string &node_table_name)
    {
        std::ofstream node_file(node_table_name, std::ios::out);

        if (node_file.is_open())
        {
            node_file << ",";
            for (int i = 0; i < num_samples; i++)
                node_file << "sample_" << i;
            node_file << "\n";

            for (int i = 1; i <= nodes.size(); i++)
            {
                node_file << i << ",";
                for (int j = 1; j <= num_samples; j++)
                    node_file << j << (j < num_samples ? "," : "");
                node_file << "\n";
            }
            node_file.close();
        }
    }
    void write_edge_table(const std::string &edge_table_name)
    {
        std::set<std::string> table_index{};
        for (int i = 1; i < intervals.size(); i++)
        {
            if (intervals[i - 1].cid != intervals[i].cid)
                continue;
            auto u = ivl_to_node[i - 1]->name;
            auto v = ivl_to_node[i]->name;
            auto ori_u = get_orientation_of(i - 1);
            auto ori_v = get_orientation_of(i);
            auto line = u + "(" + ori_u + ") -> " + v + "(" + ori_v + ")";
                        
            table_index.insert(line);
        }
        // std::sort(table_index.begin(), table_index.end()); // the set is already sorted?
        std::ofstream edge_file(edge_table_name, std::ios::out);
        if (edge_file.is_open())
        {
            edge_file << " ";
            for(int j = 0; j < num_samples; j++)
                edge_file << "sample_" << j << ",";
            edge_file << "\n";
            for (auto row : table_index)
            {
                edge_file << row;
                for (int col = 1; col <= num_samples; col++)
                {
                    edge_file << "," << col;
                }
                edge_file << "\n";
            }
            edge_file.close();
        }
    }
    void write_fasta(const std::string &fasta_name, bool remove_identical = true)
    {
        std::vector<pElement> fastas{};
        for (auto n : nodes)
        {
            for (auto i : n->GetIntervalIndices())
            {
                auto e = build_elem(n->name, get_cid_of(i), get_start_of(i), get_end_of(i), n->GetFullName(),
                                    get_orientation_of(i), get_sid_of(i), "fasta");
                fastas.push_back(e);
            }
        }
        std::set<std::string> tokens{};
        for (auto &e : fastas)
        {
            tokens.insert(e->ToToken());
        }
        std::ofstream fasta_file(fasta_name, std::ios::out);
        if (fasta_file.is_open())
        {
            for (auto &e : fastas)
            {
                auto tk = e->ToToken();
                if (tokens.find(tk) != tokens.end())
                {
                    fasta_file << e->ToString() << "\n";
                    tokens.erase(tk);
                }
                else if (!remove_identical)
                    fasta_file << e->ToString() << "\n";
            }
            fasta_file.close();
        }
    }
    void write_gfa(const std::string &gfa_name, bool remove_identical = true)
    {
        std::vector<pElement> seg_type_list{};
        std::cout << "nodes.size(): " << nodes.size() << "\n";
        for (auto &n : nodes)
        {
            auto intervals = n->GetIntervalIndices();
            for (int i : intervals)
            {
                auto e = build_elem(n->name, get_cid_of(i), get_start_of(i), get_end_of(i),
                                    n->GetFullName(), get_orientation_of(i), get_sid_of(i), "gfa");
                seg_type_list.push_back(e);
            }
        }
        
        std::set<std::string> tokens{};
        for (auto &e : seg_type_list)
            tokens.insert(e->ToToken());
        std::cout << "seg-type-list.size(): " << seg_type_list.size() << std::endl;
        std::vector<pElement> out_list{};

        for (auto &e : seg_type_list)
        {
            auto tk = e->ToToken();
            if (tokens.find(tk) != tokens.end())
            {
                out_list.push_back(e);
                tokens.erase(tk);                
            }
            else if (!remove_identical)
            {
                out_list.push_back(e);
            }
        }
        
        for (int i = 1; i < intervals.size(); i++)
        {
            if (intervals[i - 1].cid != intervals[i].cid)
                continue;
            auto &u = ivl_to_node[i - 1];
            auto &v = ivl_to_node[i];
            auto ori_u = get_orientation_of(i - 1);
            auto ori_v = get_orientation_of(i);
            GFALine gfa_link(u->GetName(), ori_u, v->name, ori_v, std::string{"0"});
            out_list.push_back(std::make_shared<GFALine>(gfa_link));
        }
       
        std::ofstream gfa_file(gfa_name, std::ios::out);
        if (gfa_file.is_open())
        {
            for (auto e : out_list)
            {
                gfa_file << e->ToString() << "\n";
            }

            gfa_file.close();
        }
    }

    std::shared_ptr<Element> build_elem(std::string &nid, int cid, unsigned int start, unsigned int end, const std::string &full_name,
                                        char ori, int sid, const std::string &format)
    {
        auto &extended_contig = contig_container->at(cid);
        auto start_of_contig = extended_contig.GetBwdExt();
        auto end_of_contig = extended_contig.GetBwdExt() + extended_contig.GetUnExtLen() + 1;
        bool is_tag = extended_contig.GetIsTag();
        std::string seq = extended_contig.GetSeq()->substr(start, end-start);
        std::string s_start, s_end;

        if (end < start_of_contig)
        {
            s_start = "0-" + Util::convert_to_string(abs(start - start_of_contig));
            s_end = "0-" + Util::convert_to_string(abs(end - start_of_contig));
            is_tag = true;
        }
        else if (start < start_of_contig && start_of_contig <= end && end < end_of_contig)
        {
            s_start = "0-" + Util::convert_to_string(abs(start - start_of_contig));
            s_end = Util::convert_to_string(end - start_of_contig);
            is_tag = true;
        }
        else if (start_of_contig <= start && start < end_of_contig &&
                 start_of_contig <= end && end < end_of_contig)
        {
            s_start = Util::convert_to_string(start - start_of_contig);
            s_end = Util::convert_to_string(end - start_of_contig);
        }
        else if (start_of_contig <= start && start < end_of_contig && end >= end_of_contig)
        {
            s_start = Util::convert_to_string(start - start_of_contig);
            s_end = Util::convert_to_string(extended_contig.GetUnExtLen()) + "+" +
                    Util::convert_to_string(end - end_of_contig);
            is_tag = true;
        }
        else if (start >= end_of_contig)
        {
            s_start = Util::convert_to_string(extended_contig.GetUnExtLen()) + "+" +
                      Util::convert_to_string(start - end_of_contig);
            s_end = Util::convert_to_string(extended_contig.GetUnExtLen()) + "+" +
                    Util::convert_to_string(end - end_of_contig);
            is_tag = true;
        }
        else if (start < start_of_contig && end >= end_of_contig)
        {
            s_start = "0-" + Util::convert_to_string(start - start_of_contig);
            s_end = Util::convert_to_string(extended_contig.GetUnExtLen()) +
                    "+" + Util::convert_to_string(end - end_of_contig);
            is_tag = true;
        }
        else
            throw("Edge case escaped.");

        Element *e;

        if (format == "fasta")
        {
            auto cns = CopanNodeSeq(nid, full_name, sid, cid, *extended_contig.GetName(), start, end, ori, is_tag, seq);
            return std::make_shared<CopanNodeSeq>(cns);
        }
        else if (format == "gfa")
        {
            Rest rest;
            rest.contig_name = *extended_contig.GetName();
            rest.end = s_end;
            rest.full_name = full_name;
            rest.is_tag = is_tag;
            rest.nid = nid;
            rest.ori = ori;
            rest.sid = sid;
            rest.start = s_start;

            auto gfas = GFASegment(GFAType::S, nid, seq, rest);
            return std::make_shared<GFASegment>(gfas);
        }
        else
            throw "No current format othen than fasta or gfa";
    }
};

#endif