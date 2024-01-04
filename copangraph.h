#ifndef COPAN_GRAPH_H
#define COPAN_GRAPH_H
#include "aux_functions.h"
#include <string>
#include <vector>
#include <map>
#include <exception>
#include <set>
#include <iostream>
#include <ostream>
#include <functional>
#include <memory>
#include "ext_contig.h"
#include "gluepoints.h"
#include "pe_ext_reader.h"


template <>
struct std::hash<std::pair<uint32_t, uint32_t>>
{
    std::size_t operator()(const std::pair<const uint32_t, const uint32_t> &v) const
    {
        size_t h = (size_t(v.first) << 32) + size_t(v.second);
        h *= 1231231557ull; // "random" uneven integer.
        h ^= (h >> 32);
        return h;
    }
};

std::pair<uint32_t, uint32_t> get_key(uint32_t x, uint32_t y)
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
    std::string ToString() override { return "\tINFO:" + nid + ":" + full_name + ":" + Util::to_str<int>(sid) + ":" + contig_name + 
        ":" + start + ":" + end + ":" + ori + ":" + Util::to_str(is_tag); }
    std::string ToToken() override { return nid + ":" + Util::to_str<int>(sid) + ":" + contig_name + ":" + 
        start + ":" + end + ":" + ori + ":" +  Util::to_str(is_tag); }
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
        return m_rest.nid + ":" + Util::to_str<int>(m_rest.sid) + ":" + m_rest.contig_name + ":" + (m_rest.start) + ":" +
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
    CopanNodeSeq(std::string &nid, const std::string &full_name, uint32_t sid, uint32_t cid, const std::string &contig_name,
                 uint32_t lp, uint32_t rp, char ori, bool is_tag, const std::string &seq, int fold = 80) : m_nid{nid}, m_full_name{full_name},
                                                                                                                   m_cid{cid}, m_sid{sid}, m_contig_name{contig_name}, m_lp{lp}, m_rp{rp}, m_ori{ori}, m_is_tag{is_tag}, m_seq{seq}, m_fold{fold} {}

    std::string ToString() override
    {
        auto hdr = ">" + m_nid + ":" + m_full_name + ":" + Util::to_str<uint32_t>(m_sid) + ":" + Util::to_str<uint32_t>(m_cid) +
                   ":" + m_contig_name + ":" + Util::to_str<uint32_t>(m_lp) + ":" + Util::to_str<uint32_t>(m_rp) + ":" + m_ori + ":" + Util::to_str(m_is_tag);
        hdr += "\n" + m_seq;
        return hdr;
    }
    std::string ToToken() override
    {
        return m_nid + ":" + Util::to_str<uint32_t>(m_sid) + ":" + Util::to_str<uint32_t>(m_cid) +
               Util::to_str<uint32_t>(m_lp) + ":" + Util::to_str<uint32_t>(m_rp) + ":" + m_ori;
    }

private:
    std::string m_nid;
    std::string m_full_name;
    uint32_t m_sid;
    uint32_t m_cid;
    std::string m_contig_name;
    uint32_t m_lp;
    uint32_t m_rp;
    char m_ori;
    std::string m_seq;
    bool m_is_tag;
    int m_fold;
};

class Node
{
public:
    Node() {}
    Node(int i, uint32_t alpha, uint32_t beta, uint32_t gamma, uint32_t delta) : uid{++uid_counter}, intervals{i}, init_alpha{alpha}, init_beta{beta}, init_gamma{gamma}, init_delta{delta}
    {
        oris[i] = '+';
        name = Util::to_str<int>(uid);
        auto key1 = get_key(init_alpha, init_beta);
        auto key2 = get_key(init_gamma, init_delta);
        common_numerals.insert(key1);
        common_numerals.insert(key2);
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
        sorted = n.sorted;
        oris = n.oris;
    }
    std::string GetFullName()
    {
        std::string fullName = "[";
        for (auto p : common_numerals)
        {
            fullName.append("(");
            fullName.append(Util::to_str<uint32_t>(p.first));
            fullName.append(",");
            fullName.append(Util::to_str<uint32_t>(p.second));
            fullName.append(")");
            fullName.append(",");
        }
        if (fullName.find_last_of(',') == fullName.size() - 1)
            fullName.erase(fullName.find_last_of(','));
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
        return "Node(a=" + Util::to_str(init_alpha) + ",b=" + Util::to_str(init_beta) +
               ",g=" + Util::to_str(init_gamma) + ",d=" + Util::to_str(init_delta) + ", name=" +
               GetName() + ",ivl=" + Util::convert_to_csv(intervals) + ",oris=" + Util::convert_to_csv(oris) + ")";
    }

    void add(int i, uint32_t alpha, uint32_t beta, uint32_t gamma, uint32_t delta,
             std::map<std::pair<uint32_t, uint32_t>, Node *> &lookup)
    {
        intervals.push_back(i);

        auto key_a_b = get_key(alpha, beta);
        auto key_g_d = get_key(gamma, delta);

        std::unordered_set<std::pair<uint32_t, uint32_t>> intersection{};

        set_intersection(intersection, key_a_b, key_g_d);

        common_numerals = intersection;

        auto key_a_b_init = get_key(init_alpha, init_beta);
        if (common_numerals.find(key_a_b_init) == common_numerals.end())
            lookup[key_a_b_init] = nullptr;

        auto key_g_d_init = get_key(init_gamma, init_delta);
        if (common_numerals.find(key_g_d_init) == common_numerals.end())
            lookup[key_g_d_init] = nullptr;
        
        if ((alpha == init_alpha && beta == init_beta) || (alpha == init_gamma && beta == init_delta) ||
            (gamma == init_alpha && delta == init_beta) || (gamma == init_gamma && delta == init_delta))
            oris[i] = '+';
        else if ((alpha == init_beta && beta == init_alpha) || (alpha == init_delta && beta == init_gamma) ||
                 (gamma == init_beta && delta == init_alpha) || (gamma == init_delta && delta == init_gamma))
            oris[i] = '-';
        else
            throw "No relation between interval numerals and node numerals.";
    }

private:
    void set_intersection(std::unordered_set<std::pair<uint32_t, uint32_t>>& intersection, 
        std::pair<uint32_t, uint32_t>& key_a_b, std::pair<uint32_t, uint32_t>& key_g_d)
    {
        for (const std::pair<uint32_t, uint32_t>& p: common_numerals) {
            if (p.first == key_a_b.first && p.second == key_a_b.second)
                intersection.insert(key_a_b);
            if (p.first == key_g_d.first && p.second == key_g_d.second)
                intersection.insert(key_g_d);
        }   
    }
    friend class CoPanGraph;
    std::string name;
    std::vector<int> intervals;
    std::unordered_map<int, char> oris;
    uint32_t init_alpha, init_beta, init_gamma, init_delta;
    int uid;
    std::unordered_set<std::pair<uint32_t, uint32_t>> common_numerals;
    bool sorted = false;
    static int uid_counter;
};

int Node::uid_counter = 0;
typedef std::map<std::pair<uint32_t, uint32_t>, Node *> NodeMap;

class CoPanGraph
{
private:
    std::vector<SequenceInterval> intervals;
    ContigContainerPtr contig_container;
    IdMapPtr cid_to_sid;
    const unsigned num_samples;
    std::vector<Node*> nodes{};
    std::vector<Node*> ivl_to_node{};
    NodeMap node_lookup{};

public:
    CoPanGraph(std::vector<SequenceInterval> &seq_intervals, ContigContainerPtr contig_container_ptr, IdMapPtr contig_to_sample_map, uint32_t num_of_samples) : intervals{seq_intervals}, contig_container{contig_container_ptr}, cid_to_sid{contig_to_sample_map}, num_samples{num_of_samples}
    {
    }
    uint32_t get_start_of(int i)
    {
        return intervals[i].start;
    }
    uint32_t get_end_of(int i)
    {
        return intervals[i].end;
    }
    char get_orientation_of(int i)
    {
        return ivl_to_node[i]->oris[i];
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
                ivl_to_node.push_back(node);
            }
            else
            {
                node = new Node(i, ivl.alpha, ivl.beta, ivl.gamma, ivl.delta);
                ivl_to_node.push_back(node);
                nodes.push_back(node);

                node_lookup[key_a_b] = node;
                node_lookup[key_g_d] = node;
            }
        }

        logger.Debug("Node_lookup size: " + Util::to_str<unsigned long>(node_lookup.size()));
        logger.Debug("ivl_to_node size: " + Util::to_str<unsigned long>(ivl_to_node.size()));
    }

    void output_graph()
    {
        auto rg_pref = config.GetValue<std::string>("graph_pref");
        auto output_dir = config.GetValue<std::string>("output_dir");
        auto file_name = output_dir + "/" + rg_pref;
        auto gfa_file = file_name+ config.GetValue<std::string>("gfa_file_ext");
        auto fasta_file = file_name + config.GetValue<std::string>("fasta_file_ext");
        auto node_table_file = file_name + config.GetValue<std::string>("node_color_file_ext");
        auto edge_table_file = file_name + config.GetValue<std::string>("edge_color_file_ext");

        logger.Info("writing to files"); 
        
        write_node_table(node_table_file);
        write_edge_table(edge_table_file);
        write_gfa(gfa_file);
        write_fasta(fasta_file);
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
                auto cid = get_cid_of(i);
                
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

        for (const pElement &e : seg_type_list)
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

    std::shared_ptr<Element> build_elem(std::string &nid, int cid, uint32_t start, uint32_t end, const std::string &full_name,
                                        char ori, int sid, const std::string &format)
    {
        auto &extended_contig = contig_container->at(cid);
        auto start_of_contig = extended_contig.GetBwdExt();
        auto end_of_contig = extended_contig.GetBwdExt() + extended_contig.GetUnExtLen() + 1;
        bool is_tag = extended_contig.GetIsTag();
        std::string seq = extended_contig.Seq.substr(start, end-start);
        std::string s_start, s_end;
       
        if (end < start_of_contig)
        {
            s_start = "0-" + Util::to_str<uint32_t>(abs(start - start_of_contig));
            s_end = "0-" + Util::to_str<uint32_t>(abs(end - start_of_contig));
            is_tag = true;
        }
        else if (start < start_of_contig && start_of_contig <= end && end < end_of_contig)
        {
            s_start = "0-" + Util::to_str<uint32_t>(abs(start - start_of_contig));
            s_end = Util::to_str<uint32_t>(end - start_of_contig);
            is_tag = true;
        }
        else if (start_of_contig <= start && start < end_of_contig &&
                 start_of_contig <= end && end < end_of_contig)
        {
            s_start = Util::to_str<uint32_t>(start - start_of_contig);
            s_end = Util::to_str<uint32_t>(end - start_of_contig);
        }
        else if (start_of_contig <= start && start < end_of_contig && end >= end_of_contig)
        {
            s_start = Util::to_str<uint32_t>(start - start_of_contig);
            s_end = Util::to_str<uint32_t>(extended_contig.GetUnExtLen()) + "+" +
                    Util::to_str<uint32_t>(end - end_of_contig);
            is_tag = true;
        }
        else if (start >= end_of_contig)
        {
            s_start = Util::to_str<uint32_t>(extended_contig.GetUnExtLen()) + "+" +
                      Util::to_str<uint32_t>(start - end_of_contig);
            s_end = Util::to_str<uint32_t>(extended_contig.GetUnExtLen()) + "+" +
                    Util::to_str<uint32_t>(end - end_of_contig);
            is_tag = true;
        }
        else if (start < start_of_contig && end >= end_of_contig)
        {
            s_start = "0-" + Util::to_str<uint32_t>(start - start_of_contig);
            s_end = Util::to_str<uint32_t>(extended_contig.GetUnExtLen()) +
                    "+" + Util::to_str<uint32_t>(end - end_of_contig);
            is_tag = true;
        }
        else
            throw("Edge case escaped.");

        Element *e;

        if (format == "fasta")
        {
            auto cns = CopanNodeSeq(nid, full_name, sid, cid, extended_contig.Name, start, end, ori, is_tag, seq);
            return std::make_shared<CopanNodeSeq>(cns);
        }
        else if (format == "gfa")
        {
            Rest rest;
            rest.contig_name = extended_contig.Name;
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