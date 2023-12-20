#ifndef GLUEPOINTS_H
#define GLUEPOINTS_H

#include "alignment.h"
#include "async_logger.h"
#include <vector>
#include "aux_functions.h"
#include "ext_contig.h"
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <memory>
#include <ostream>
#include "set.h"
#include "global_copan.h"

enum PointType
{
    START_RIGHT = 0,
    START_LEFT = 1,
    END_RIGHT = 2,
    END_LEFT = 3
};

struct Point
{
    Point(unsigned int q_cid, unsigned int s_cid, unsigned int q_pos, unsigned int s_pos, const PyAlignment &alignment, PointType pt) : q_cid(q_cid), s_cid(s_cid), q_pos(q_pos), s_pos(s_pos), alignment(alignment), pointType{pt} {}
    unsigned int q_cid, s_cid, q_pos, s_pos;
    const PyAlignment &alignment;
    PointType pointType;
    std::string ToString() const
    {
        return std::to_string(q_cid) + "," + std::to_string(s_cid) + "," + std::to_string(q_pos) + "," + std::to_string(s_pos) +
            std::to_string(pointType);
    }
};

struct SequenceInterval
{
    SequenceInterval(uint64_t a, uint64_t b, uint64_t g, uint64_t d, uint32_t c, uint64_t s, uint32_t e) : alpha(a), beta(b), gamma(g), delta(d), cid(c), start(s), end(e) {}
    uint64_t alpha, beta, gamma, delta;
    uint32_t cid, start, end;

    void print()
    {
        std::cout << "a: " << alpha << " b: " << beta << " g: " << gamma << " d: " << delta << " start: " << start << " end: " << end << std::endl;
    }

    std::string ToString()
    {
        return "a:" + std::to_string(alpha) + " b:" + std::to_string(beta) + " g:" + std::to_string(gamma) +
               " d:" + std::to_string(delta);
    }
};

template <class T>
struct Cluster
{
    Cluster(T const &d, bool cluster_id) : data(d), cluster_id(cluster_id){};
    bool cluster_id;
    T data;
    std::string ToString()
    {
        return cluster_id ? "true" : "false";
    }
};
struct Breakpoint
{
    Breakpoint(uint32_t cid, uint32_t pos, PointType type, bool q_derived) : cid(cid), pos(pos), pointType(type), q_derived(q_derived) {}
    uint32_t cid, pos;
    bool q_derived;
    PointType pointType;
    void print()
    {
        std::cout << "cid: " << cid << " pos: " << pos << " type: " << pointType << " q_derived: " << q_derived << std::endl;
    }
    std::string ToString()
    {
        return "cid: " + std::to_string(cid) + " pos: " + std::to_string(pos) + " type: " + std::to_string(pointType) + " q_derived: " + std::to_string(q_derived);
    }
};
struct Gluepoint
{
    uint64_t sr, sl, er, el;
    uint32_t cid, pos;

    Gluepoint() {}
    Gluepoint(uint64_t init_gpid, uint32_t _cid, uint32_t _pos) : cid(_cid), pos(_pos), sr(0), sl(0), er(0), el(0)
    {
        switch (init_gpid % 4)
        {
        case PointType::START_RIGHT:
            sr = init_gpid;
            break;
        case PointType::START_LEFT:
            sl = init_gpid;
            break;
        case PointType::END_RIGHT:
            er = init_gpid;
            break;
        case PointType::END_LEFT:
            el = init_gpid;
            break;
        }
    }

    void print() const
    {
        std::cout << sr << ' ' << sl << ' ' << er << ' ' << el << ' ' << pos << ' ' << cid << std::endl;
    }
    std::string ToString()
    {
        return std::to_string(cid) + ": " + std::to_string(sr) + "," + std::to_string(sl) + "," + std::to_string(er) + "," + std::to_string(el);
    }
};

uint32_t abs(uint32_t a, uint32_t b)
{
    return (a > b) * (a - b) + (a <= b) * (b - a);
}
struct ClusterGenerator
{

    ClusterGenerator(std::vector<Cluster<Point const *>> const &c) : clusters(c)
    {
        current = clusters.begin();
    }
    std::vector<Cluster<Point const *>> clusters;
    std::vector<Cluster<Point const *>>::iterator current;

    void get_cluster(std::vector<Cluster<Point const *>>::iterator &left, std::vector<Cluster<Point const *>>::iterator &right)
    {
        left = current;
        right = left;
        while (right < clusters.end() && left->cluster_id == right->cluster_id)
        {
            right++;
        }
        current = right;
    }
};

template <typename Iterator, typename Lambda>
unsigned int median(Iterator const &begin, Iterator const &end, Lambda &&extract_val)
{
    assert((end - begin) != 0);
    size_t n = end - begin;
    if (n % 2 == 1)
    { // odd
        return extract_val(begin + n / 2);
    }
    else
    {
        unsigned int a = extract_val(begin + n / 2);
        unsigned int b = extract_val(begin + n / 2 - 1);
        return (a + b) / 2;
    }
};

typedef std::vector<Gluepoint> GluepointList;
typedef std::unordered_map<unsigned int, GluepointList> GluepointMap;
typedef std::shared_ptr<GluepointMap> GluepointMapPtr;

auto same_cluster_q = [](Point const & a, Point const & b, uint32_t max_separation) -> bool {
    return ((a.q_pos > b.q_pos) && (a.q_pos - b.q_pos <= max_separation)) || ((a.q_pos <= b.q_pos) && (b.q_pos - a.q_pos <= max_separation));
};

auto same_cluster_s = [](Point const & a, Point const & b, uint32_t max_separation) -> bool {
    return (a.s_cid == b.s_cid) && (((a.s_pos > b.s_pos) && (a.s_pos - b.s_pos <= max_separation)) || ((a.s_pos <= b.s_pos) && (b.s_pos - a.s_pos <= max_separation)));
};

unsigned int abs_diff(unsigned int a, unsigned int b)
{
    return (a > b) * (a - b) + (a <= b) * (b - a);
}

std::shared_ptr<std::unordered_map<unsigned int, std::vector<Point>>> construct_endpoints(std::vector<PyAlignment> &alignments)
{
    std::unordered_map<unsigned int, std::vector<Point>> endpoints;
    endpoints.reserve(alignments.size());
    for (const auto &a : alignments)
    {

        unsigned int q_cid = a.q_cid;
        unsigned int s_cid = a.s_cid;
        unsigned int q_begin = a.q_begin;
        unsigned int s_begin = a.s_begin;
        unsigned int q_end = a.q_end;
        unsigned int s_end = a.s_end;
        
        // endpoints[s_cid].emplace_back(s_cid, q_cid, s_begin, q_begin, a, not is_end);
        // endpoints[s_cid].emplace_back(s_cid, q_cid, s_end, q_end, a, is_end);

        if (s_begin < s_end)
        {
            endpoints[q_cid].emplace_back(q_cid, s_cid, q_begin, s_begin, a, PointType::START_RIGHT);
            endpoints[q_cid].emplace_back(q_cid, s_cid, q_end, s_end, a, PointType::END_LEFT);
        }
        else
        { // then we have reverse complement
            endpoints[q_cid].emplace_back(q_cid, s_cid, q_begin, s_begin, a, PointType::END_RIGHT);
            endpoints[q_cid].emplace_back(q_cid, s_cid, q_end, s_end, a, PointType::START_LEFT);
        }
        
    }
    return std::make_shared<std::unordered_map<unsigned int, std::vector<Point>>>(endpoints);
}

template <typename Lambda>
void cluster_points(std::vector<Cluster<Point const *>> &clusters, bool &cluster_id, unsigned int max_separation,
                    Lambda &&same_cluster, bool complete_link)
{
    // loop through sorted endpoints and cluster
    size_t left_bound = 0;
    for (size_t i = 1; i < clusters.size(); i++) {
        // if distance less than max_separation, point belongs to same cluster, otherwise assign it to new cluster
        if (clusters[i].data->pointType == clusters[left_bound].data->pointType && same_cluster(*(clusters[i].data), *(clusters[left_bound].data), max_separation) ) {
            clusters[i].cluster_id = cluster_id;
            if (not complete_link) {
                left_bound = i;
            }
        }
        else {
            cluster_id = not cluster_id; // flip to denote a new boundary
            clusters[i].cluster_id = cluster_id;
            left_bound = i;
        }
    }
    cluster_id = not cluster_id;
}

std::vector<Cluster<Point const *>> cluster_endpoints(std::unordered_map<unsigned int, std::vector<Point>> &endpoints,
                                                      unsigned int max_separation = 75, bool complete_link = true)
{
    std::vector<Cluster<Point const *>> all_clusters;
    // get total number of endpoints for reserve operation
    size_t num_endpoints = 0;
    for (auto const &[k, v] : endpoints)
    {
        num_endpoints += v.size();
    }
    all_clusters.reserve(num_endpoints);
    // iterate through endpoints and cluster
    bool cluster_bound = false;
    for (auto const &[cid, eps] : endpoints)
    {
        // Set up the cluster vector. Each Cluster.data points to a different Point in eps - the endpoints for contig cid
        std::vector<Cluster<Point const *>> clusters(eps.size(), Cluster<Point const *>(eps.data(), cluster_bound));
        for (size_t i = 0; i < clusters.size(); i++)
        {
            clusters[i].data += i;
        }
        std::sort(clusters.begin(), clusters.end(), [](Cluster<Point const *> &a, Cluster<Point const *> &b)
                  {
                if (a.data->pointType != b.data->pointType) {
                    return a.data->pointType < b.data->pointType; 
                }
                else {
                    return a.data->q_pos < b.data->q_pos; // then by q_pos
                } });
        cluster_points(clusters, cluster_bound, max_separation, same_cluster_q, complete_link);
        all_clusters.insert(all_clusters.end(), clusters.begin(), clusters.end());
    }
    return all_clusters;
}

void initialize(
    ContigContainer &contigs,
    std::unordered_map<unsigned int, std::vector<SetNode<Breakpoint> *>> &lookup)
{

    // add start and end breakpoints to gluepoints
    bool is_end = true;
    lookup.reserve(contigs.size());

    for (const auto &c : contigs)
    {
        unsigned int cid = c.GetId();
        auto start = new SetNode<Breakpoint>(Breakpoint(cid, 0, PointType::START_RIGHT, true));
        auto end = new SetNode<Breakpoint>(Breakpoint(cid, c.GetLen(), PointType::END_LEFT, true));

        auto &bp_on_contig = lookup[cid];
        bp_on_contig.emplace_back(start);
        bp_on_contig.emplace_back(end);
    }
}

std::vector<Point const *> find_overlapping_alignments(Breakpoint const &bp, std::vector<Point> const &endpoints,
                                                       unsigned int max_separation)
{
    std::vector<Point const *> overlapping_alignments;
    for (size_t i = 0; i < endpoints.size(); i += 2) {
        
        if (
            (abs(endpoints[i + 1].q_pos, bp.pos) > max_separation) &&
            (abs(bp.pos, endpoints[i].q_pos) > max_separation) &&
            (((endpoints[i].q_pos < bp.pos - 1) && (bp.pos - 1 < endpoints[i+1].q_pos)) ||
            ((endpoints[i].q_pos < bp.pos + 1) && (bp.pos + 1 < endpoints[i+1].q_pos)))
        ) {
            overlapping_alignments.emplace_back(&endpoints[i]);
        }
    }
    return overlapping_alignments;
}

unsigned int get_split_position(uint32_t q_pos, PyAlignment const &a)
{
    unsigned int q_begin = a.q_begin;
    unsigned int s_begin = a.s_begin;
    unsigned int q_end = a.q_end;
    unsigned int s_end = a.s_end;

     int complement_handler = (s_begin <= s_end) ? 1 : -1;
    if (complement_handler == 1) {
        if (q_pos <= q_begin)
            return s_begin;
        if (q_pos >= q_end)
            return s_end;
    }
    else {
        if (q_pos <= q_begin)
            return s_end;
        if (q_pos >= q_end)
            return s_begin;
    }
    uint32_t projected_pos;
    if (true /*q_info.size == 0*/) {
        double stretch_ratio = abs(s_end, s_begin) / static_cast<double> (q_end - q_begin);
        projected_pos = s_begin + complement_handler * stretch_ratio * (q_pos - q_begin);
        if (complement_handler == 1) {
            return std::max(s_begin, std::min(projected_pos, s_end));
        }
        else {
            return std::max(s_end, std::min(projected_pos, s_begin));
        }
    }
    else {
        //auto ptr = upper_bound(
        //    q_match_pos, q_match_pos + q_info.size,
        //    q_pos, [](uint32_t const & value, uint32_t const & element) {return value <= element;}
        //);
        //size_t i  = ptr - q_match_pos;
        //assert (q_match_pos < ptr && ptr < q_match_pos + q_info.size);
        //double stretch_ratio = static_cast<double>(s_match_pos[i] - s_match_pos[i-1]) / (q_match_pos[i] - q_match_pos[i-1]);
        //uint32_t projected_pos = s_match_pos[i-1] + (q_pos - q_match_pos[i-1]) * stretch_ratio;
        //return max(s_match_pos[i-1], min(projected_pos, s_match_pos[i]));
   }
}

void merge_gp(std::vector<SetNode<Breakpoint> *> &gluepoints,
              std::unordered_map<unsigned int, std::vector<SetNode<Breakpoint> *>> &lookup, unsigned int max_separation, PointType pointType)
{
    // get the first point with type is_end

    for (auto bp : gluepoints) {
        // only process the breakpoints of they match is_end
        if (bp->data.pointType != pointType) {
            continue;
        }
        auto & bp_on_contig = lookup[bp->data.cid];
        auto p = upper_bound(
            bp_on_contig.begin(), bp_on_contig.end(),
            bp->data, [](Breakpoint const & value, SetNode<Breakpoint> * const & element) {
                if (value.pointType != element->data.pointType) {
                    return value.pointType < element->data.pointType;
                }
                else {
                    return value.pos  <= element->data.pos;
                }
            }
        );

        if (p > bp_on_contig.begin() && abs(bp->data.pos, (*(p-1))->data.pos) <= max_separation && bp->data.pointType == (*(p-1))->data.pointType) {
            unionSetMerge(bp, *(p-1));
            //bp->data.q_derived = bp->data.q_derived || (*(p-1))->data.q_derived;
            //(*(p-1))->data.q_derived = bp->data.q_derived;
        }
        if (p < bp_on_contig.end() && abs((*p)->data.pos, bp->data.pos) <= max_separation && bp->data.pointType == (*p)->data.pointType) {
            unionSetMerge(bp, *p);
            //bp->data.q_derived = bp->data.q_derived || (*p)->data.q_derived;
            //(*p)->data.q_derived = bp->data.q_derived;
        }

        bp_on_contig.insert(p, bp);
    }

    auto union_ptr = gluepoints.begin();
    for (; union_ptr != gluepoints.end(); union_ptr++) {
        if ((*union_ptr)->data.pointType == pointType) {
            break;
        }
    }

    // return early if theres no points of type
    if  (union_ptr == gluepoints.end()) {
        return;
    }

    // otherwise merge the points with the first identified point
    for (auto p = gluepoints.begin(); p < gluepoints.end(); p++) {
        // only merge if the end types are the same, and both match is_end
        if ((*p)->data.pointType == pointType) {
            unionSetMerge(*union_ptr, *p);
        }
    };
}
void replace_with_root_pointer(std::unordered_multimap<SetNode<Breakpoint> *, SetNode<Breakpoint> *> &equivalent_points)
{
     // get the unique set of keys
    std::unordered_set<SetNode<Breakpoint> *> keys;
    for (auto const & [k, v] : equivalent_points) {
        keys.insert(k);
    }

    std::unordered_multimap<SetNode<Breakpoint> *, SetNode<Breakpoint> *> testing;

    // for each key:
    //  - get the set of k,v pairs associated with it
    //  - remove the k,v pairs from equiv points
    //  - replace each k,v pair with the set of k, v
    for (auto k : keys) {

        auto range = equivalent_points.equal_range(k);

        // get k, v pairs
        std::vector<std::pair<SetNode<Breakpoint> *, SetNode<Breakpoint> *>> kv_pairs;
        for (auto it = range.first; it != range.second; it++) {
            kv_pairs.emplace_back(*it);
        }

        // remove all kv pairs
        equivalent_points.erase(range.first, range.second);

        // replace with root set
        for (auto [k, v] : kv_pairs) {
            auto root_k = findSetMerge(k);
            auto root_v = findSetMerge(v);
            bool exists = false;
            auto range = equivalent_points.equal_range(root_k);
            for (auto it = range.first; it != range.second; it++) {
                if (it->second == root_v) {
                    exists = true;
                }
            }
            if (!exists) {
                equivalent_points.insert({root_k, root_v});
            }
        }
    }
}

void update_equivalent_points(SetNode<Breakpoint> *inducing_point,
                              PointType equivalent_point_type, std::vector<SetNode<Breakpoint> *> const &gluepoint,
                              std::unordered_set<SetNode<Breakpoint> *> &equivalent_points)
{

    // Scan breakpoints to find  one example of point type if available. Choice doesn't matter, as they'll
    // all of the same type will be in the same set
    SetNode<Breakpoint> * equivalent_point = nullptr;
    for (auto bp : gluepoint) {
        if (bp->data.pointType == equivalent_point_type) {
            equivalent_point = bp; break;
        }
    }
    unionSetEqvPt(inducing_point, equivalent_point);
    equivalent_points.insert(inducing_point);
    equivalent_points.insert(equivalent_point);



    // check if alrealy exists, and if not, add. This is purely to keep the size down.
    //equivalent_points.insert({inducing_point, equivalent_point});
    //equivalent_points.insert({equivalent_point, inducing_point});
    //bool exists = false;
    //auto range = equivalent_points.equal_range(inducing_point);
    //for (auto it = range.first; it != range.second; it++) {
    //    if (it->second == equivalent_point) {
    //        exists = true;
    //    }
    //}
    //if (!exists) {
    //    equivalent_points.insert({inducing_point, equivalent_point});
    //}

    //exists = false;
    //range = equivalent_points.equal_range(equivalent_point);
    //for (auto it = range.first; it != range.second; it++) {
    //    if (it->second == inducing_point) {
    //        exists = true;
    //    }
    //}
    //if (!exists) {
    //    equivalent_points.insert({equivalent_point, inducing_point});
    //}
}

void connect_roots_through_equivalent_points(std::unordered_set<SetNode<Breakpoint> *> &equivalent_points)
{
    // std::vector<SetNode<Breakpoint> *> eq_nodes_vector;
    // eq_nodes_vector.reserve(equivalent_points.size());
    
    // for (auto& ep: equivalent_points)
    //     eq_nodes_vector.push_back(ep);
    
    // std::sort(eq_nodes_vector.begin(), eq_nodes_vector.end(), [] (SetNode<Breakpoint>* lhs, SetNode<Breakpoint>*rhs)
    // {
    //     return lhs->data.cid < rhs->data.cid || lhs->data.cid == rhs->data.cid && lhs->data.pointType < rhs->data.pointType ||
    //         lhs->data.cid == rhs->data.cid && lhs->data.pointType == rhs->data.pointType && lhs->data.pos < rhs->data.pos ||
    //         lhs->data.cid == rhs->data.cid && lhs->data.pointType == rhs->data.pointType && lhs->data.pos == rhs->data.pos &&
    //             lhs->data.q_derived < rhs->data.q_derived;
    // }
    // );

    for (auto ptr : equivalent_points) 
    {  
        auto root_m = findSetMerge(ptr);
        auto root_m_ep = findSetMerge(findSetEqvPt(ptr));
        if (root_m < root_m_ep)
        {
            root_m_ep->parent_merge = root_m;
        }
        else {
            root_m->parent_merge = root_m_ep;
        }       
        
    }
}

void process_overlap(
    SetNode<Breakpoint> *const inducing_point,
    Point const *const a,
    std::unordered_set<SetNode<Breakpoint> *> &equivalent_points,
    std::unordered_map<uint32_t, std::vector<SetNode<Breakpoint> *>> &lookup,
    uint32_t max_separation)
{
    
    // determine alignment complement, position of inducing point in s, and s_cid
    auto is_fc = [](const PyAlignment &alignment)
    {
        return alignment.s_begin <= alignment.s_end;
    };
    bool is_fwd = is_fc(a->alignment);
    uint32_t s_pos = get_split_position(inducing_point->data.pos, a->alignment);
    uint32_t s_cid = static_cast<uint32_t>(a->alignment.s_cid);
    
    // compute merge and equivalence point under split logic
    if (is_fwd) {
        auto qsr = new SetNode<Breakpoint>(Breakpoint(inducing_point->data.cid, inducing_point->data.pos, PointType::START_RIGHT, true));
        auto qel = new SetNode<Breakpoint>(Breakpoint(inducing_point->data.cid, inducing_point->data.pos, PointType::END_LEFT, true));
        auto ssr = new SetNode<Breakpoint>(Breakpoint(s_cid, s_pos, PointType::START_RIGHT, false));
        auto sel = new SetNode<Breakpoint>(Breakpoint(s_cid, s_pos, PointType::END_LEFT, false));
        std::vector<SetNode<Breakpoint> *> split_breakpoints = {qsr, qel, ssr, sel};
        switch (inducing_point->data.pointType) {
            case PointType::START_RIGHT:
            case PointType::END_LEFT:
                split_breakpoints.push_back(inducing_point);
                break;
            case PointType::START_LEFT:
                update_equivalent_points(inducing_point, PointType::END_LEFT, split_breakpoints,equivalent_points);
                break;
            case PointType::END_RIGHT:
                update_equivalent_points(inducing_point, PointType::START_RIGHT, split_breakpoints,equivalent_points);
                break;
        }
        merge_gp(split_breakpoints, lookup, max_separation, PointType::START_RIGHT);
        merge_gp(split_breakpoints, lookup, max_separation, PointType::END_LEFT);

    } else {
        auto qsl = new SetNode<Breakpoint>(Breakpoint(inducing_point->data.cid, inducing_point->data.pos, PointType::START_LEFT, true));
        auto qer = new SetNode<Breakpoint>(Breakpoint(inducing_point->data.cid, inducing_point->data.pos, PointType::END_RIGHT, true));
        auto ssr = new SetNode<Breakpoint>(Breakpoint(s_cid, s_pos, PointType::START_RIGHT, false));
        auto sel = new SetNode<Breakpoint>(Breakpoint(s_cid, s_pos, PointType::END_LEFT, false));
        std::vector<SetNode<Breakpoint> *> split_breakpoints = {qsl, qer, ssr, sel};
        switch (inducing_point->data.pointType) {
            case PointType::START_RIGHT:
                update_equivalent_points(inducing_point, PointType::END_RIGHT, split_breakpoints,equivalent_points);
                update_equivalent_points(inducing_point, PointType::END_LEFT, split_breakpoints,equivalent_points);
                break;
            case PointType::END_LEFT:
                update_equivalent_points(inducing_point, PointType::START_LEFT, split_breakpoints,equivalent_points);
                update_equivalent_points(inducing_point, PointType::START_RIGHT, split_breakpoints,equivalent_points);
                break;
            case PointType::START_LEFT:
                split_breakpoints.emplace_back(inducing_point);
                update_equivalent_points(inducing_point, PointType::START_RIGHT, split_breakpoints,equivalent_points);
                break;
            case PointType::END_RIGHT:
                split_breakpoints.emplace_back(inducing_point);
                update_equivalent_points(inducing_point, PointType::END_LEFT, split_breakpoints,equivalent_points);
                break;
        }
        update_equivalent_points(qsl, PointType::START_RIGHT, split_breakpoints, equivalent_points);
        update_equivalent_points(qer, PointType::END_LEFT, split_breakpoints, equivalent_points);
        merge_gp(split_breakpoints, lookup, max_separation, PointType::START_RIGHT);
        merge_gp(split_breakpoints, lookup, max_separation, PointType::START_LEFT);
        merge_gp(split_breakpoints, lookup, max_separation, PointType::END_RIGHT);
        merge_gp(split_breakpoints, lookup, max_separation, PointType::END_LEFT);        
    }

    
}

void group_breakpoints_into_gluepoints(std::unordered_map<uint32_t, std::vector<Point>> &endpoints,
                                       std::unordered_set<SetNode<Breakpoint> *> &equivalent_points,
                                       std::vector<Cluster<Point const *>> const &clusters, std::vector<PyAlignment> &alignments,
                                       std::unordered_map<uint32_t, std::vector<SetNode<Breakpoint> *>> &lookup, int32_t const max_separation, bool complete_link)
{
    // loop through each cluster construct gluepoints
    ClusterGenerator q_gen(clusters);
    std::vector<Cluster<Point const *>>::iterator q_cluster_l, q_cluster_r;
    q_gen.get_cluster(q_cluster_l, q_cluster_r);
    
    while (q_cluster_r - q_cluster_l > 0) {

        // construct inducing point on q
        std::vector<SetNode<Breakpoint> *> gluepoint;
        uint32_t median_point = median(q_cluster_l, q_cluster_r, [](std::vector<Cluster<Point const *>>::const_iterator const & a) -> uint32_t {return a->data->q_pos;});
        uint32_t q_cid = q_cluster_l->data->q_cid;   
        
        auto inducing_point = new SetNode<Breakpoint>(Breakpoint(q_cid, median_point, q_cluster_l->data->pointType, true));
        gluepoint.emplace_back(inducing_point);

        // construct projected points on s
        bool cluster_bound = false;
        std::vector<Cluster<Point const *>> s_clusters;
        s_clusters.reserve((q_cluster_r - q_cluster_l) + 128);
        for (auto p = q_cluster_l; p < q_cluster_r; p++) {
            s_clusters.emplace_back(p->data, cluster_bound);
        }
        std::sort(
            s_clusters.begin(), s_clusters.end(), [](Cluster<Point const *> &a, Cluster<Point const *> &b) {
                if (a.data->s_cid != b.data->s_cid) {
                    return a.data->s_cid < b.data->s_cid;
                }
                return a.data->s_pos < b.data->s_pos;
            }
        );

        cluster_points(s_clusters, cluster_bound, max_separation, same_cluster_s, complete_link);
        ClusterGenerator s_gen(s_clusters);
        std::vector<Cluster<Point const *>>::iterator s_cluster_l, s_cluster_r;
        s_gen.get_cluster(s_cluster_l, s_cluster_r);
        
        while ((s_cluster_r - s_cluster_l) > 0) {
            median_point = median(s_cluster_l, s_cluster_r, [](std::vector<Cluster<Point const *>>::const_iterator const & a) -> uint32_t {return a->data->s_pos;});
            uint32_t s_cid = s_cluster_l->data->s_cid;

            // if the alignment was the reverse complement of q, the point type must be inverted when projected to s
            PointType type;
            if (q_cluster_l->data->pointType == PointType::START_LEFT) {
                type = PointType::START_RIGHT;
            }
            else if (q_cluster_l->data->pointType == PointType::END_RIGHT) {
                type = PointType::END_LEFT;
            }
            else {
                type = q_cluster_l->data->pointType;
            }
            auto s_bp = new SetNode<Breakpoint>(Breakpoint(s_cid, median_point, type, false));
            gluepoint.emplace_back(s_bp);
            s_gen.get_cluster(s_cluster_l, s_cluster_r);
        }
        // merge points, and set equivalence
        switch(inducing_point->data.pointType) {
            case PointType::START_RIGHT:
                merge_gp(gluepoint, lookup, max_separation, PointType::START_RIGHT);
                break;
            case PointType::START_LEFT:
                merge_gp(gluepoint, lookup, max_separation, PointType::START_RIGHT);
                merge_gp(gluepoint, lookup, max_separation, PointType::START_LEFT);
                update_equivalent_points(inducing_point, PointType::START_RIGHT, gluepoint, equivalent_points);
                break;
            case PointType::END_RIGHT:
                merge_gp(gluepoint, lookup, max_separation, PointType::END_RIGHT);
                merge_gp(gluepoint, lookup, max_separation, PointType::END_LEFT);
                update_equivalent_points(inducing_point, PointType::END_LEFT, gluepoint, equivalent_points);
                break;
            case PointType::END_LEFT:
                merge_gp(gluepoint, lookup, max_separation, PointType::END_LEFT);
                break;
        }
        
        // split alignments (solve for in-transitivity)
        auto overlapping_alignments = find_overlapping_alignments(inducing_point->data, endpoints[inducing_point->data.cid], max_separation);
        for (auto a : overlapping_alignments) {
            process_overlap(inducing_point, a, equivalent_points, lookup, max_separation);
        }

        q_gen.get_cluster(q_cluster_l, q_cluster_r);
    }
}

void add_consensus_gluepoint(std::vector<SetNode<Breakpoint> *>::const_iterator left_arr,
                             std::vector<SetNode<Breakpoint> *>::const_iterator right_arr, std::vector<Gluepoint> &consensus_gps,
                             unsigned int max_separation, uint64_t gp_id)
{
    size_t cluster_sz = (*(right_arr - 1))->data.pos - (*left_arr)->data.pos;
    unsigned int cid = (*left_arr)->data.cid;
    
    if (cluster_sz > max_separation)
    {
        consensus_gps.emplace_back(gp_id, cid, (*left_arr)->data.pos);
        unsigned int repeats = std::floor(cluster_sz / max_separation);
        unsigned int mode = cluster_sz / repeats;
        for (size_t i = 1; i < repeats; i++)
        {
            unsigned int pos = (*left_arr)->data.pos + i * mode;
            consensus_gps.emplace_back(gp_id, cid, pos);
        }
        consensus_gps.emplace_back(gp_id, cid, (*(right_arr - 1))->data.pos);
        
    }
    else
    {
        unsigned int consensus_pos = median(left_arr, right_arr, [](std::vector<SetNode<Breakpoint> *>::const_iterator const &a) -> unsigned int
                                            { return (*a)->data.pos; });
        consensus_gps.emplace_back(gp_id, cid, consensus_pos);
        
    }
}

void merge_s_only_clusters(std::unordered_map<uint32_t, std::vector<Point>> &endpoints,
                           std::unordered_map<uint32_t, std::vector<SetNode<Breakpoint> *>> &lookup,
                           std::unordered_set<SetNode<Breakpoint> *> &equivalent_points, uint32_t max_separation)
{
    // set q-derived variable across all points
    for (auto [cid, bps] : lookup) {
        auto left = bps.begin();
        auto right = bps.begin();
        while (right != bps.end()) {
            bool q_derived = (*right)->data.q_derived;
            while (right != bps.end() && findSetMerge(*left) == findSetMerge(*right)) {
                q_derived = q_derived || (*right)->data.q_derived;
                right++;
            }

            // if at least one element in the cluster-group is q-derived
            // set all of them to be q-derived
            if (!q_derived) {
                auto overlapping_alignments = find_overlapping_alignments((*left)->data, endpoints[cid], max_separation);
                for (auto& a : overlapping_alignments) {
                    process_overlap(*left, a, equivalent_points, lookup, max_separation);
                }
                
            }
            left = right;
        }
    }
}

std::vector<Gluepoint> build_consensus_gluepoints(
    std::unordered_map<uint32_t, std::vector<SetNode<Breakpoint> *>> &lookup,
    std::unordered_map<SetNode<Breakpoint> *, uint64_t> &set_to_id,
    std::unordered_map<uint64_t, SetNode<Breakpoint> *> &id_to_set,
    uint32_t max_separation)
{
    std::vector<Gluepoint> consensus_gps;
    uint64_t id_gen_sr = 4;
    uint64_t id_gen_sl = 5;
    uint64_t id_gen_er = 6;
    uint64_t id_gen_el = 7;    

    uint64_t *id_gen_ptr;
    for (auto & [cid, bp_on_contigs] : lookup) {
        auto left_arr = bp_on_contigs.begin(), right_arr = bp_on_contigs.begin() + 1;
        // group and process the breakpoints
        while (left_arr != bp_on_contigs.end()) {
            // find group bound by right_arr extension --- is this not interacting badly with complete linkage?

            while (
                right_arr != bp_on_contigs.end() &&
                (*right_arr)->data.pos - (*(right_arr - 1))->data.pos <= max_separation &&
                (*right_arr)->data.pointType == (*(right_arr - 1))->data.pointType
            ) {
                right_arr++;
            }

            // Get the id of the gluepoint. If it's a new point, then the dictionary lookup will be 0
            // in which case, we need to assign a new gp_id to the glue point.
            auto set_ptr = findSetMerge(*left_arr);
            uint64_t & gp_id = set_to_id[set_ptr];
            if (gp_id == 0) {
                switch ((*left_arr)->data.pointType) {
                    case PointType::START_RIGHT:
                        id_gen_ptr = &id_gen_sr; break;
                    case PointType::START_LEFT:
                        id_gen_ptr = &id_gen_sl; break;
                    case PointType::END_RIGHT:
                        id_gen_ptr = &id_gen_er; break;
                    case PointType::END_LEFT:
                        id_gen_ptr = &id_gen_el; break;
                }
                gp_id = *id_gen_ptr;
                *id_gen_ptr += 4;
            }
            id_to_set[gp_id] = set_ptr;
            add_consensus_gluepoint(left_arr, right_arr, consensus_gps, max_separation, gp_id);

            // continue iteration
            left_arr = right_arr;
            right_arr++;
        }
        // empty the vector
        bp_on_contigs.clear();
    }


    // construct the final gluepoints. Gluepoints of different type, but at distance < ms away from one another
    // will be grouped.
    auto min_breakpoint_comp = [&consensus_gps] (std::vector<Gluepoint>::iterator a, std::vector<Gluepoint>::iterator b) {
        if (a == consensus_gps.end() and b != consensus_gps.end()) {
            return false;
        }
        if (a != consensus_gps.end() and b == consensus_gps.end()) {
            return true;
        }
        if (a == consensus_gps.end() and b == consensus_gps.end()) {
            return false;
        }
        return a->pos < b->pos;
    };

    auto abs_diff = [] (uint32_t const & a, uint32_t const & b) {
        return (a>b)*(a - b) + (b>=a)*(b - a);
    };

    auto contig_start = consensus_gps.begin(), contig_end = consensus_gps.begin();
    while (contig_start != consensus_gps.end()) {

        // Move pointers to bound the start/end of the next contig
        while (contig_end != consensus_gps.end() && contig_start->cid == contig_end->cid) {
            contig_end++;
        }

        // Find the bounds of the types within the contig bounds, type order is sr, sl, er, el
        auto sr_end = contig_start, sl_end = contig_start, er_end = contig_start, el_end = contig_end;
        while (sr_end != contig_end && (sr_end->sr != 0)) sr_end++;
        sl_end = sr_end;
        while (sl_end != contig_end && (sl_end->sl != 0)) sl_end++;
        er_end = sl_end;
        while (er_end != contig_end && (er_end->er != 0)) er_end++;

        // begin assignment walk. Algorithm makes a pass through each point. If points of different type are
        // <= max_separation away the smallest point then we group data into the same breakpoint and move the
        // iterator forward
        auto sr_it = contig_start, sl_it = sr_end, er_it = sl_end, el_it = er_end;
        while (sr_it != consensus_gps.end() || sl_it != consensus_gps.end() || er_it != consensus_gps.end() || el_it != consensus_gps.end()) {

            // if a pointer is invalid, i.e, it's reached it's end, set it to nullptr
            // such that min_breakpoint_comp will never select it

            sr_it = (sr_it < sr_end) ? sr_it : consensus_gps.end();
            sl_it = (sl_it < sl_end) ? sl_it : consensus_gps.end();
            er_it = (er_it < er_end) ? er_it : consensus_gps.end();
            el_it = (el_it < el_end) ? el_it : consensus_gps.end();


            // get iterator to the gluepoint with the smallest pos value
            auto min_gp = std::min(
                std::min(sr_it, sl_it, min_breakpoint_comp),
                std::min(er_it, el_it, min_breakpoint_comp),
                min_breakpoint_comp
            );
            if (min_gp == consensus_gps.end()) {
                break;
            }

            // if there are any points less than max_sep away from the minimum:
            // 1) copy over the gluepoint number
            // 2) if the point isn't the min, zero out the gluepoint
            // 3) increment the iterator
            if (abs_diff(min_gp->pos, sr_it->pos) <= max_separation && sr_it < sl_end) {
                min_gp->sr = sr_it->sr;
                sr_it->sr = (min_gp != sr_it) ? 0 : sr_it->sr;
                sr_it++;
            }
            if (abs_diff(min_gp->pos, sl_it->pos) <= max_separation && sl_it < sl_end) {
                min_gp->sl = sl_it->sl;
                sl_it->sl = (min_gp != sl_it) ? 0 : sl_it->sl;
                sl_it++;
            }
            if (abs_diff(min_gp->pos, er_it->pos) <= max_separation && er_it < er_end) {
                min_gp->er = er_it->er;
                er_it->er = (min_gp != er_it) ? 0 : er_it->er;
                er_it++;
            }
            if (abs_diff(min_gp->pos, el_it->pos) <= max_separation && el_it < el_end) {
                min_gp->el = el_it->el;
                el_it->el = (min_gp != el_it) ? 0 : el_it->el;
                el_it++;
            }
        }

        // unique sections (i.e interval bounds with no alignments) need to be given a non-zero label to identify them.
        // We pass through the gluepoints, and for all gluepoints with some non-zero entry, we specify a unique number

        for (auto it = contig_start; it < contig_end; it++) {
            if (it->sr || it->sl || it->er || it->el) {
                if (it->sr == 0) {
                    it->sr = id_gen_sr;
                    id_gen_sr += 4;
                }
                if (it->sl == 0) {
                    it->sl = id_gen_sl;
                    id_gen_sl += 4;
                }
                if (it->er == 0) {
                    it->er = id_gen_er;
                    id_gen_er += 4;
                }
                if (it->el == 0) {
                    it->el = id_gen_el;
                    id_gen_el += 4;
                }
            }
        }

        std::stable_sort(
            contig_start, contig_end, [](Gluepoint const & a, Gluepoint const & b) {
                return a.pos < b.pos;
            }
        );
        // process next contig
        contig_start = contig_end;
    }


    // Remove empty gluepoints
    std::stable_sort(
        consensus_gps.begin(), consensus_gps.end(), [](Gluepoint const & a, Gluepoint const & b) {
            bool a_clear = !(a.sr || a.sl || a.er || a.el);
            bool b_clear = !(b.sr || b.sl || b.er || b.el);
            return a_clear < b_clear;
        }
    );
    auto true_back = consensus_gps.end()-1;
    while (true_back > consensus_gps.begin() && !(true_back->sr || true_back->sl || true_back->er || true_back->el)) {
        true_back--;
    }
    consensus_gps.erase(true_back+1, consensus_gps.end());



    return consensus_gps;
}

uint64_t minimize_point(uint64_t next_p, uint64_t min_p, std::unordered_multimap<uint64_t, uint64_t> const &id_eqv_pts,
                        std::unordered_set<uint64_t> &tracked,
                        std::unordered_map<uint64_t, uint64_t> const &min_lut)
{
    // already found the minimum of the graph in previous round.
    auto memo_min = min_lut.find(next_p);
    if (memo_min != min_lut.end())
    {
        return memo_min->second;
    }

    if (next_p < min_p)
    {
        min_p = next_p;
    }

    // if the point is unique, return it
    auto range = id_eqv_pts.equal_range(next_p);
    if (range.first == id_eqv_pts.end())
    {
        return min_p;
    }

    // otherwise, we need to basically dfs search
    for (auto it = range.first; it != range.second; it++)
    {

        if (tracked.find(it->second) != tracked.end())
        {
            // then we've already explored this point, so don't explore
            continue;
        }

        // track it
        tracked.insert(it->second);

        // explore it
        min_p = minimize_point(it->second, min_p, id_eqv_pts, tracked, min_lut);
    }
    return min_p;
}

void minimize_equivalent_points(std::vector<Gluepoint> &gps,
                                std::unordered_map<SetNode<Breakpoint> *, uint64_t> &set_to_id,
                                std::unordered_map<uint64_t, SetNode<Breakpoint> *> &id_to_set)
{
    for (auto it = gps.begin(); it < gps.end(); it++)
    {
        auto set_sr = id_to_set[it->sr];
        auto set_sl = id_to_set[it->sl];
        auto set_er = id_to_set[it->er];
        auto set_el = id_to_set[it->el];

        if (set_sr)
        {
            auto root_sr = findSetMerge(set_sr);
            it->sr = set_to_id[root_sr];
        }
        if (set_sl)
        {
            auto root_sl = findSetMerge(set_sl);
            it->sl = set_to_id[root_sl];
        }
        if (set_er)
        {
            auto root_er = findSetMerge(set_er);
            it->er = set_to_id[root_er];
        }
        if (set_el)
        {
            auto root_el = findSetMerge(set_el);
            it->el = set_to_id[root_el];
        }
    }

}
std::vector<SequenceInterval> build_sequence_intervals(std::vector<Gluepoint> const &gps)
{ // loop through gluepoint array and construct intervals
    std::vector<SequenceInterval> seq_ivls;
    seq_ivls.reserve(gps.size() - 1);

    for (auto it = gps.begin() + 1; it != gps.end(); it++)
    {
        if (it->cid != (it - 1)->cid)
        {
            continue;
        }        
        seq_ivls.emplace_back((it - 1)->sr, (it)->el, (it - 1)->er, it->sl, it->cid, (it - 1)->pos, it->pos);
    }
    return seq_ivls;
}

std::vector<SequenceInterval> cpp_gen_gp(std::vector<PyAlignment> &alignments, ContigContainerPtr contigs)
{
    auto max_separation = config.GetValue<unsigned int>("max_separation");
    std::cout << "construct endpoints from alignments\n";
    auto endpoints = construct_endpoints(alignments); // same output
    std::cout << "cluster endpoints: " << endpoints->size() << std::endl;

    auto clusters_q_projection = cluster_endpoints(*endpoints, max_separation, true); // same output     

    std::cout << "clusters_q_projection size: " << clusters_q_projection.size() << std::endl;
    std::unordered_map<uint32_t, std::vector<SetNode<Breakpoint> *>> bp_lookup;

    std::cout << "Grouping breakpoints into breakpoint sets..." << std::endl;
    initialize(*contigs, bp_lookup); // same output
    std::cout << "bp_lookup size: " << bp_lookup.size() << std::endl;

    std::unordered_set<SetNode<Breakpoint> *> equivalent_points;
    group_breakpoints_into_gluepoints(*endpoints, equivalent_points, clusters_q_projection, alignments, bp_lookup, max_separation, true);
    
    std::cout << "merging s-only clusters\n";
    std::cout << "equivalent_points size: " << equivalent_points.size() << std::endl;
    
    merge_s_only_clusters(*endpoints, bp_lookup, equivalent_points, max_separation);
    
    std::cout << "SIZE OF EQUIVALENT POINTS: " << equivalent_points.size() << std::endl;
    std::cout << "Building consensus gluepoints from breakpoint sets...\n";
    std::unordered_map<SetNode<Breakpoint> *, uint64_t> set_to_id;
    std::unordered_map<uint64_t, SetNode<Breakpoint> *> id_to_set;
    auto consensus_gluepoints = build_consensus_gluepoints(bp_lookup, set_to_id, id_to_set, max_separation);
    std::cout << "connecting roots through equivalent points...\n";
    connect_roots_through_equivalent_points(equivalent_points);

    std::cout << "Minimize equivalent points...\n";
    minimize_equivalent_points(consensus_gluepoints, set_to_id, id_to_set);

    std::cout << "Build sequence intervals...\n";
    auto seq_ivls = build_sequence_intervals(consensus_gluepoints);

    std::cout << "seq_ivls.size(): " << seq_ivls.size() << "\n";
    return seq_ivls;
}

#endif