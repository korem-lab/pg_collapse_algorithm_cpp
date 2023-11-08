#ifndef SET_H
#define SET_H

#include <vector>
#include <unordered_map>

template <class T>
struct SetNode
{
    SetNode(T data): parent_merge(this), parent_ep(this), rank(0), data(data) {}

    SetNode* parent_merge;
    SetNode* parent_ep;
    uint64_t rank;
    T data;
};

template <class T>
inline SetNode<T>* findSetMerge(SetNode<T>* elem)
{
    if (elem->parent_merge != elem)
    {
        elem->parent_merge = findSetMerge(elem->parent_merge);
        return elem->parent_merge;
    }
    else
    {
        return elem;
    }
}

template <class T>
void unionSetMerge(SetNode<T>* const node1, SetNode<T>* const node2)
{
    SetNode<T>* root1 = findSetMerge(node1);
    SetNode<T>* root2 = findSetMerge(node2);
    if (root1 == root2) return;

    if (root1->rank > root2->rank)
    {
        root2->parent_merge = root1;
    }
    else
    {
        root1->parent_merge = root2;
        if (root1->rank == root2->rank)
        {
            ++root2->rank;
        }
    }
}
template <class T>
SetNode<T>* findSetEqvPt(SetNode<T>* elem)
{
    if (elem->parent_ep != elem)
    {
        elem->parent_ep = findSetEqvPt(elem->parent_ep);
        return elem->parent_ep;
    }
    else
    {
        return elem;
        return elem;
    }
}
template <class T>
void unionSetEqvPt(SetNode<T>* node1, SetNode<T>* node2)
{
    SetNode<T>* root1 = findSetEqvPt(node1);
    SetNode<T>* root2 = findSetEqvPt(node2);
    if (root1 == root2) return;
    if (root1 < root2)
    {
        root2->parent_ep = root1;
    }
    else
    {
        root1->parent_ep = root2;
    }
}

template <typename T>
std::unordered_map<SetNode<T>*, std::vector<T>>
groupBySet(const std::vector<SetNode<T>*>& sets)
{
    std::unordered_map<SetNode<T>*, std::vector<T>> groups;
    for (auto& setNode : sets)
    {
        groups[findSetMerge(setNode)].push_back(setNode->data);
    }
    return groups;
}

//Vector that stores the set nodes and automatically deletes
//them in the end. Does not have virtual table - do not use polymorphism!
//All set elements should be pushed before any set operations are
//applied (like union) - do not modify the vector afterwards
template <typename T>
class SetVec : public std::vector<SetNode<T>*>
{
public:
    ~SetVec()
    {
        for (auto& x : *this) delete x;
    }
};

#endif 