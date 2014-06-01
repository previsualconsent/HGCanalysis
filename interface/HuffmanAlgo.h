#ifndef _huffman_algo_
#define _huffman_algo_

#include "TH1F.h"

#include <iostream>
#include <queue>
#include <map>
#include <iterator>
#include <algorithm>


/**
   @short A node in the Huffman tree containing, frequency, code and neighbours
   examples can be found in xhttp://rosettacode.org/wiki/Huffman_coding
*/
class HuffmanTreeNode 
{
 public:
  
    HuffmanTreeNode(int cval, int fval) : c(cval), f(fval), left(0), right(0), internal(false) { }
    HuffmanTreeNode(HuffmanTreeNode* c0, HuffmanTreeNode* c1) : c(0), f(c0->f + c1->f), left(c0), right(c1), internal(true) {}
    ~HuffmanTreeNode() { if(left) delete left; if(right) delete right; }
    void print() { std::cout << "c=" << c << " f=" << f << std::endl; }
    const int c;    
    const int f;
    HuffmanTreeNode *const left;
    HuffmanTreeNode *const right;
    bool internal;
};

/**
   @short returns the node with higher frequency
*/
struct HuffmanTreeNodeCmp
{
  bool operator()(const HuffmanTreeNode* lhs, const HuffmanTreeNode* rhs) const { return lhs->f > rhs->f; }
};


/**
   @short a value with the number of bits
 */
class HuffmanCode 
{
 public:
    HuffmanCode() : val(0), nbits(0) {}
    HuffmanCode(unsigned int v,unsigned int n) : val(v), nbits(n) { }
    unsigned int val, nbits;
};
typedef std::map<int, HuffmanCode> HuffmanCodeMap;


HuffmanTreeNode* BuildHuffmanTree(TH1 *h,bool checkNullEntries);
void GenerateHuffmanCodes(const HuffmanTreeNode* node, const HuffmanCode prefix, HuffmanCodeMap& outCodes);
HuffmanCodeMap getHuffmanCodesFrom(TH1F *h,bool checkNullEntries=true);
void testHuffmanAlgo();



#endif
