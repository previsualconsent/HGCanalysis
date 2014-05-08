#include "UserCode/HGCanalysis/interface/HuffmanAlgo.h"

#include <bitset>

using namespace std;

//builds the Huffman tree and returns it's root node
HuffmanTreeNode* BuildHuffmanTree(TH1 *h)
{
  if(h==0) return 0;
  
  //the first element is the one with highest priority
  std::priority_queue<HuffmanTreeNode*, std::vector<HuffmanTreeNode*>, HuffmanTreeNodeCmp> trees;
  
  //build the frequency ordered tree
  for(int i = 1; i <= h->GetXaxis()->GetNbins(); i++)
    {
      int cts=h->GetBinContent(i);
      if(cts) trees.push( new HuffmanTreeNode(i,cts) );
    }
  
  //build the Huffman tree setting parent -> (child1, child2) nodes
  while (trees.size() > 1)
    {
      HuffmanTreeNode* childR = trees.top();
      trees.pop();
      
      HuffmanTreeNode* childL = trees.top();
      trees.pop();
      
      HuffmanTreeNode* parent = new HuffmanTreeNode(childR, childL);
      trees.push(parent);
    }
  
  //all done here
  return trees.top();
}


//
void GenerateHuffmanCodes(const HuffmanTreeNode* node, const HuffmanCode prefix, HuffmanCodeMap& outHuffmanCodes)
{
  if( node->internal )
    {
      HuffmanCode leftPrefix(prefix.val<<1,prefix.nbits+1);
      GenerateHuffmanCodes(node->left, leftPrefix, outHuffmanCodes);
      
      HuffmanCode rightPrefix(prefix.val<<1 | 1,prefix.nbits+1);
      GenerateHuffmanCodes(node->right, rightPrefix, outHuffmanCodes);
    }
  else
    {
      outHuffmanCodes[node->c] = prefix;
    }
}


//
void testHuffmanAlgo()
{

  //historam a string
  //char SampleString[] = "this is an example for huffman encoding";
  char SampleString[] = "aaaaaabbbbbbccccdddeeeffghhhhiiiiiijjjkkllmn";
  std::map<char,int> strCts;
  for(size_t i=0; i<sizeof(SampleString)-1; i++)
    {
      if( strCts.find( SampleString[i] ) == strCts.end() ) strCts[ SampleString[i] ]=0;
      strCts[ SampleString[i] ]++;
    }

  TH1F *h=new TH1F("hstr",TString("String histo for:")+SampleString,strCts.size(),0,strCts.size());
  int xbin(1);
  for(std::map<char,int>::iterator it = strCts.begin(); it!= strCts.end(); it++, xbin++)
    {
      h->GetXaxis()->SetBinLabel(xbin,TString(it->first));
      h->SetBinContent(xbin,it->second);
    }
  h->Draw();

  //encode
  HuffmanTreeNode *rootNode = BuildHuffmanTree(h);
  HuffmanCodeMap codes;
  GenerateHuffmanCodes(rootNode, HuffmanCode(), codes);
  delete rootNode;

  for(HuffmanCodeMap::iterator it = codes.begin();
      it!=codes.end();
      it++)
    {
      std::bitset<8> binVal(it->second.val);
      cout << it->first << " " << h->GetXaxis()->GetBinLabel(it->first) << " ";
      for(unsigned int ibit=0; ibit<it->second.nbits; ibit++)
	cout << binVal[it->second.nbits-ibit-1];
      cout << endl;
    }
  
}


