#include "UserCode/HGCanalysis/interface/HuffmanAlgo.h"

#include <bitset>

using namespace std;

//builds the Huffman tree and returns it's root node
HuffmanTreeNode* BuildHuffmanTree(TH1 *h,bool checkNullEntries)
{
  if(h==0) return 0;
  
  //the first element is the one with highest priority
  std::priority_queue<HuffmanTreeNode*, std::vector<HuffmanTreeNode*>, HuffmanTreeNodeCmp> trees;
  
  //build the frequency ordered tree
  for(int i = 1; i <= h->GetXaxis()->GetNbins(); i++)
    {
      int cts=h->GetBinContent(i);
      if(checkNullEntries && cts==0) cts=1; 
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
HuffmanCodeMap getHuffmanCodesFrom(TH1F *h, bool checkNullEntries)
{
  //encode
  HuffmanTreeNode *rootNode = BuildHuffmanTree(h,checkNullEntries);
  HuffmanCodeMap codes;
  GenerateHuffmanCodes(rootNode, HuffmanCode(), codes);
  delete rootNode;

  for(HuffmanCodeMap::iterator it = codes.begin();
      it!=codes.end();
      it++)
    {
      std::bitset<8> binVal(it->second.val);
      cout << it->first << " " << h->GetXaxis()->GetBinLowEdge(it->first) << " " << h->GetBinContent(it->first) << " ";
      for(unsigned int ibit=0; ibit<it->second.nbits; ibit++) cout << binVal[it->second.nbits-ibit-1];
      cout << " " << it->second.nbits << endl;
    }

  return codes;
}

//
int getTriggerBits(float nMips,float eta)
{
  int nBits(2);
  if( fabs(eta)<2.0 )
    {
      if(nMips<10)      nBits=1;
      else if(nMips<96) nBits+=4;
      else              nBits+=8;
    }
  else if(fabs(eta)<2.5)
    {
      if(nMips<10)      nBits=1;
      else if(nMips<96) nBits+=4;
      else              nBits+=8;
    }      
  else 
    {
      if(nMips<25)       nBits=1;
      else if(nMips<192) nBits+=4;
      else               nBits+=8;
    }      
  return nBits;
}

//
int getReadoutBits(float nMips,float eta)
{
  int nBits(2);
  if( fabs(eta)<2.0 )
    {
      if(nMips<0.4)      nBits=1;
      else if(nMips<6.4) nBits+=6;
      else               nBits+=10;
    }
  else if(fabs(eta)<2.5)
    {
      if(nMips<0.4)      nBits=1;
      else if(nMips<6.4) nBits+=6;
      else               nBits+=10;
    }      
  else 
    {
      if(nMips<0.4)      nBits=1;
      else if(nMips<6.4) nBits+=6;
      else               nBits+=10;
    }      
  return nBits;
}

//
std::map<string,float> testCompressionAlgos(TH1F *readoutH,TH1F *triggerH,float eta)
{
  std::map<std::string,float> results;

  std::map<std::string,TH1F *> histos;
  histos["readout"]=readoutH;
  histos["trigger"]=triggerH;
  for(std::map<std::string,TH1F *>:: iterator it=histos.begin(); it!= histos.end(); it++)
    {
      TH1F *h=it->second;
      bool isTrigger( it->first=="trigger" );

      //how many bits per word needed to fully describe this spectrum
      int nbitsPerWord=TMath::Log(h->GetXaxis()->GetNbins())/TMath::Log(2)+1;
      cout << nbitsPerWord << endl;

      //get the Huffman codes
      HuffmanCodeMap codes=getHuffmanCodesFrom(h);

      //sum up
      int totalBits(0), totalHuffmanBits(0), totalHGCBits(0);
      for(HuffmanCodeMap::iterator it = codes.begin();
	  it!=codes.end();
	  it++)
	{
	  float thr=h->GetXaxis()->GetBinLowEdge(it->first);
	  float cts=h->GetBinContent(it->first);
	  totalBits        += cts*nbitsPerWord;
	  totalHuffmanBits += cts*it->second.nbits;
	  totalHGCBits     += cts*(isTrigger ? getTriggerBits(thr,eta) : getReadoutBits(thr,eta));
	}


      results["Huffman " + it->first] = ((float)totalHuffmanBits/(float)totalBits);
      results["HGC " + it->first]     = ((float)totalHGCBits/(float)totalBits);
    }

  //all done here
  return results;
}

