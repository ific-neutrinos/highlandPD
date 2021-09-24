#ifndef cosmicsCreateMiniTree_h
#define cosmicsCreateMiniTree_h

#include "SimpleLoopBase.hxx"
#include "pdDataClasses.hxx"
#include "InputManager.hxx"
#include <set>

class cosmicsCreateMiniTree: public SimpleLoopBase {
 public:

  cosmicsCreateMiniTree(int argc, char *argv[]);
  virtual ~cosmicsCreateMiniTree(){}

protected:

  //---- These are mandatory functions
  virtual bool Initialize     ();
  virtual bool InitializeSpill(){return true;}

  virtual void DefineOutputTree();

  virtual void FinalizeSpill(){}
  virtual void Finalize     ();

  virtual bool Process();
  //--------------------

  virtual void FillMiniTree();

  bool RecoCandidateExists();
  bool TruthCandidateExists();

  enum miniTreeIndex{
    minitree = OutputManager::enumSpecialTreesLast+1
  };


protected:

  Int_t _totalEntries;
  Int_t _savedEntries;

  bool _saveCandidates;

  AnaSpill* _spill;
};

#endif