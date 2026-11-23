#ifndef ToyBoxMuCand_h
#define ToyBoxMuCand_h

#include "ToyBoxPD.hxx"
#include "muCandidateUtils.hxx"

/// Caja de la selección de candidatos a mu+ de K+ en reposo (secondaryKaonProductionXSSelection): los candidatos con
/// sus variables (muCandUtils::Candidate) y, para cada uno, hasta qué corte de la cadena de su población ha llegado.
class ToyBoxMuCand: public ToyBoxPD{
public:

  ToyBoxMuCand();
  virtual ~ToyBoxMuCand(){}

  /// Obligatorio en las clases derivadas (no hace nada; ver ToyBoxB)
  virtual void Reset();

  /// Vacía los candidatos
  virtual void ResetBase();

public:

  /// Candidatos que pasan la preselección (ver GetMuCandidatesAction)
  std::vector<muCandUtils::Candidate> Candidates;

  /// Para cada candidato: población (0 en la beam slice, 1 fuera) y número de cortes de la cadena de su población
  /// que pasa seguidos (0: ninguno; el primero es el de la población)
  std::vector<Int_t> Population;
  std::vector<Int_t> Level;

  /// Candidatos antes de la preselección geométrica del módulo
  Int_t NPreselected;

  /// Vértice usado por las variables de conexión: final del primario de la beam slice (sin él, HasVertex = false)
  bool HasVertex;
  Float_t VertexPos[3];
};

#endif
