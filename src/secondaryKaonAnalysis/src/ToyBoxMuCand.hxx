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

  /// Eventos con volcado completo (FullDump, MC): cada candidato está en una de las dos vistas del evento. 0 la vista
  /// normal (lo que el módulo guardaría en un evento normal, muCandUtils::MakeNormalView): la de la selección; 1 la
  /// vista completa (todo lo guardado): pasa por los mismos cortes, pero no cuenta para que el evento pase ni se
  /// guarda en cand_*; sirve para la eficiencia con toda la información (truemu_level_full). En los eventos normales
  /// solo hay vista 0
  std::vector<Int_t> View;

  /// Cómo guardaría el módulo el candidato en un evento normal (muCandUtils::NormalSave: 2 completo, 1 resumen)
  std::vector<Int_t> SaveMode;

  /// Evento con volcado completo
  bool FullDump;

  /// Copias en resumen de la vista normal (las partículas de los candidatos de vista 0 que no son del evento)
  std::vector<std::shared_ptr<AnaParticlePD>> SummaryCopies;

  /// Candidatos antes de la preselección geométrica del módulo
  Int_t NPreselected;

  /// Vértice usado por las variables de conexión: final del primario de la beam slice (sin él, HasVertex = false)
  bool HasVertex;
  Float_t VertexPos[3];

  /// MC: vértice perfecto en el sistema de la reconstrucción (último hit de la partícula true del haz; ver
  /// muCandUtils::EventContext::tvtx)
  bool HasTrueVertex;
  Float_t TrueVertexPos[3];
};

#endif
