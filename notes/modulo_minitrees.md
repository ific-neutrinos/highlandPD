# Módulo CreateHighlandMiniTree y minitrees ampliados (phase2.md)

Movido de CLAUDE.md el 23-nov-2026 (sin cambios de texto). Compilación en FNAL, grid y producción: notes/grid_produccion.md.

## Fase 2 de phase2.md: minitrees ampliados (en curso, 16-nov-2026)

Plan en `phase2.md` (raíz del repo, sin seguir por git a propósito). Ramas `phase2/minitrees` en highlandPD (desde `develop`) y en `../highlandAnalysis` (repo `ific-neutrinos/highlandAnalysis`, desde `master`, que va 6 commits por delante de su `develop`). Sin push. En highlandAnalysis, el tag `pre-phase2` es el `master` anterior (con `BeamPartIDEff`, `TracksForCali`, `.todo` y `setup.backup`, ya borrados); los tags del repo se nombran por release, y el siguiente sería `dunesw_v09_56_00d00` cuando compile en FNAL.

### Módulo (`highlandAnalysis/CreateHighlandMiniTree_module.cc`)

- **Aquí no se puede compilar**: no hay LArSoft. Las firmas se han comprobado con las cabeceras de dunesw v09_56_00d00 bajadas de GitHub (protoduneana y duneprototypes v09_56_00d00, larsim v09_24_03, larreco v09_12_08, larpandora v09_13_06, lardataobj v09_09_00, lardata v09_09_00, larcorealg v09_06_05, nusimdata v1_25_02). Hay que compilarlo en dunegpvm con esa versión. Revisado de principio a fin (commit `7c4c290`).
- `analyze()`: `BuildEventCaches` (productos, asociaciones, backtracking de todos los hits en MC), `FillEventInfo`, `FillDQInfo`, `Beam` (siempre), truth, bunch (beam slice), hits no asociados, cósmicos y cuentas de hits true.
- Eventos marcados, guardados para la normalización: DQ malo (`GoodDaq = false`, con `Beam` y sin bunch: los análisis no los ven y `CreateMiniTree` de highland los cuenta en el siguiente spill guardado), sin truth del haz, sin beam slice (`HasPandora = false`, bunch sin partículas de la beam slice). `Beam->SpillsSincePreviousSavedSpill = 1` y `POTSincePreviousSavedSpill = 0` en cada evento (antes -999).
- Hits: `Hits[3]` solo con calorimetría (como antes, para no cambiar los análisis); el resto de hits del PFParticle va a `HitsNoCalo[3]`. La calorimetría descarta sobre todo los hits 2D sin space point (comprobado en `LArPandoraTrackCreation` y `Calorimetry_module`). Pandora crea un space point por hit 2D, así que `AnaSpacePointPD` guarda un solo `HitID`.
- Flags (`CreateHighlandMiniTree.fcl`): `SaveHitsAllPlanes`, `SaveSpacePoints`, `SaveTrackSpacePoints` (falso), `CosmicsMode` 0-4 (cada nivel incluye el anterior: 1 resumen, 2 completos cerca de los extremos de la beam slice, 3 también si cruzan el volumen del haz, 4 todos; `Cosmics: true` antiguo = 4), `CosmicsRadius`, `BeamVolume{Start,Direction}{Data,MC}` (distintos en datos y MC; puntos de entrada de los cortes de calidad de 6 GeV, dirección nominal pendiente de medir), `BeamVolumeRadius`, `BeamVolumeLength`, `CosmicsTruthMode`, `SaveUnassociatedHits`, `UnassociatedHitsRadius`, `TrueTrajPDGs`, `SaveForcedTracks` (pandora2Track/pandora2cali para PFParticles sin traza), `SaveTrueHitCounts`, `SaveHitSignal`, `FullDumpPrescale`, `DataUtils` y `UnstableHVTimeRanges` (`@local::HVfilter.TimeRanges`).
- `AnaDataQualityPD`: lo que miran los filtros de FEMBs y HV inestable, sin filtrar. Los FEMBs y la consistencia de timestamps necesitan los `RawDigit` y sus `RDTimeStamp`, que no suelen estar en los ficheros reco: entonces quedan a -999.
- Decidido con el usuario: nada de Geant4Reweight por ahora.

### highlandPD

- Constructores de copia arreglados; miembros y clases nuevas (`AnaSpacePointPD`, `AnaTrueTrajPointPD`, `AnaForcedTrackPD`, `AnaDataQualityPD`, `HitsNoCalo`, etc.). Test de regresión idéntico a `fixA/` en cada paso.
- Protecciones para spills sin bunch o sin `Beam` en las correcciones, `GetBeamTrueParticle` y los `CreateMiniTree*` (acumulan la cuenta de spills antes de saltar).
- Limpieza de `pdDataClasses` (decidida con el usuario): fuera los miembros sin uso (posiciones y direcciones SCE de `AnaParticlePD`, `CALO`, las variables de PID obsoletas, `dEdx_calib` y los dE/dx y rangos intermedios de las calibraciones, `Pi0_decay_ID`, `LengthInTPC`, `MomentumInTPC`, `BeamMomentumInTPC`, `AnaWireCNN`/`CNNwires`), los que solo rellenaba el conversor de `pduneana` (`truncLibo_dEdx`, `RangeMomentum_alt`, `Length_alt`) y `vtx_CNN_michelscore`/`vtx_CNN_NHits` con sus 10 ramas de `kaonTree` (con minitrees del módulo valían siempre -999; se pueden recalcular con los hits y su CNN). Se quedan `isBeamPart` (lo lee `brokenTracksAnalysis`), `Distance_to_closest_particle` y `forced_daughter*` (los rellenan las selecciones). Duplicados para la Fase 4: `TPCid`/`PlaneID` frente a `WireID`, y `BeamMomentum` frente a `BeamParticle->Momentum`.
- `test/phase2/run_test_minitree.C`: minitree sintético con los miembros nuevos y eventos marcados; comprueba la ida y vuelta tras `Clone()`. Los análisis y los `CreateMiniTree*` corren sobre él, y `secondaryKaon` da la misma salida que sobre las entradas originales (solo cambia `entry`).

### Cómo funcionan los minitrees (comprobado)

- Cadena: módulo (`highlandana/MiniTree`, rama `Spill` de clase `AnaSpillPD`) → `CreateMiniTree*` de highlandPD (conversor `minitree` de `CreateMiniTreePD`) → minitree de análisis (rama `AnaSpill`, con árbol `config`). highland lee las dos con un `AnaSpill*`: ROOT acepta un puntero a la clase base, y dentro los objetos son `AnaBunchPD`, `AnaParticlePD`, etc.
- Tras `GetEntry` hay que llamar a `CopyVectorsIntoArrays()` y `RedoLinks()`, como `MiniTreeConverter::GetSpill`.
- `Clone()` de un spill **no clona la truth**: el clon comparte `TrueParticles` y `TrueVertices` con el original (`isClone = true`) y solo el original la borra. Hay que borrar el clon antes que el original. Antes el módulo nunca borraba el original y perdía la truth en cada evento.
- Las clases de `pdDataClasses` no tienen `ClassDef`: ROOT las versiona por checksum, así que añadir o quitar miembros no rompe la lectura de ficheros antiguos.


## Pendiente para la próxima versión del módulo (después de `d01`; lista abierta, 25-nov-2026)

Decisión del usuario (25-nov): no se para más la producción de MC. Se corre con `d01` una buena parte (al menos
`v1_00`, quizá la mitad de los bloques) para avanzar en la selección, y lo que salga de la revisión de la nota se
acumula aquí para una versión siguiente (`d02`) antes de los bloques que queden y de los datos.

- **Truth del pileup del haz completa (MC).** `FillBeamTrueInfo` solo guarda la good particle de `(*mcTruths)[0]` y su
  cadena; las demás partículas del haz del mismo `MCTruth` entran solo como ancestros de algo emparejado a un
  PFParticle guardado (`SaveTrueAncestors`). Falta: recorrer las primarias del `ParticleList` (madre 0) de ese
  `MCTruth` y guardarlas sin recursión si no están (Generation -999, Origin 4, como los ancestros; `KOrigin` no
  cambia). Flag propuesto `SaveBeamPileupPrimaries`. Coste estimado ~1 % (350 B por partícula). Da la tasa de pileup
  por trigger, cuántas entran en el TPC, su composición, momento y tiempo: el modelo del fondo de pileup de la medida B.
- **Instrumentación del haz con varias partículas (datos).** `FillBeamInfo` guarda `nTracks`, `nMomenta`, `nFibers[3]`,
  un TOF (`GetTOF`), los PDG compatibles y solo la traza 0 y el momento 0. Falta: todas las trazas del haz (inicio,
  final, dirección), todos los momentos, todos los TOF (y canales) y los tiempos de los Cherenkov si hay varios, para
  estudiar en datos los triggers con más de una partícula en la BI. Comprobar antes la ventana de lectura de la BI
  frente a la del TPC (el pileup del TPC está a -1/+3 ms del trigger; la BI, probablemente mucho más corta: por
  confirmar en la documentación de ProtoDUNE). En MC la BI se emula solo de la partícula que dispara: no hay
  equivalente.
