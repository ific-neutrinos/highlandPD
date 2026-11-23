# Fase 4: reestructuración de highlandPD (después de reestructurar highland)

Movido de CLAUDE.md el 23-nov-2026 (sin cambios de texto).


- **Análisis con varios candidatos por evento** (dice el usuario, 23-nov): highland trabaja por evento, y `secondaryKaon` necesitó apaños (ramas por candidato, `FinalizeConfiguration` propio para el `accum_level` y el árbol `truth`). El hotfix de `accum_level` de `secondaryProton`/`secondaryKaonProductionXS` es temporal: en la reestructuración, repensar cómo highland trata selecciones por candidato y sus eficiencias en la truth.
- Separar el núcleo `pd*` (una librería) de los análisis (una librería y un ejecutable cada uno, opcionales vía CMake).
- **Convención de nombres de las clases** (decidido con el usuario el 20-nov): unificar `xxxPD` / `pdXxx` en las clases internas (`EventBoxPD`, `ToyBoxPD`, `CreateMiniTreePD` frente a `pdBaseAnalysis`, `pdBaseConverter`, `pdAnaUtils`...). Nunca las clases de datos (`AnaSpillPD`, `AnaParticlePD`...): ROOT guarda sus nombres en los minitrees. Los paquetes, las clases de análisis y selección, los ficheros de parámetros y los ejecutables ya siguen la convención de highland (paquetes y clases de análisis en minúscula, como `baseAnalysis`).
- Factorizar el copy-paste:
  - El bloque que recalcula los χ² (`Chi2Proton`, `Chi2Muon`, `Chi2ndf`) está en 4 variaciones y su restauración en 3 `UndoSystematic`: helper común (sirve también para el bug A1 de la auditoría).
  - El bucle "partículas + `Original`" se repite en las 11 correcciones y el de la `SystBox` en las 6 variaciones.
  - `ParticlePositionSCECorrection` repite la parte geométrica de `SCEGeometricVariation`; `SCECorrection` = `HitPositionSCECorrection` + `HitPitchSCECorrection`.
  - Los filtros de `CreateMiniTreeCandidates` y `CreateMiniTreeEventCandidates` son iguales y con números cableados.
  - `secondaryKaonAnalysis` vs `secondaryKaonProductionXSAnalysis`: 73 % iguales.
  - `CreateMiniTreeCandidates` vs `CreateMiniTreeEventCandidates`: 76 % iguales.
- `pdBaseAnalysis` debería registrar las correcciones y sistemáticos estándar, activables por parámetros. `pdExampleAnalysis` no usa sistemáticos.
- `EventBoxPD::FillLongTracks` usa `anaUtils::GetAllTracksUsingDet`, que hace `static_cast<AnaTrackB*>` sobre `AnaParticlePD`. Ver el `CLAUDE.md` de highland.
- `pdDataClasses` (ojo: `phase2.md` decide no quitar nada de las clases por ahora, y necesita `EmptyEvent`, `HasPandora` y `BeamMomentum`):
  - Miembros marcados como prueba u obsoletos: `Distance_to_closest_particle` (l.224), `forced_daughter`/`forced_daughter_matched` (l.241-242), `AnaHitPD::dEdx_calib`.
  - `EmptyEvent` y `HasPandora` "should probably go in DataQuality".
  - `CNNwires` está en `AnaBunchPD` y en `AnaEventPD`.
  - Las posiciones SCE de `AnaParticlePD` son `double`; el resto de la jerarquía usa `Float_t`.
- `AnaBeamPD::BeamMomentum` duplica `BeamParticle->Momentum`: `BeamMomSmearingCorrection` iguala uno al otro.
- Del `.todo` de los autores (borrado el 16-nov), sigue vigente:
  - Las selecciones no funcionan si se añaden splits antes de un primer corte (highland); documentarlo o arreglarlo.
  - Los índices de cortes de las cabeceras de las selecciones no se usan.
  - Constantes de `pdDataClasses.hxx` que solo usa `standardPDTree`: moverlas a su namespace.
  - Los constructores de los análisis PD y de `baseAnalysis` no hacen nada; ¿añadir `DefineCategories` a `AnalysisAlgorithm`?
  - Unificar los nombres `xxxPD` / `pdXxx`. ¿Por qué `pdAnaUtils` y no extender `anaUtils`?
  - `ToyBoxPD::Reset` repite lo de `ToyBoxB`. `EventBoxPD` con varios grupos de objetos es difícil de usar.
  - Cortes de calidad del haz por energía y run; versionado (tag por presentación); dependencias entre paquetes y ficheros de parámetros.
  - Los diccionarios salen en `build/` y no en `build/dict` (CMake no acepta `/` en los nombres).
- Bug de highland (`psycheEventModel/src/BaseDataClasses.cxx:1583`): `ConvertProcess` espera `"kaon.Inelastic"` en vez de `"kaon-Inelastic"`, así que el inelástico de K- no se reconoce. Afecta a `phase2.md` (cadena del K±).
- Coordinación con highland: `EventBoxId` de highland define `kEventBoxPD`, `kEventBoxKaon` y `kEventBoxdEdx`. Cuando highland permita registrar IDs desde fuera, moverlos aquí.

