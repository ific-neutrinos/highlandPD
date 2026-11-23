# Servicios (SCE, calorimetría, PID) y ficheros de datos

Movido de CLAUDE.md el 23-nov-2026 (sin cambios de texto).

## Servicios de calorimetría y SCE (hecho, 13-nov-2026)

Decidido con el usuario (primero el bloque de diseño, antes que los bugs sueltos):

- **`pdServiceManager`** en `pdUtils`, con acceso global `pd::services()` al estilo de `ND::params()`. Tiene que ser global porque las correcciones y los sistemáticos no tienen acceso a `ana()`. Usa métodos con tipo (`SCE()`, `Calo()`, `PID()`), creación perezosa (nunca en la inicialización estática), objetos `const` y lista de servicios en el log.
- **Auditoría:** todas las `Calorimetry` y `SpaceCharge` tienen la misma configuración base. Solo cambia estado por evento o por toy: lifetime (`LifetimeVariation`, `SCEVariation`), ModBox A/B (`RecombinationVariation`) y qué SCE usar. Basta con una instancia nominal de cada una, y ese estado viaja como contexto en cada llamada.
- **SCE en los toys sin copias del mapa.** La interpolación (spline natural `"b2e2"` en un eje y bilineal en los otros) es lineal en los valores de los nodos. Variación global: `D_toy = f·D_nominal`. Por vóxel: funciones base del spline (`B_i(x)`, calculadas una vez por eje) por los factores del toy, que ya guarda `ToyExperiment`. Hay que mantener los dos modos (por vóxel era el objetivo; nunca se midieron incertidumbres por vóxel, por eso se usa la global). El sistemático de SCE cuesta hoy ~2 GB y ~10 s (100 mapas de ~20 MB).
- **Validación:** el servicio debe dar salida idéntica bit a bit. El SCE global por factor no saldrá idéntico (hoy el mapa escalado se guarda en `TH3F` float); se acepta con tolerancia, enseñando antes al usuario el tamaño de las diferencias.
- Pasos: (1) servicio con SCE y Calorimetry nominales + contexto, migrando las clases vivas; (2) SCE global por factor; (3) SCE por vóxel; (4) unificar `SCEVariation` y `SCEGeometricVariation`.
- **Paso 2 hecho** (commit `5ae2d90`): `SCEScaling` (por ahora solo `Global`) como argumento opcional de las consultas de `SpaceCharge`, y `CaloContext::SCEScale`. Las variaciones calculan el factor con `GetToyScaling(toy)` y ya no hay mapas por toy ni límite de 100 toys. Memoria 2,44 GB -> 364 MB, inicialización ~9 s -> ~1 s. Se quitó la API que modificaba el mapa y la variación por vóxel antigua (`VarySCEMapLocally`, no usada, que necesitaba una copia del mapa por toy); el paso 3 la reimplementa sin copias.
- **Paso 1 hecho** (commit `9f8d1d9`): `pdServiceManager` con `pd::services().SCE()` y `pd::services().Calo()`. El método se llama `Calo()` y no `Calorimetry()` porque un método con el nombre de la clase la oculta dentro de `pdServiceManager`. `Calorimetry` es inmutable, recibe un `CaloContext` (`NominalContext()` + cambios) y lee sus constantes de `pdUtils.parameters.dat`. Test idéntico. Las clases que solo usan las calibraciones (`SCECorrection`, `HitPositionSCECorrection`, `HitPitchSCECorrection`, `LifetimeCorrection`, `LifetimeVariation`, `SCEVariation`, `dQdxXCalibration`) solo se han validado compilando: no hay minitree para `dEdxCalibration.exe`.
- **Pasos 3 y 4 hechos** (commits `1ae7c28`, `56dbba0`):
  - `SCEVariationBase` reúne lo común (mapa del servicio y `GetToyScaling`). Las dos derivadas mantienen efectos distintos a propósito: `SCEGeometricVariation` (kaones) solo cambia la geometría, porque la calorimetría la cubre `dEdxVariation`; `SCECalorimetricVariation` (antes `SCEVariation`, calibraciones) cambia las posiciones de los hits y el pitch de dQ/dx. No unificarlas.
  - Modo con `pdSystematics.SCE.VariationMode` en `src/pdSystematics/parameters/pdSystematics.parameters.dat`: 0 global (`data/SCE.dat`, sigma 8 %), 1 por vóxel (`data/SCE2.dat`, 21.793 cajas con sigma 5 % de prueba). El modo por vóxel usa un lector propio porque `BinnedParams` tiene `NMAXSYSTSOURCEBINS = 50`.
  - Por vóxel, `SpaceCharge` evalúa el spline con la base cardinal (`_basis`) y los factores de cada nodo, sin copias. Validado con `root -l -b -q test/sce/run_test_sce_voxel.C` frente al método antiguo (copia escalada y splines nuevos): diferencias relativas <= 4e-8.
  - Coste del modo por vóxel en `secondaryKaon.exe` (1000 eventos, 10 toys): 52 s frente a 22 s en global; memoria igual.
  - `SCE2.dat` tiene cajas independientes para los dos volúmenes (42.439); los nodos de x = 0 (cátodo) comparten caja. La implementación antigua repartía los factores por índice de bin entre volúmenes, casi seguro sin querer.
- Las variables de toys solo se rellenan en eventos con super candidato (el corte 3 borra los candidatos si no lo hay): por eso el test usa todo el MC.
- Hecho: quitadas las 6 clases muertas (`CalorimetryCalibration`, `RecombinationMCCorrection`, `dQdx{Cal,Norm,XCal,YZCal}Variation`) y sus `.dat`. Copia en el backup.

Encontrado en la auditoría, pendiente:

- `LifetimeCorrection` (solo en las calibraciones) usa siempre la lifetime de MC. Según el usuario, las calibraciones no se corren sobre datos (solo hay un rastro comentado en `dQdxXCalibration.cxx:90`), así que no se cambia.
- Campo E: `Calorimetry` usa 0,4867 kV/cm y `RecombinationDataCorrection` usa 0,553 (con sus propios ModBox). El usuario no recuerda si es a propósito. Al pasarlos a parámetros, mantener los dos valores con nombres distintos para no cambiar resultados.
- `dQdxXCalibration` lee `$HIGHLANDPDPATH/yz_corrections_mc_lifetime_mypc.root`, un fichero personal que no está en el repo.

## Servicios y ficheros de datos

- `pd::services()` (`pdUtils/src/pdServiceManager`): `SCE()`, `Calo()` y `PID()` (`PIDTables`: plantillas de dE/dx por PDG, energía cinética frente a rango por partícula y tablas NIST de muón y protón; un único servicio porque todo depende de la hipótesis de partícula). Se crean la primera vez que se piden, también dentro del bucle de eventos; el guardia restaura `gDirectory` y evita que los histogramas se cuelguen del directorio actual.
- `pd::DataFile(paquete, fichero)` y `pd::DataDir(paquete)` (`pdUtils/src/pdDataFiles`): único punto para las rutas de datos. `$<PKG>ROOT/data`, o antes `$HIGHLANDPD_DATADIR/<paquete>` si está definida (para tener los datos fuera del repo). Error claro y `exit(1)` si falta algo. No usar `getenv` para rutas de datos.
- Nombres de fichero en `pdUtils.parameters.dat`: `pdUtils.SCE.MapFile`, `pdUtils.Calorimetry.{Norm,XCal,YZCal}File`, `pdUtils.PID.{dEdxTemplatesFile,KEvsRangeFile}`.
- En highland y psyche no hay nada equivalente (solo los singletons `ND::params()` y `ND::versioning()`, y `BinnedParams`). `BinnedParams` hace `exit(0)` si no encuentra su fichero: anotar para la reestructuración de highland. El gestor de servicios y `pd::DataFile` podrían subir a highland como algo genérico.
- Masas: constantes `pdAnaUtils::k{Muon,Pion,Kaon,Proton}Mass` (PDG 2022) y `pdAnaUtils::GetParticleHypothesis(pdg, nombre, masa)`. Las funciones que dependen de la hipótesis avisan por `std::cerr` y devuelven -999 (o -1 en `ComputeCSDARange`) si el PDG no está soportado; no usar valores locales.
- Pendiente: `dQdxXCalibration` lee `$HIGHLANDPDPATH/yz_corrections_mc_lifetime_mypc.root`, que no está en el repo. Las macros de `dEdxCalibrationAnalysis` tienen rutas personales de Fermilab.

