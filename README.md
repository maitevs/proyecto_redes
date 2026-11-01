# Atribución especie-resuelta de biotransformaciones multi-paso en comunidades microbianas

Proyecto semestral — optimización de redes metabólicas
Pontificia Universidad Católica de Chile, DIQB

---

## Problema

Los modelos de comunidad basados en AGORA2 predicen si un microbioma metaboliza
un fármaco, pero entregan un resultado agregado: sí o no a nivel comunitario. No
identifican qué especie ejecuta cada paso de una cascada de biotransformación ni
qué metabolito intermedio se traspasa entre miembros.

El abordaje actual (FBA sobre el modelo de comunidad, complementado con shadow
prices) depende de una solución óptima particular, y la atribución que produce
no es independiente de esa solución.

## Objetivo

Desarrollar un método de selección de comunidades mínimas **ancladas a un
miembro fijo**, capaz de identificar qué organismos son estructuralmente
necesarios para una biotransformación multi-paso y cuáles son funcionalmente
redundantes, con validación posterior de viabilidad de crecimiento.

El método es independiente del caso: funciona para cualquier organismo ancla,
cualquier función metabólica objetivo y cualquier comunidad de partida. La
biotransformación de fármacos es el caso de prueba, no el método.

## Aporte metodológico

`mincom` (Metage2Metabo) busca cualquier comunidad mínima que produzca los
targets. Aquí se impone la restricción de que un organismo determinado esté
presente en toda solución, cambiando la pregunta de *"¿cuál es la comunidad
mínima?"* a *"¿cuál es la comunidad mínima alrededor de este organismo?"*.

---

## Pipeline

| Etapa | Herramienta | Entrada | Salida |
|---|---|---|---|
| Modelos | AGORA2 v2.01 | — | 15 GSMNs (SBML) |
| **Preprocesamiento** | `preparar_sbml_para_m2m.py` | GSMNs | GSMNs sin reacciones de frontera |
| Validación | Seed2LP | GSMN + target | conjunto mínimo de seeds |
| Medio | derivado de los modelos | GSMNs | lista de seeds |
| Selección cualitativa | Metage2Metabo (`mincom`) | seeds + targets + pool | simbiontes esenciales / alternativos |
| Validación cuantitativa | SteadyCom | comunidad seleccionada | viabilidad de crecimiento |

La etapa de preprocesamiento no estaba en el diseño original. Se incorporó tras
el hallazgo descrito más abajo.

---

## Caso de estudio

Cascada de levodopa, documentada experimentalmente (Maini Rekdal et al., 2019):

```
levodopa  --TyrDC-->  dopamina  --Dadh-->  m-tiramina
       (E. faecalis)          (E. lenta)
```

La conversión reduce la biodisponibilidad del fármaco en pacientes con
Parkinson. Ningún organismo conocido ejecuta ambos pasos, lo que la hace un
control adecuado para un método de atribución.

**Organismo ancla:** *Enterococcus faecalis* ATCC 29212 (paso 1).

Nota: ATCC 29212 es una cepa de referencia de laboratorio, no una cepa
probiótica. El método requiere que el ancla sea un miembro fijo, no que sea
probiótico. Una extensión natural es anclar en *E. faecalis* Symbioflor 1
(DSM 16431), que sí tiene uso probiótico documentado.

**Pool:** 15 cepas, con redundancia intencional en ambos pasos de la cascada
(3 *E. lenta*, 2 *Enterococcus*) para evaluar la capacidad del método de
distinguir simbiontes esenciales de alternativos, más 9 comensales intestinales
como control negativo.

---

## Resultados

### 1. Verificación de la cascada en AGORA2

Se descargaron 15 reconstrucciones de AGORA2 v2.01 (release 2026-09-18) y se
inspeccionó la presencia de los metabolitos de la cascada.

**IDs encontrados:** `levdopa`, `dpmn` (dopamina), `mtym` (m-tiramina), en
compartimentos `_c`, `_e` y `_p`.

| Organismo | levdopa | dpmn | mtym | Rol |
|---|:---:|:---:|:---:|---|
| *E. faecalis* ATCC 29212 | sí | sí | — | ancla, paso 1 |
| *E. faecium* 504 | sí | sí | — | redundancia paso 1 |
| *E. lenta* A2 | sí | sí | **sí** | paso 2 |
| *E. lenta* DSM 11767 | sí | sí | **sí** | paso 2 |
| *E. lenta* DSM 2243 | — | sí | **sí** | paso 2 |
| *E. coli* 042 | sí | sí | — | control impuro |
| Resto de controles (8) | — | — | — | control |

La cascada está anotada y distribuida entre organismos. *E. lenta* DSM 2243
carece de levodopa pero conserva dopamina y m-tiramina: fenotipo de
dehidroxilador puro, coherente con la observación experimental de que solo una
fracción de las cepas de *Eggerthella* convierte dopamina cuantitativamente.

Los controles contienen `dopa4glcur_e` y `dopa3glcur_e` (glucurónidos de
dopamina, sustratos de beta-glucuronidasa), que no pertenecen a la cascada.

### 2. Validación del método con Seed2LP

Se inquirió a *E. lenta* A2 qué nutrientes requiere para exportar m-tiramina,
bajo tres configuraciones:

| Dopamina | Modo | Objetivo | Resultado | Tiempo |
|---|---|---|---|---|
| prohibida | hybrid | m-tiramina + biomasa | sin solución | 2714 s |
| permitida | reasoning | m-tiramina | 4 seeds (incl. `M_dpmn_e`) | 5,7 s |
| prohibida | reasoning | m-tiramina | 5 seeds (incl. `M_levdopa_c`) | 4,8 s |

Seed2LP identificó dopamina como seed mínimo sin información previa,
recuperando la dependencia interespecie descrita experimentalmente.

El modo `hybrid` impone adicionalmente flujo en la reacción de biomasa
(`R_bio1`), lo que explica la inviabilidad de la primera configuración y su
costo computacional (~45 min frente a ~5 s en `reasoning`).

### 3. Hallazgo: las reacciones de frontera invalidan el seed en m2m

La primera ejecución de `mincom` devolvió una comunidad mínima de **un solo
organismo** (*E. lenta* DSM 2243), sin intercambios (`exchanged: []`) y sin
simbiontes esenciales. *E. faecalis* no aparecía.

La inspección de las redes descartó las explicaciones iniciales:

- *E. lenta* A2 contiene una sola reacción con levodopa, `3HLYTCL`
  (`levdopa_c → dpmn_c + co2_c`), y **ningún transportador**: no puede
  absorberla del medio.
- *E. lenta* DSM 2243 ni siquiera posee `3HLYTCL`. Su ruta es
  `DOPAENT4tc` (importa dopamina) → `DOPADH` → `MTYMt2r` (exporta m-tiramina).
  Requiere dopamina extracelular, que no estaba en el seed.

**Control decisivo.** Se ejecutó `mincom` con un seed reducido a agua, protón y
fosfato. m-tiramina resultó producible — estequiométricamente imposible, al no
haber fuente de carbono ni de nitrógeno.

**Causa.** En formato COBRA una reacción de intercambio se escribe sin
productos:

```
EX_dopa_e_ | dpmn_e <=>
```

Network expansion activa una reacción cuando todos sus sustratos están
disponibles. En la dirección inversa esta reacción no tiene sustratos, de modo
que la condición se cumple vacuamente y el metabolito queda disponible con
independencia del seed.

Seed2LP elimina estas reacciones por defecto (`Import reaction: Removed`, con
flag `--keep-import-reactions` para conservarlas). Metage2Metabo no lo hace,
porque fue diseñado para redes reconstruidas con Pathway Tools, que no las
incluyen. AGORA2 sí las trae.

**Consecuencia.** Sin preprocesamiento, el seed no restringe nada y cualquier
resultado de `mincom` sobre reconstrucciones AGORA2 carece de validez. Esto
explica retroactivamente la ausencia de cross-feeding en las corridas
iniciales.

**Solución.** `preparar_sbml_para_m2m.py` elimina las reacciones de frontera
(`EX_`, `DM_`, `SK_`, `sink_`) conservando los transportadores, de modo que el
único acceso al medio sea a través del seed.

### 4. Medio derivado de los modelos

Se descartó la lista de seeds publicada por m2m para microbiota intestinal
(93 metabolitos del promedio europeo de VMH) por estar en nomenclatura MetaCyc,
incompatible con el espacio de nombres VMH de AGORA2.

En su lugar se derivó el medio de los propios modelos: metabolitos con reacción
de intercambio presentes en al menos 8 de los 15 organismos. El criterio
garantiza consistencia de identificadores por construcción.

Exclusiones aplicadas y su justificación:

| Excluido | Motivo |
|---|---|
| `biomass_c` | no es nutriente |
| `o2_e` | los organismos son anaerobios intestinales |
| `ac_e`, `for_e`, `succ_e`, `lac_L_e`, `lac_D_e`, `fum_e` | productos de fermentación; sembrarlos contamina la interpretación del cross-feeding |
| `cd2_e`, `hg2_e`, `pb2_e` | metales tóxicos; el transportador existe para expulsión |
| fármacos (`digoxin_e`, `chlphncl_e`, `czp_e`, `r406_e`, `dfdcytd_e`, …) | reacciones de metabolismo de fármacos de AGORA2, ajenas a la dieta |
| `dpmn_e`, `mtym_e` | pertenecen a la cascada; sembrarlos anula el experimento |

Medio resultante: 95 metabolitos más `levdopa_e`.

**Limitación declarada:** no es una dieta fisiológica, sino el conjunto de
metabolitos que los modelos esperan encontrar en su ambiente. Una dieta VMH
formal queda como refinamiento.

### 5. Obstáculos técnicos documentados

- La descarga individual de reconstrucciones desde
  `vmh.life/files/reconstructions/AGORA2/version2.01/sbml_files/individual_reconstructions/`
  devuelve una página HTML (~12 KB) en lugar del SBML; problema reportado desde
  octubre de 2026. Se usó la descarga desde la interfaz web.
- *Clostridium sporogenes* ATCC 15579 no se deja leer con COBRApy
  (`CobraSBMLError`). Excluido del análisis.
- `-psf` (possible seeds) en Seed2LP solo opera con `-m minimize`.
- `m2m seeds` escribe siempre en `seeds.sbml`; seeds y targets deben generarse
  en directorios distintos.

---

## Estructura del repositorio

```
proyecto_redes/
├── agora2/
│   └── sbml_files.txt              # índice de nombres de AGORA2 v2.01
├── sbml/                           # 14 reconstrucciones en SBML
├── sbml_sin_ex/                    # reconstrucciones preprocesadas
├── m2m/                            # seeds y targets
├── m2m_targets/
├── inputs_Elenta-A2/               # entradas de Seed2LP
├── results/                        # salidas
├── preparar_sbml_para_m2m.py
├── Build_seed2lp.ipynb
└── mis15.txt
```

Los archivos SBML no se versionan (ver `.gitignore`); se descargan desde
https://www.vmh.life/reconstructionDownload

## Requisitos

```
Python >= 3.11
cobra
seed2lp
metage2metabo
Gurobi (licencia académica)
gcc >= 7
```

## Reproducción

```bash
# 1. Preprocesar los modelos
python preparar_sbml_para_m2m.py sbml/ sbml_sin_ex/

# 2. Generar seeds y targets en SBML
m2m seeds --metabolites m2m/seeds.txt   -o m2m/
m2m seeds --metabolites m2m/targets.txt -o m2m_targets/

# 3. Control negativo: medio mínimo, el target no debe ser producible
m2m mincom -n sbml_sin_ex/ -s m2m_vacio/seeds.sbml -t m2m_targets/seeds.sbml \
  -o results/m2m_mincom_vacio_limpio

# 4. Selección de comunidad mínima
m2m mincom -n sbml_sin_ex/ -s m2m/seeds.sbml -t m2m_targets/seeds.sbml \
  -o results/m2m_mincom_limpio
```

## Próximos pasos

1. Verificar el control negativo sobre los modelos preprocesados.
2. Implementar la restricción de ancla sobre `mincom`.
3. Validación con SteadyCom y cuantificación de la atrición entre la etapa
   cualitativa y la cuantitativa.
4. Evaluar el modo `community` de Seed2LP como alternativa a `mincom`.

---

## Referencias

- Heinken A. et al. (2023). Genome-scale metabolic reconstruction of 7,302
  human microorganisms for personalized medicine. *Nature Biotechnology* 41(9),
  1320–1331. https://doi.org/10.1038/s41587-022-01628-0
- Maini Rekdal V. et al. (2019). Discovery and inhibition of an interspecies
  gut bacterial pathway for levodopa metabolism. *Science* 364, eaau6323.
  https://doi.org/10.1126/science.aau6323
- Belcour A. et al. (2020). Metage2Metabo, microbiota-scale metabolic
  complementarity for the identification of key species. *eLife* 9, e61968.
  https://doi.org/10.7554/eLife.61968
- Ghassemi Nedjad C. et al. (2025). Seed2LP: seed inference in metabolic
  networks for reverse ecology applications. *Bioinformatics* 41(4), btaf140.
  https://doi.org/10.1093/bioinformatics/btaf140
- Chan S.H.J. et al. (2017). SteadyCom: Predicting microbial abundances while
  ensuring community stability. *PLoS Computational Biology* 13(5), e1005539.
  https://doi.org/10.1371/journal.pcbi.1005539
- Burgard A.P. et al. (2004). Flux coupling analysis of genome-scale metabolic
  network reconstructions. *Genome Research* 14(2), 301–312.
- Saa P. et al. (2022). Modeling approaches for probing cross-feeding
  interactions in the human gut microbiome. *CSBJ* 20, 79–89.
  https://doi.org/10.1016/j.csbj.2021.12.006
