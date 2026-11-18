# Atribución especie-resuelta de biotransformaciones multi-paso en comunidades microbianas

Proyecto semestral — Optimización de Redes Metabólicas (IIQ3823)
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
redundantes, garantizando a la vez la viabilidad de crecimiento de la comunidad.

El método es independiente del caso: funciona para cualquier organismo ancla,
cualquier función metabólica objetivo y cualquier comunidad de partida. La
biotransformación de fármacos es el caso de prueba, no el método.

## Aporte metodológico

`mincom` (Metage2Metabo) busca cualquier comunidad mínima que produzca los
targets mediante expansión de red, un criterio booleano que no considera
estequiometría. Este trabajo extiende **SteadyCom** (Chan et al., 2017), que
evalúa una comunidad dada, a un **MILP que selecciona la comunidad**:

- **Variables binarias de selección acopladas a la biomasa** (`ε·y_k ≤ X_k ≤ y_k`).
  Como las cotas de flujo escalan con `X_k`, la biomasa actúa como interruptor sin
  constantes big-M adicionales (las propias cotas cumplen ese rol; ver sección 9).
- **Tasa de crecimiento común fija** (`v_bio,k = μ·X_k`), lo que mantiene el problema
  lineal, con un barrido de μ para verificar que el resultado no depende de su valor.
- **Exportación mínima del producto terminal**, sin intermediarios en el medio.
- **Restricción de ancla**: el organismo ancla está presente en toda solución, lo que
  cambia la pregunta de *"¿cuál es la comunidad mínima?"* a *"¿cuál es la comunidad
  mínima alrededor de este organismo?"*. Variante **funcional**: el ancla debe además
  portar flujo en una reacción específica de la ruta.
- **Enumeración de todas las comunidades minimales** con cortes de no superconjunto,
  que excluyen cada solución encontrada y todas las que la contienen.
- **Ranking e interpretación de roles con SteadyComFVA**: crecimiento, producción,
  rangos de abundancia, qué paso realiza cada especie y qué metabolitos se traspasan.

---

## Pipeline

| Etapa | Herramienta | Entrada | Salida |
|---|---|---|---|
| Modelos | AGORA2 v2.01 | — | 15 GSMNs (SBML); 14 legibles |
| **Etapa I (Entrega 1): análisis cualitativo** | | | |
| Preprocesamiento | `preparar_sbml_para_m2m.py` | GSMNs | GSMNs sin reacciones de frontera |
| Validación | Seed2LP | GSMN + target | conjunto mínimo de seeds |
| Medio | derivado de los modelos | GSMNs | lista de seeds |
| Selección cualitativa | Metage2Metabo (`mincom`) | seeds + targets + pool | simbiontes esenciales / alternativos |
| **Etapa II (Entrega 2): análisis estequiométrico** | | | |
| Modelo de comunidad | `src/comunidad.py` | GSMNs originales + medio | modelo SteadyCom con μ fijo |
| Selección y enumeración | `src/milp.py` | modelo de comunidad + ancla | comunidades minimales verificadas |
| Barrido, ranking y roles | `analisis_entrega2.py` | comunidades minimales | tablas en `resultados_e2/` |

La Etapa II usa las reconstrucciones **originales** (`sbml/`), no las
preprocesadas: el modelo estequiométrico necesita las reacciones de intercambio
para construir el pool común.

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

## Resultados — Etapa I (Entrega 1)

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

---

## Resultados — Etapa II (Entrega 2)

Configuración común: medio `data/medio_ampliado.csv`, captación máxima de 1
mmol h⁻¹ por nutriente (biomasa total = 1 gDW), producción mínima de m-tiramina
θ = 0,1 mmol h⁻¹, Gurobi (licencia académica).

### 5. Preprocesamiento para el modelo estequiométrico

Con las reconstrucciones originales, ningún organismo crece. Se identificaron y
corrigieron tres causas, todas en el código (los SBML no se modifican):

| Problema | Corrección | Justificación |
|---|---|---|
| `pe180_c`, `peai17_c`, `pei17_c`, `pg180_c` solo aparecen en la reacción de biomasa (sin ruta de síntesis) | se retiran de la biomasa (`reparar_biomasa`) | callejones sin salida de la reconstrucción; sin ellos el modelo puede crecer |
| reacciones de demanda (`DM_`) reversibles actúan como fuente | se fija `lb = 0` (`cerrar_demandas_reversibles`) | una demanda no debe suministrar metabolito |
| `medio.csv` no permite crecimiento (μ* = 0 para *E. faecalis* y *E. lenta*) | se usa `medio_ampliado.csv` | agrega 12 metabolitos: `12dgr180_e`, `26dap_M_e`, `2obut_e`, `cit_e`, `cys_L_e`, `glyc3p_e`, `hsulfid_e`, `lanost_e`, `nmn_e`, `ocdca_e`, `pydx_e`, `q8_e` |

**Reacciones sink.** Ambos organismos dependen de `sink_PGPm1_c_`, que suministra
un precursor de lípidos de membrana sin ruta de síntesis en las reconstrucciones:

| Escenario | *E. faecalis* μ* | *E. lenta* DSM 2243 μ* | Par (produciendo) μ* |
|---|---|---|---|
| sinks abiertas | 1,032 | 0,024 | 2,063 |
| sin `sink_PGPm1` de *E. lenta* | 1,032 | 0 | 0 |
| sin `sink_thissh` de *E. lenta* | 1,032 | 0,024 | 2,063 |
| sin `sink_PGPm1` de *E. faecalis* | 0 | 0,024 | 0 |

μ* no cambia al variar la cota de las sinks entre −0,5 y −10: la sink cubre un
requerimiento cualitativo pero no limita el crecimiento. **Supuesto:** los
organismos pueden obtener los precursores lipídicos que sus reconstrucciones no
sintetizan; se mantienen las sinks de AGORA2.

### 6. Validación de la cascada con SteadyCom (par *E. faecalis* + *E. lenta* DSM 2243)

| Prueba | Resultado |
|---|---|
| *E. faecalis* sola | crece (μ* = 1,032), no produce m-tiramina |
| *E. lenta* DSM 2243 sola | crece (μ* = 0,024), no produce m-tiramina |
| Par, exigiendo m-tiramina | μ* = 2,063 |
| Producción máxima a 0,5·μ* | 1,0 (= captación máxima de levodopa: conversión 1:1) |
| Contrafactual sin levodopa | producción = 0 |

Flujos de intercambio: *E. faecalis* consume levodopa y secreta dopamina;
*E. lenta* consume dopamina y secreta m-tiramina.

### 7. Comunidades minimales (MILP sobre 14 organismos)

| Modo | Comunidades (todas de tamaño 2) |
|---|---|
| Con ancla | *E. faecalis* + {*E. lenta* A2, DSM 11767, DSM 2243} |
| Ancla funcional (`3HLYTCL` en *E. faecalis*) | las mismas 3 |
| Sin ancla | {*E. faecalis*, *E. faecium*} × {*E. lenta* A2, DSM 11767, DSM 2243} = 6 |

- **Ninguna cepa completa la ruta sola** (no existe comunidad de tamaño 1).
- **Los 9 controles negativos no aparecen en ninguna comunidad**, incluido
  *E. coli* 042, que M2M clasificaba como simbionte esencial.
- La enumeración termina sin más soluciones: toda combinación de los 14 organismos
  que produce m-tiramina contiene un *Enterococcus* y una *E. lenta*.
- **Esenciales / alternativos (con ancla):** ninguna cepa individual es esencial;
  las 3 cepas de *E. lenta* son alternativas. A nivel de especie, *E. lenta* es
  indispensable. Sin ancla, *E. faecium* es alternativa a *E. faecalis*.
- **Costo del ancla = 0**: con o sin ancla la comunidad mínima tiene 2 miembros.

**Contraste con la Etapa I.** Por expansión de red (con semillas citosólicas),
*E. lenta* A2 y DSM 11767 parecían completar la ruta solas, y *E. coli* 042
aparecía como esencial. El modelo estequiométrico corrige ambos resultados y
coincide con la evidencia experimental (Maini Rekdal et al., 2019).

### 8. Barrido de μ

Las mismas comunidades aparecen para μ ∈ {0,01; 0,05; 0,1; 0,2; 0,5; 1,0} h⁻¹:
3 con ancla y 6 sin ancla, todas de tamaño 2. **El resultado no depende del μ
fijado.** (`resultados_e2/barrido_mu.csv`)

### 9. Hallazgo numérico: big-M implícito

Con θ = 0,001 y las tolerancias por defecto del solver, la enumeración sobre 14
organismos devolvió comunidades de un solo miembro (p. ej. *Roseburia* o
*Akkermansia* produciendo m-tiramina), biológicamente imposibles.

**Causa.** Una binaria "casi cero" dentro de la tolerancia de integralidad
(y ≈ 10⁻⁵) permite X ≈ 10⁻⁵ y, con cotas de flujo de 1000·X, flujos de hasta 0,01
> θ: un organismo trabaja sin ser contado. Las cotas de flujo actúan como big-M
implícito.

**Solución.** (1) Tolerancias de integralidad y factibilidad de 10⁻⁹;
(2) θ = 0,1, muy superior a cualquier fuga posible; (3) cada comunidad se
**verifica** re-resolviendo con las binarias fijas, y la enumeración se detiene si
una no es factible.

### 10. Ranking (SteadyComFVA)

Todas las comunidades evaluadas bajo las mismas condiciones; FVA a μ = f·μ*.
"Producción garantizada" = mínimo de exportación de m-tiramina sin exigirla.

| Comunidad | μ* | Producción [min, max] (f = 0,9) | X_min (f = 0,9) |
|---|---|---|---|
| *E. faecalis* + *E. lenta* DSM 2243 | 2,063 | [0, 1] | Efa 0,903 · Ele 0,0009 |
| *E. faecium* + *E. lenta* DSM 2243 | 1,689 | [0, 1] | Efm 0,881 · Ele 0,0009 |
| *E. faecalis* + *E. lenta* DSM 11767 | 1,645 | [0, 1] | Efa 0,582 · Ele 0,303 |
| *E. faecalis* + *E. lenta* A2 | 1,638 | [0, 1] | Efa 0,585 · Ele 0,300 |
| *E. faecium* + *E. lenta* DSM 11767 | 1,334 | [0, 1] | Efm 0,576 · Ele 0,310 |
| *E. faecium* + *E. lenta* A2 | 1,328 | [0, 1] | Efm 0,579 · Ele 0,307 |

- **Producción**: idéntica en todas (máximo limitado por la levodopa disponible;
  producción garantizada 0, es decir, no acoplada al crecimiento). No discrimina.
- **Crecimiento**: las parejas con *E. faecalis* superan a las con *E. faecium*;
  A2 y DSM 11767 son prácticamente equivalentes.
- **Abundancias**: a f = 0,5 *E. lenta* puede tener abundancia 0 (solo se requiere
  para producir); a f = 0,9 con A2 o DSM 11767 debe ser ≥ ~30 % (el *Enterococcus*
  la requiere para crecer).
- **Cautela con DSM 2243**: lidera por μ* con solo 0,09 % de abundancia, y en su
  solución aparece `o2_e` como metabolito traspasado (*E. lenta* es anaerobia).
  Se interpreta como artefacto de la reconstrucción (una biomasa mínima puede
  mover flujos de hasta 1000·X), no como la mejor comunidad.

### 11. Roles (FVA cualitativo)

En **todas** las comunidades y a f = 0,5 y 0,9:

| Organismo | Paso 1 (`3HLYTCL`) | Paso 2 (`DOPADH`) | Levodopa | Dopamina | m-tiramina |
|---|---|---|---|---|---|
| *Enterococcus* (faecalis o faecium) | obligatorio | no la tiene | siempre consume | siempre secreta | — |
| *E. lenta* A2 / DSM 11767 | **no ocurre** (FVA [0, 0]) | obligatorio | — | siempre consume | siempre secreta |
| *E. lenta* DSM 2243 | no la tiene | obligatorio | — | siempre consume | siempre secreta |

**Hallazgo.** AGORA2 asigna la descarboxilasa (`3HLYTCL`) a *E. lenta* A2 y
DSM 11767, pero estos modelos solo contienen `levdopa_c`: carecen de levodopa
extracelular y de su transportador. La reacción nunca porta flujo y el paso 1
recae siempre en *Enterococcus*, en concordancia con la evidencia experimental.
Esto resuelve la "divergencia de anotación" de la Entrega 1, cuya atribución
usaba semillas citosólicas que omitían el transporte.

**Cross-feeding recurrente** (`resultados_e2/intercambios.csv`): dopamina (en
todas), 2-oxobutanoato, purinas (adenosina, hipoxantina), lactato, succinato y
aminoácidos (leucina, prolina, asparagina, glutamato, alanina, tirosina).
Los metabolitos con rango simétrico ±97 (`czp_e`, `anzp_e`, `chlphncl_e`,
`nzp_e`, …) corresponden a ciclos sin costo y se descartan.

---

## Estructura del repositorio

```
proyecto_redes/
├── agora2/
│   └── sbml_files.txt              # índice de nombres de AGORA2 v2.01
├── sbml/                           # reconstrucciones originales (Etapa II)
├── sbml_sin_ex/                    # reconstrucciones preprocesadas (Etapa I)
├── m2m/                            # seeds y targets
├── m2m_targets/
├── inputs_Elenta-A2/               # entradas de Seed2LP
├── results/                        # salidas de la Etapa I
├── data/
│   ├── medio.csv                   # medio derivado (Etapa I)
│   └── medio_ampliado.csv          # medio usado en la Etapa II
├── src/
│   ├── comunidad.py                # modelo SteadyCom con μ fijo
│   └── milp.py                     # MILP de selección, cortes y verificación
├── ejemplo_entrega2.py             # validación del par + enumeración
├── analisis_entrega2.py            # barrido de μ, ranking y roles
├── resultados_e2/                  # salidas de la Etapa II
├── preparar_sbml_para_m2m.py
├── Build_seed2lp.ipynb
└── mis15.txt
```

Si los SBML no están en `sbml/`, se descargan desde
https://www.vmh.life/reconstructionDownload

### Salidas de la Etapa II (`resultados_e2/`)

| Archivo | Una fila es… | Contenido |
|---|---|---|
| `barrido_mu.csv` | un μ en un modo (con/sin ancla) | comunidades minimales encontradas |
| `ranking.csv` | una comunidad a una fracción f de μ* | μ*, producción [min, max], abundancias mínimas, posición |
| `abundancias.csv` | un organismo en una comunidad | rango [X_min, X_max] |
| `roles_ruta.csv` | un organismo × reacción/intercambio de la ruta | rango FVA y rol (obligatorio, posible, no ocurre…) |
| `intercambios.csv` | un metabolito traspasado × organismo | rango FVA (+ secreta, − consume) y rol |

## Requisitos

```
Python >= 3.11
cobra
seed2lp
metage2metabo
Gurobi (licencia académica)   # alternativa libre: highspy + osqp (SOLVER = "hybrid")
gcc >= 7
```

## Reproducción

### Etapa I

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

### Etapa II

```bash
# Validación del par + MILP con los 14 organismos (con ancla, ancla funcional, sin ancla)
python ejemplo_entrega2.py --todos

# Ranking y roles (~2 min con Gurobi)
python analisis_entrega2.py --solo-ranking

# Barrido de μ (MILP 12 veces sobre los 14 organismos)
python analisis_entrega2.py --solo-barrido
```

## Próximos pasos

1. Investigar el caso *E. lenta* DSM 2243: el μ* elevado con abundancia mínima y la
   aparición de oxígeno como metabolito traspasado; evaluar un ranking con
   abundancia mínima por miembro (p. ej. ≥ 1 %).
2. Identificar qué nutriente limita el crecimiento de *E. lenta* (μ* = 0,024 sola).
3. Evaluar la sensibilidad al medio (`medio.csv` vs `medio_ampliado.csv` vs dieta VMH).
4. Revisar la evidencia que sustenta la actividad de *E. lenta* DSM 2243 en AGORA2
   frente a la literatura experimental.
5. Aplicar el método a un segundo caso de biotransformación.

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
- Raghu A.K., Palanikumar I., Raman K. (2024). Designing function-specific minimal
  microbiomes from large microbial communities. *npj Systems Biology and
  Applications* 10, 46. https://doi.org/10.1038/s41540-024-00373-1
- Burgard A.P. et al. (2004). Flux coupling analysis of genome-scale metabolic
  network reconstructions. *Genome Research* 14(2), 301–312.
- Saa P. et al. (2022). Modeling approaches for probing cross-feeding
  interactions in the human gut microbiome. *CSBJ* 20, 79–89.
  https://doi.org/10.1016/j.csbj.2021.12.006
