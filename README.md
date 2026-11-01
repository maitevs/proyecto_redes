# Atribución especie-resuelta de biotransformaciones multi-paso en comunidades microbianas

Proyecto semestral — optimización de redes metabólicas
Pontificia Universidad Católica de Chile, DIQB

## Aporte metodológico

`mincom` (Metage2Metabo) busca cualquier comunidad mínima que produzca los
targets. Aquí se impone la restricción de que un organismo determinado esté
presente en toda solución. Esto cambia la pregunta de *"¿cuál es la comunidad
mínima?"* a *"¿cuál es la comunidad mínima alrededor de este organismo?"*.

De esa restricción se derivan las demás decisiones: el medio se infiere desde los
requerimientos del organismo ancla (Seed2LP) en vez de fijarse genéricamente, los
targets son metabolitos de fármaco en vez de biomasa, y se agrega un filtro
estequiométrico porque la selección cualitativa no garantiza crecimiento.

---

## Pipeline

| Etapa | Herramienta | Entrada | Salida |
|---|---|---|---|
| Modelos | AGORA2 v2.01 | — | 15 GSMNs (SBML) |
| Inferencia de medio | Seed2LP | GSMN + target | conjunto de seeds |
| Selección cualitativa | Metage2Metabo (`mincom`) | seeds + targets + pool | simbiontes esenciales / alternativos |
| Validación cuantitativa | SteadyCom | comunidad seleccionada | viabilidad de crecimiento |

Las dos primeras etapas operan sobre una abstracción booleana (network
expansion): un metabolito es producible si es alcanzable desde los seeds. No hay
estequiometría ni balance de masa. La última etapa sí es estequiométrica.

---

## Referencias

- Heinken A. et al. (2023). Genome-scale metabolic reconstruction of 7,302 human
  microorganisms for personalized medicine. *Nature Biotechnology* 41(9),
  1320–1331. https://doi.org/10.1038/s41587-022-01628-0
- Maini Rekdal V. et al. (2019). Discovery and inhibition of an interspecies gut
  bacterial pathway for levodopa metabolism. *Science* 364, eaau6323.
  https://doi.org/10.1126/science.aau6323
- Belcour A. et al. (2020). Metage2Metabo, microbiota-scale metabolic
  complementarity for the identification of key species. *eLife* 9, e61968.
  https://doi.org/10.7554/eLife.61968
- Ghassemi Nedjad C. et al. (2025). Seed2LP: seed inference in metabolic networks
  for reverse ecology applications. *Bioinformatics* 41(4), btaf140.
  https://doi.org/10.1093/bioinformatics/btaf140
- Chan S.H.J. et al. (2017). SteadyCom: Predicting microbial abundances while
  ensuring community stability. *PLoS Computational Biology* 13(5), e1005539.
  https://doi.org/10.1371/journal.pcbi.1005539
- Burgard A.P. et al. (2004). Flux coupling analysis of genome-scale metabolic
  network reconstructions. *Genome Research* 14(2), 301–312.
- Saa P. et al. (2022). Modeling approaches for probing cross-feeding
  interactions in the human gut microbiome. *CSBJ* 20, 79–89.
  https://doi.org/10.1016/j.csbj.2021.12.006