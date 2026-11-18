"""
Ejemplo de uso de src/comunidad.py y src/milp.py (Entrega 2).

    python ejemplo_entrega2.py            # par + 5 organismos de la ruta (rápido)
    python ejemplo_entrega2.py --todos    # los 14 modelos legibles (usar Gurobi)

Requisitos: cobra (pip install cobra) y un solver:
  - Gurobi (recomendado, licencia académica): SOLVER = "gurobi"
  - HiGHS (libre):  pip install highspy osqp  ->  SOLVER = "hybrid"
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "src"))

from comunidad import Comunidad, leer_medio
from milp import SeleccionMinima

SOLVER = "gurobi"          # cambiar a "hybrid" si no hay Gurobi
MEDIO = "data/medio_ampliado.csv"
TARGET = "mtym_e"          # m-tiramina
THETA = 0.1               # exportación mínima de m-tiramina (mmol/h); >> 1000 x tolerancia del solver
MU = 0.1                   # tasa de crecimiento común fija (1/h)
ANCLA = "Enterococcus_faecalis_ATCC_29212"


def ruta(nombre):
    return f"sbml/{nombre}_AGORA2_v2.01.xml"


def prueba_par():
    print("=== 1. SteadyCom del par E. faecalis + E. lenta DSM 2243 ===")
    com = Comunidad.desde_sbml(
        {"Efa": ruta("Enterococcus_faecalis_ATCC_29212"),
         "Ele": ruta("Eggerthella_lenta_DSM_2243")},
        medio=leer_medio(MEDIO), mu=MU, solver=SOLVER)
    print("Reparaciones aplicadas:", com.reparaciones)

    for k in ["Efa", "Ele"]:
        com.presentes([k])
        print(f"  {k} sola: mu* = {com.mu_max():.4f}, "
              f"mu* produciendo m-tiramina = {com.mu_max(producir={TARGET: THETA}):.4f}")
    com.presentes(["Efa", "Ele"])
    mu_star = com.mu_max(producir={TARGET: THETA})
    print(f"  Par: mu* produciendo m-tiramina = {mu_star:.4f}")

    com.fijar_mu(0.5 * mu_star)
    e_max = com.optimizar(com.ex_com[TARGET].flux_expression, "max")
    print(f"  max e_mtym a 0.5 mu* = {e_max:.4f}; X = {com.abundancias()}")
    for k, met, v in com.resumen_intercambios():
        if met in ("levdopa_e", "dpmn_e", "mtym_e"):
            print(f"    {k:4s} {met:10s} {v:+.4f}   (+ secreta, - consume)")

    com.fijar_captacion("levdopa_e", 0.0)
    print("  Contrafactual sin levodopa: max e_mtym =",
          com.optimizar(com.ex_com[TARGET].flux_expression, "max"))


def enumerar(nombres):
    print(f"\n=== 2. MILP sobre {len(nombres)} organismos ===")
    t = time.time()
    com = Comunidad.desde_sbml({n: ruta(n) for n in nombres},
                               medio=leer_medio(MEDIO), mu=MU, solver=SOLVER)
    print(f"  construido en {time.time() - t:.0f} s")
    com.exigir_produccion(TARGET, THETA)
    sel = SeleccionMinima(com, ancla=ANCLA)

    print("-- ancla presente")
    sel.enumerar()
    sel.limpiar_cortes()

    print("-- ancla funcional (E. faecalis debe hacer 3HLYTCL)")
    sel.ancla_funcional("3HLYTCL", 1e-3)
    sel.enumerar()
    sel.limpiar_cortes()
    sel.quitar_ancla_funcional("3HLYTCL")

    print("-- sin ancla (para el costo del ancla)")
    sel.y[ANCLA].lb = 0
    sel.ancla = None
    sel.enumerar()


if __name__ == "__main__":
    prueba_par()
    if "--todos" in sys.argv:
        nombres = [p.name.replace("_AGORA2_v2.01.xml", "") for p in sorted(Path("sbml").glob("*.xml"))
                   if "Clostridium" not in p.name]   # SBML inválido
    else:
        nombres = [ANCLA, "Enterococcus_faecium_504", "Eggerthella_lenta_DSM_2243",
                   "Eggerthella_lenta_A2", "Eggerthella_lenta_DSM_11767"]
    enumerar(nombres)