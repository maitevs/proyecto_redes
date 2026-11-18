from pathlib import Path
from ejemplo_entrega2 import Comunidad, leer_medio, SOLVER

BASE = Path(__file__).resolve().parent

def ruta(nombre):
    return BASE / "sbml" / f"{nombre}_AGORA2_v2.01.xml"
com = Comunidad.desde_sbml(
    {
        "Efa": ruta("Enterococcus_faecalis_ATCC_29212"),
        "Ele": ruta("Eggerthella_lenta_DSM_2243"),
    },
    medio=leer_medio(BASE / "data" / "medio_ampliado.csv"),
    mu=0.1,
    solver="glpk",
)

# Exigir presencia de ambas especies
for X in com.X.values():
    X.lb = 0.001

# Exigir producción mínima de m-tiramina
com.exigir_produccion("mtym_e", 0.1)

print("\n=== VALIDACIÓN FINAL ===")

if com.factible():
    print("Comunidad viable: SÍ")
    print("Abundancias:", com.abundancias())
    print("Producción:", round(com.flujo("EXcom_mtym_e"), 4))
    print("3HLYTCL:", round(com.flujo("3HLYTCL__Efa"), 4))
    print("DOPADH:", round(com.flujo("DOPADH__Ele"), 4))
else:
    print("Comunidad viable: NO")

# Control negativo sin levodopa
com.liberar_produccion("mtym_e")
com.fijar_captacion("levdopa_e", 0)

produccion_sin = com.optimizar(
    com.ex_com["mtym_e"].flux_expression, "max"
)

print("\nProducción sin levodopa:", produccion_sin)