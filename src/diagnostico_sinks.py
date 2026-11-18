from pathlib import Path
import cobra
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
SBML = BASE / "sbml"
DATA = BASE / "data"

# Modelo original de E. lenta
ruta = SBML / "Eggerthella_lenta_DSM_2243_AGORA2_v2.01.xml"

modelo = cobra.io.read_sbml_model(str(ruta))

print("=== REACCIONES SINK Y DEMAND ===")

resultados = []

for rxn in modelo.reactions:

    if rxn.id.lower().startswith(("sink_", "sk_", "dm_")):

        resultados.append({
            "reaccion": rxn.id,
            "ecuacion": rxn.reaction,
            "lb": rxn.lower_bound,
            "ub": rxn.upper_bound,
            "reversible": rxn.reversibility
        })

df = pd.DataFrame(resultados)

print("Reacciones encontradas:", len(df))

if not df.empty:
    print(df.to_string(index=False))
    df.to_csv(
        DATA / "diagnostico_sinks_Elenta.csv",
        index=False
    )

print("\n=== IMPACTO DE SINKS EN CRECIMIENTO ===")

lipidos_bloqueados = [
    "pe180_c",
    "peai17_c",
    "pei17_c",
    "pg180_c"
]

sinks_sospechosas = [
    "sink_PGPm1_c_",
    "sink_thissh_c_"
]


def preparar_modelo_diagnostico(original):
    m = original.copy()

    # Medio abierto, solo para diagnóstico
    for rxn in m.exchanges:
        if len(rxn.metabolites) == 1:
            coef = next(iter(rxn.metabolites.values()))

            if coef == -1:
                rxn.lower_bound = -1000.0

    # Biomasa modificada, solo para diagnóstico
    bio = m.reactions.get_by_id("biomass467")

    for met_id in lipidos_bloqueados:
        met = m.metabolites.get_by_id(met_id)

        if met in bio.metabolites:
            coef = bio.metabolites[met]
            bio.add_metabolites({met: -coef})

    m.objective = bio
    m.objective_direction = "max"

    return m


escenarios = {
    "Ambas sinks abiertas": [],
    "Sin sink_PGPm1": ["sink_PGPm1_c_"],
    "Sin sink_thissh": ["sink_thissh_c_"],
    "Ambas sinks cerradas": sinks_sospechosas
}

for nombre, cerrar in escenarios.items():

    m = preparar_modelo_diagnostico(modelo)

    for rxn_id in cerrar:
        rxn = m.reactions.get_by_id(rxn_id)

        # Impedir generación artificial del metabolito.
        # Mantener posible su eliminación.
        rxn.lower_bound = 0.0

    sol = m.optimize()

    print(f"\n{nombre}")
    print("Estado:", sol.status)

    if sol.status == "optimal":
        print("Crecimiento:", round(sol.objective_value, 6))

        for rid in sinks_sospechosas:
            print(
                rid, ":",
                round(sol.fluxes[rid], 6)
            )