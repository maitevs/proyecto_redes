from pathlib import Path
import cobra
import pandas as pd

# Directorio principal del proyecto
BASE = Path(__file__).resolve().parent.parent

DATA = BASE / "data"
SBML = BASE / "sbml"


def cargar_modelo(nombre_archivo):
    """
    Carga un modelo metabólico AGORA2 sin modificarlo.

    Parámetros
    ----------
    nombre_archivo : str
        Nombre del archivo SBML.

    Retorna
    -------
    cobra.Model
        Modelo metabólico cargado.
    """
    ruta = SBML / nombre_archivo

    if not ruta.is_file():
        raise FileNotFoundError(f"No existe el modelo: {ruta}")

    return cobra.io.read_sbml_model(str(ruta))


def cargar_medio(nombre_archivo="medio.csv"):
    """
    Carga el medio metabólico preliminar.
    """
    ruta = DATA / nombre_archivo

    if not ruta.is_file():
        raise FileNotFoundError(f"No existe el medio: {ruta}")

    return pd.read_csv(ruta)

def aplicar_medio(modelo, medio):
    """Aplica el medio a una copia del modelo."""
    m = modelo.copy()

    intercambios = {}

    for rxn in m.exchanges:
        if len(rxn.metabolites) == 1:
            met = next(iter(rxn.metabolites))
            intercambios[met.id] = rxn

    # Cerrar todas las captaciones
    for rxn in m.exchanges:
        rxn.lower_bound = max(0, rxn.lower_bound)

    # Abrir los nutrientes autorizados
    for _, fila in medio.iterrows():
        met_id = fila["metabolito"]

        if met_id in intercambios:
            rxn = intercambios[met_id]
            coef = next(iter(rxn.metabolites.values()))

            if coef != -1:
                raise ValueError(f"Revisar signo de {rxn.id}")

            rxn.lower_bound = float(fila["lb"])
            rxn.upper_bound = float(fila["ub"])

    return m

from cobra.flux_analysis import flux_variability_analysis


def analizar_reacciones(modelo, reacciones):
    """
    Calcula los flujos mínimos y máximos de las
    reacciones de interés mediante FVA.

    No exige crecimiento positivo.
    """

    presentes = [
        rid for rid in reacciones
        if rid in modelo.reactions
    ]

    if not presentes:
        return pd.DataFrame()

    fva = flux_variability_analysis(
        modelo,
        reaction_list=presentes,
        fraction_of_optimum=0.0
    )

    fva.index.name = "reaccion"
    return fva.reset_index()

def construir_comunidad(modelos):
    """
    Combina modelos manteniendo separados sus metabolitos
    y reacciones mediante prefijos.

    modelos: diccionario {prefijo: modelo}
    """
    comunidad = cobra.Model("comunidad_levodopa")
    reacciones_nuevas = []

    for prefijo, modelo in modelos.items():
        metabolitos = {}

        for met in modelo.metabolites:
            nuevo_met = cobra.Metabolite(
                id=f"{prefijo}__{met.id}",
                name=met.name,
                compartment=f"{prefijo}_{met.compartment}"
            )
            metabolitos[met.id] = nuevo_met

        for rxn in modelo.reactions:
            nueva_rxn = cobra.Reaction(
                id=f"{prefijo}__{rxn.id}",
                name=rxn.name,
                lower_bound=rxn.lower_bound,
                upper_bound=rxn.upper_bound
            )

            nueva_rxn.add_metabolites({
                metabolitos[met.id]: coef
                for met, coef in rxn.metabolites.items()
            })

            reacciones_nuevas.append(nueva_rxn)

    comunidad.add_reactions(reacciones_nuevas)

    return comunidad


if __name__ == "__main__":

    # PRUEBA 1: E. FAECALIS
    print("=== VALIDACIÓN DE ANALYSIS.PY ===")

    modelo = cargar_modelo(
        "Enterococcus_faecalis_ATCC_29212_AGORA2_v2.01.xml"
    )

    medio = cargar_medio()
    modelo_controlado = aplicar_medio(modelo, medio)

    print("Modelo:", modelo.id)
    print("Reacciones:", len(modelo.reactions))
    print("Metabolitos del medio:", len(medio))

    solucion = modelo_controlado.optimize()

    print("Estado:", solucion.status)
    print("Crecimiento:", solucion.objective_value)

    # PRUEBA 2: E. LENTA
    print("\n=== VALIDACIÓN CON E. LENTA ===")

    lenta = cargar_modelo(
        "Eggerthella_lenta_DSM_2243_AGORA2_v2.01.xml"
    )

    lenta_controlada = aplicar_medio(lenta, medio)

    print("Modelo:", lenta.id)
    print("Reacciones:", len(lenta.reactions))

    intercambios_lenta = {
        next(iter(r.metabolites)).id: r
        for r in lenta_controlada.exchanges
        if len(r.metabolites) == 1
    }

    for met_id in ["levdopa_e", "dpmn_e", "mtym_e"]:
        rxn = intercambios_lenta.get(met_id)

        if rxn is not None:
            print(met_id, ":", rxn.bounds)
        else:
            print(met_id, ": sin intercambio")

    sol_lenta = lenta_controlada.optimize()

    print("Estado:", sol_lenta.status)
    print("Crecimiento:", sol_lenta.objective_value)

    print("\n=== FVA DE BIOTRANSFORMACIÓN ===")

    reacciones_objetivo = [
        "3HLYTCL",
        "DOPADH",
        "DOPA34DHPHEt",
        "DOPA34DHPHEt6",
        "DOPAENT4tc",
        "MTYMt2r"
    ]

    for nombre, m in [
        ("E. faecalis", modelo_controlado),
        ("E. lenta", lenta_controlada)
    ]:
        print(f"\n--- {nombre} ---")

        resultado = analizar_reacciones(
            m,
            reacciones_objetivo
        )

        print(resultado.to_string(index=False))

    print("\n=== CONTROL POSITIVO: DOPAMINA EN E. LENTA ===")

    # Trabajar sobre una copia
    lenta_con_dopamina = lenta_controlada.copy()

    # Identificar intercambio de dopamina por metabolito
    ex_dopamina = next(
        (
            r for r in lenta_con_dopamina.exchanges
            if any(
                met.id == "dpmn_e"
                for met in r.metabolites
            )
        ),
        None
    )

    if ex_dopamina is None:
        raise ValueError("No se encontró intercambio de dopamina")

    # Permitir captación de dopamina como control diagnóstico
    ex_dopamina.lower_bound = -10.0

    # Evaluar reacciones relevantes
    resultado_control = analizar_reacciones(
        lenta_con_dopamina,
        ["DOPADH", "DOPAENT4tc", "MTYMt2r"]
    )

    print(resultado_control.to_string(index=False))

    print("\n=== PRODUCCIÓN DE M-TIRAMINA ===")

    # Copia del modelo con dopamina disponible
    m_prod = lenta_con_dopamina.copy()

    # Identificar intercambio de m-tiramina
    ex_mtym = next(
        r for r in m_prod.exchanges
        if any(met.id == "mtym_e" for met in r.metabolites)
    )

    # Maximizar la exportación de m-tiramina
    m_prod.objective = ex_mtym
    m_prod.objective_direction = "max"

    sol_prod = m_prod.optimize()

    print("Estado:", sol_prod.status)

    if sol_prod.status == "optimal":
        print(
            "Exportación máxima m-tiramina:",
            sol_prod.objective_value
        )

        for rid in ["DOPADH", "DOPAENT4tc", "MTYMt2r"]:
            print(rid, ":", sol_prod.fluxes[rid])

        print(
            "Intercambio dopamina:",
            sol_prod.fluxes[ex_dopamina.id]
        )

        print(
            "Intercambio m-tiramina:",
            sol_prod.fluxes[ex_mtym.id]
        )

        print("\n=== PRODUCCIÓN DE DOPAMINA EN E. FAECALIS ===")

    # Usar una copia con el medio controlado
    m_dopa = modelo_controlado.copy()

    # Buscar intercambio de dopamina
    ex_dopa = next(
        r for r in m_dopa.exchanges
        if any(met.id == "dpmn_e" for met in r.metabolites)
    )

    # Maximizar secreción de dopamina
    m_dopa.objective = ex_dopa
    m_dopa.objective_direction = "max"

    sol_dopa = m_dopa.optimize()

    print("Estado:", sol_dopa.status)

    if sol_dopa.status == "optimal":
        print(
            "Exportación máxima dopamina:",
            round(sol_dopa.objective_value, 6)
        )

        for rid in [
            "3HLYTCL",
            "DOPA34DHPHEt",
            "DOPA34DHPHEt6"
        ]:
            print(
                rid, ":",
                round(sol_dopa.fluxes[rid], 6)
            )

        print(
            "Intercambio dopamina:",
            round(sol_dopa.fluxes[ex_dopa.id], 6)
        )

        print("\n=== GUARDANDO RESULTADOS ===")

    salida = DATA / "resultados_validacion_individual.csv"

    resultados = pd.DataFrame([
        {
            "organismo": "E. faecalis",
            "sustrato": "levodopa",
            "producto": "dopamina",
            "reaccion": "3HLYTCL",
            "flujo_reaccion": sol_dopa.fluxes["3HLYTCL"],
            "flujo_exportacion": sol_dopa.fluxes[ex_dopa.id],
            "crecimiento": solucion.objective_value
        },
        {
            "organismo": "E. lenta",
            "sustrato": "dopamina",
            "producto": "m-tiramina",
            "reaccion": "DOPADH",
            "flujo_reaccion": sol_prod.fluxes["DOPADH"],
            "flujo_exportacion": sol_prod.fluxes[ex_mtym.id],
            "crecimiento": sol_lenta.objective_value
        }
    ])

    resultados.to_csv(salida, index=False)

    print("Archivo guardado:", salida)

    print("\n=== VALIDACIÓN COMUNITARIA ===")

    comunidad = construir_comunidad({
        "EF": modelo_controlado,
        "EL": lenta_controlada
    })

    # Identificar intercambios de dopamina
    ex_dopa_ef = comunidad.reactions.get_by_id(
        "EF__EX_dopa_e_"
    )

    ex_dopa_el = comunidad.reactions.get_by_id(
        "EL__EX_dopa_e_"
    )

    # E. lenta puede recibir dopamina del pool compartido
    ex_dopa_el.lower_bound = -10.0

    # E. faecalis solo puede secretarla
    ex_dopa_ef.lower_bound = 0.0

    # Balance del pool compartido de dopamina
    balance_dopa = comunidad.problem.Constraint(
        ex_dopa_ef.flux_expression
        + ex_dopa_el.flux_expression,
        lb=0,
        ub=0,
        name="balance_dopamina_compartida"
    )

    comunidad.add_cons_vars([balance_dopa])

    # Maximizar m-tiramina exportada por E. lenta
    ex_mtym_el = comunidad.reactions.get_by_id(
        "EL__EX_mtym_e_"
    )

    comunidad.objective = ex_mtym_el
    comunidad.objective_direction = "max"

    sol_com = comunidad.optimize()

    print("Estado:", sol_com.status)

    if sol_com.status == "optimal":

        print(
            "m-tiramina exportada:",
            round(sol_com.objective_value, 6)
        )

        for rid in [
            "EF__3HLYTCL",
            "EF__EX_dopa_e_",
            "EL__EX_dopa_e_",
            "EL__DOPADH",
            "EL__EX_mtym_e_"
        ]:
            print(
                rid, ":",
                round(sol_com.fluxes[rid], 6)
            )

        # Control negativo: sin entrada de levodopa
        comunidad_sin_levodopa = comunidad.copy()

        ex_lev = comunidad_sin_levodopa.reactions.get_by_id(
            "EF__EX_34dhphe_e_"
        )

        ex_lev.lower_bound = 0.0

        sol_neg = comunidad_sin_levodopa.optimize()

        print("\n=== CONTROL SIN LEVODOPA ===")
        print("Estado:", sol_neg.status)

        if sol_neg.status == "optimal":
            print(
                "m-tiramina exportada:",
                round(sol_neg.objective_value, 6)
            )

            print("\n=== GUARDANDO VALIDACIÓN COMUNITARIA ===")

    resultados_comunidad = pd.DataFrame([
        {
            "escenario": "Con levodopa",
            "estado": sol_com.status,
            "flujo_3HLYTCL": sol_com.fluxes["EF__3HLYTCL"],
            "dopamina_exportada_EF": sol_com.fluxes["EF__EX_dopa_e_"],
            "dopamina_captada_EL": sol_com.fluxes["EL__EX_dopa_e_"],
            "flujo_DOPADH": sol_com.fluxes["EL__DOPADH"],
            "m_tiramina_exportada": sol_com.fluxes["EL__EX_mtym_e_"]
        },
        {
            "escenario": "Sin levodopa",
            "estado": sol_neg.status,
            "flujo_3HLYTCL": sol_neg.fluxes["EF__3HLYTCL"],
            "dopamina_exportada_EF": sol_neg.fluxes["EF__EX_dopa_e_"],
            "dopamina_captada_EL": sol_neg.fluxes["EL__EX_dopa_e_"],
            "flujo_DOPADH": sol_neg.fluxes["EL__DOPADH"],
            "m_tiramina_exportada": sol_neg.fluxes["EL__EX_mtym_e_"]
        }
    ])

    salida = DATA / "resultados_validacion_comunitaria.csv"
    resultados_comunidad.to_csv(salida, index=False)

    print("Archivo guardado:", salida)