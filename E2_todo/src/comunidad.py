"""
Modelo de comunidad tipo SteadyCom (Chan et al., 2017) con tasa de crecimiento fija.

Construye, a partir de varios GEMs de AGORA2 en SBML, un único modelo COBRA donde:

  - cada organismo k conserva sus reacciones (sufijo __k) y su espacio extracelular propio;
  - los intercambios EX_ de cada organismo se convierten en intercambios con un pool
    común "u" (variable v_ext[k, i], positiva si k secreta i al pool);
  - el pool intercambia con el ambiente mediante EXcom_i (negativo = entra desde el medio);
  - cada organismo tiene una biomasa X_k (abundancia relativa, sum X_k = 1);
  - las cotas de flujo escalan con la biomasa:  lb * X_k <= v <= ub * X_k;
  - la biomasa de cada organismo crece a la tasa común fija: v_bio_k = mu * X_k.

Con mu fijo todas las restricciones son lineales. Este mismo objeto es la base del MILP
(milp.py), que agrega las variables binarias de selección sobre X_k.

Uso mínimo:

    com = Comunidad.desde_sbml({"Efa": "sbml/Enterococcus...xml",
                                "Ele": "sbml/Eggerthella...xml"},
                               medio=leer_medio("data/medio.csv"))
    mu_star = com.mu_max(producir={"mtym_e": 1e-3})
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path

import cobra
from cobra import Metabolite, Model, Reaction

logging.getLogger("cobra").setLevel(logging.ERROR)

INF = 1000.0  # cota "infinita" de AGORA2


# ---------------------------------------------------------------------------
# Lectura y reparación de modelos individuales
# ---------------------------------------------------------------------------

def leer_modelo(path: str | Path) -> Model:
    return cobra.io.read_sbml_model(str(path))


def reaccion_biomasa(model: Model) -> Reaction:
    """Reacción con coeficiente distinto de cero en el objetivo del SBML."""
    rx = [r for r in model.reactions if r.objective_coefficient != 0]
    if len(rx) != 1:
        raise ValueError(f"{model.id}: se esperaba 1 reacción de biomasa, hay {len(rx)}")
    return rx[0]


def reparar_biomasa(model: Model) -> list[str]:
    """
    Quita de la biomasa los metabolitos que SOLO aparecen en ella (callejones sin salida).

    En la exportación SBML de VMH usada (release 2026-09-18), pe180_c, peai17_c, pei17_c y
    pg180_c no tienen ninguna reacción de síntesis, por lo que el modelo no puede crecer.
    Devuelve la lista de metabolitos eliminados para dejarlo documentado.
    """
    bio = reaccion_biomasa(model)
    muertos = {m: c for m, c in bio.metabolites.items() if c < 0 and len(m.reactions) == 1}
    if muertos:
        bio.subtract_metabolites(muertos)
    return sorted(m.id for m in muertos)


def cerrar_demandas_reversibles(model: Model) -> list[str]:
    """Una reacción DM_ (demanda) no debería suministrar metabolito: fija lb = 0."""
    cerradas = []
    for r in model.reactions:
        if r.id.startswith("DM_") and r.lower_bound < 0:
            r.lower_bound = 0.0
            cerradas.append(r.id)
    return cerradas


# ---------------------------------------------------------------------------
# Medio
# ---------------------------------------------------------------------------

def leer_medio(path: str | Path) -> dict[str, float]:
    """CSV con columnas: metabolito (id sin prefijo M_, p.ej. glc_D_e), captacion_max."""
    medio = {}
    with open(path, newline="", encoding="utf-8") as f:
        for fila in csv.DictReader(f):
            medio[fila["metabolito"].strip()] = float(fila["captacion_max"])
    return medio


# ---------------------------------------------------------------------------
# Comunidad
# ---------------------------------------------------------------------------

@dataclass
class Comunidad:
    model: Model                                   # modelo COBRA de la comunidad
    organismos: list[str]                          # etiquetas cortas, p.ej. ["Efa", "Ele"]
    biomasa: dict = field(default_factory=dict)    # etiqueta -> Reaction de biomasa
    v_ext: dict = field(default_factory=dict)      # (etiqueta, met_e) -> Reaction organismo<->pool
    ex_com: dict = field(default_factory=dict)     # met_e -> Reaction pool<->ambiente
    reparaciones: dict = field(default_factory=dict)
    mu: float = 0.0

    # Las variables y restricciones se buscan por nombre en cada uso: si se cambia el
    # solver (com.model.solver = "gurobi"), optlang recrea los objetos y las referencias
    # guardadas quedarían obsoletas.
    @property
    def X(self) -> dict:
        return {k: self.model.variables[f"X__{k}"] for k in self.organismos}

    @property
    def restr_mu(self) -> dict:
        return {k: self.model.constraints[f"mu__{k}"] for k in self.organismos}

    # ------------------------------------------------------------------ construcción
    @classmethod
    def desde_sbml(
        cls,
        archivos: dict[str, str | Path],
        medio: dict[str, float],
        mu: float = 0.1,
        solver: str | None = None,
        exportables: tuple[str, ...] = ("levdopa_e",),
    ) -> "Comunidad":
        """
        archivos:    {etiqueta: ruta SBML}
        medio:       {met_e: captación máxima ū_i (mmol h^-1, biomasa total = 1 gDW)}.
                     Los metabolitos fuera del medio tienen ū_i = 0 (solo pueden salir).
        exportables: metabolitos extracelulares sin EX_ en AGORA2 a los que igual se
                     conecta con el pool si el organismo los tiene (levodopa no trae EX_).
        """
        com_model = Model("comunidad")
        if solver:
            com_model.solver = solver
        com = cls(model=com_model, organismos=list(archivos))
        pool: dict[str, Metabolite] = {}

        def met_pool(met_e: str, plantilla: Metabolite) -> Metabolite:
            if met_e not in pool:
                pool[met_e] = Metabolite(f"{met_e}__u", name=plantilla.name,
                                         formula=plantilla.formula, compartment="u")
            return pool[met_e]

        for k, ruta in archivos.items():
            m = leer_modelo(ruta)
            com.reparaciones[k] = {
                "biomasa_sin": reparar_biomasa(m),
                "DM_cerradas": cerrar_demandas_reversibles(m),
            }
            bio_id = reaccion_biomasa(m).id

            # EX_biomass (salida de biomasa) no es un metabolito intercambiable: queda como
            # salida propia del organismo y no se conecta al pool.
            ex_ids = {r.id for r in m.reactions
                      if r.id.startswith("EX_") and "biomass" not in r.id}
            # metabolito extracelular de cada EX_
            ex_met = {rid: next(iter(m.reactions.get_by_id(rid).metabolites)).id for rid in ex_ids}
            # extracelulares "exportables" que no tienen EX_ (p.ej. levodopa)
            con_ex = set(ex_met.values())
            extra = [mid for mid in exportables if mid in m.metabolites and mid not in con_ex]

            # renombrar todo con sufijo de organismo
            for met in m.metabolites:
                met.id = f"{met.id}__{k}"
                met.compartment = f"{met.compartment}_{k}"
            for r in m.reactions:
                r.id = f"{r.id}__{k}"
            m.repair()

            nuevas = []
            for r in m.reactions:
                base = r.id[: -len(f"__{k}")]
                if base in ex_ids:
                    # EX_ del organismo -> intercambio organismo <-> pool
                    met_e = ex_met[base]
                    met_org = m.metabolites.get_by_id(f"{met_e}__{k}")
                    rx = Reaction(f"EXsp_{met_e}__{k}", lower_bound=-INF, upper_bound=INF)
                    rx.add_metabolites({met_org: -1, met_pool(met_e, met_org): 1})
                    nuevas.append(rx)
                    com.v_ext[(k, met_e)] = rx
                else:
                    nuevas.append(r.copy())
            for met_e in extra:
                met_org = m.metabolites.get_by_id(f"{met_e}__{k}")
                rx = Reaction(f"EXsp_{met_e}__{k}", lower_bound=-INF, upper_bound=INF)
                rx.add_metabolites({met_org: -1, met_pool(met_e, met_org): 1})
                nuevas.append(rx)
                com.v_ext[(k, met_e)] = rx

            com_model.add_reactions(nuevas)
            com.biomasa[k] = com_model.reactions.get_by_id(f"{bio_id}__{k}")

        # intercambio del pool con el ambiente: EXcom_i, lb = -ū_i
        for met_e, met in pool.items():
            rx = Reaction(f"EXcom_{met_e}", lower_bound=-medio.get(met_e, 0.0), upper_bound=INF)
            rx.add_metabolites({met: -1})
            com.ex_com[met_e] = rx
        com_model.add_reactions(list(com.ex_com.values()))
        faltan = sorted(set(medio) - set(pool))
        if faltan:
            com.reparaciones["medio_sin_pool"] = faltan  # nadie los puede usar

        com._agregar_biomasas_y_escalamiento(mu)
        com_model.objective = com_model.problem.Objective(0, direction="max")
        return com

    def _agregar_biomasas_y_escalamiento(self, mu: float) -> None:
        mdl = self.model
        prob = mdl.problem
        nuevas = []
        mdl.add_cons_vars([prob.Variable(f"X__{k}", lb=0, ub=1) for k in self.organismos])
        Xs = self.X
        for k in self.organismos:
            X = Xs[k]
            sufijo = f"__{k}"
            for r in mdl.reactions:
                if not r.id.endswith(sufijo):
                    continue
                lb, ub = r.lower_bound, r.upper_bound
                fe = r.flux_expression
                if ub != 0:
                    nuevas.append(prob.Constraint(fe - ub * X, ub=0, name=f"ub_{r.id}"))
                if lb != 0:
                    nuevas.append(prob.Constraint(fe - lb * X, lb=0, name=f"lb_{r.id}"))
            c = prob.Constraint(self.biomasa[k].flux_expression - mu * X, lb=0, ub=0,
                                name=f"mu__{k}")
            nuevas.append(c)
        nuevas.append(prob.Constraint(sum(Xs.values()), lb=1, ub=1, name="suma_X"))
        mdl.add_cons_vars(nuevas)
        self.mu = mu

    # ------------------------------------------------------------------ parámetros
    def fijar_mu(self, mu: float) -> None:
        """Cambia la tasa común: v_bio_k - mu X_k = 0."""
        Xs = self.X
        for k, c in self.restr_mu.items():
            c.set_linear_coefficients({Xs[k]: -mu})
        self.mu = mu

    def fijar_captacion(self, met_e: str, umax: float) -> None:
        self.ex_com[met_e].lower_bound = -umax

    def exigir_produccion(self, met_e: str, theta: float) -> None:
        """e_tau >= theta (salida neta del pool al ambiente)."""
        self.ex_com[met_e].lower_bound = theta

    def liberar_produccion(self, met_e: str, umax: float = 0.0) -> None:
        self.ex_com[met_e].lower_bound = -umax

    def presentes(self, etiquetas) -> None:
        """Fija X_k = 0 para los organismos fuera de 'etiquetas' (comunidad dada)."""
        for k, X in self.X.items():
            X.ub = 1 if k in etiquetas else 0

    # ------------------------------------------------------------------ resolución
    limite_tiempo: int = 120   # segundos por resolución (evita que un LP se quede pegado)

    def _resolver(self) -> str:
        self.model.solver.configuration.timeout = self.limite_tiempo
        self.model.slim_optimize()
        estado = self.model.solver.status
        if estado not in ("optimal", "infeasible"):
            logging.getLogger(__name__).warning("Estado del solver: %s (mu=%s)", estado, self.mu)
            print(f"AVISO: el solver terminó con estado '{estado}' (mu={self.mu})")
        return estado

    def factible(self) -> bool:
        self.model.objective = self.model.problem.Objective(0, direction="max")
        return self._resolver() == "optimal"

    def optimizar(self, expr, sentido: str = "max") -> float | None:
        self.model.objective = self.model.problem.Objective(expr, direction=sentido)
        if self._resolver() != "optimal":
            return None
        return self.model.solver.objective.value

    def mu_max(self, mu_hi: float = 5.0, tol: float = 1e-4,
               producir: dict[str, float] | None = None) -> float:
        """
        Crecimiento común máximo por bisección (como SteadyCom), opcionalmente exigiendo
        producción: producir={"mtym_e": theta}. Devuelve 0 si no es factible ni con mu ~ 0.
        """
        mu_original = self.mu
        respaldo = {}
        for met_e, theta in (producir or {}).items():
            respaldo[met_e] = self.ex_com[met_e].bounds
            self.exigir_produccion(met_e, theta)
        try:
            self.fijar_mu(tol)
            if not self.factible():
                return 0.0
            lo, hi = tol, mu_hi
            self.fijar_mu(hi)
            if self.factible():             # mu_hi es una cota de diseño, no se amplía
                return mu_hi
            while hi - lo > tol:
                mid = (lo + hi) / 2
                self.fijar_mu(mid)
                if self.factible():
                    lo = mid
                else:
                    hi = mid
            return lo
        finally:
            for met_e, b in respaldo.items():
                self.ex_com[met_e].bounds = b
            self.fijar_mu(mu_original)

    # ------------------------------------------------------------------ lectura
    def flujo(self, rid: str) -> float:
        return self.model.reactions.get_by_id(rid).flux

    def abundancias(self) -> dict[str, float]:
        return {k: X.primal for k, X in self.X.items()}

    def resumen_intercambios(self, umbral: float = 1e-6) -> list[tuple[str, str, float]]:
        """(organismo, metabolito, flujo) con |flujo| > umbral en la última solución."""
        out = []
        for (k, met_e), r in self.v_ext.items():
            v = r.flux
            if abs(v) > umbral:
                out.append((k, met_e, v))
        return sorted(out, key=lambda t: (t[1], t[0]))
