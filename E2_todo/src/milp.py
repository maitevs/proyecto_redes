"""
MILP de selección de comunidad mínima anclada, sobre el modelo SteadyCom de comunidad.py.

    min  sum_{k != ancla} y_k
    s.a. restricciones de SteadyCom con mu fijo           (comunidad.py)
         e_tau >= theta                                     (exigir producción)
         eps * y_k <= X_k <= y_k                            (acople binaria-biomasa)
         y_ancla = 1
         [opcional] v_{ancla, paso} >= delta * X_ancla      (ancla funcional)
         y_k en {0,1}

Si y_k = 0, entonces X_k = 0 y, por las cotas escaladas lb*X <= v <= ub*X, todos los
flujos del organismo k quedan en cero: la biomasa actúa como interruptor sin big-M.

Enumeración: tras encontrar la comunidad C, se agrega el corte de no superconjunto
    sum_{k in C, k != ancla} (1 - y_k) >= 1
que excluye C y toda comunidad que la contenga.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from comunidad import Comunidad


@dataclass
class SeleccionMinima:
    com: Comunidad
    ancla: str | None
    eps: float = 1e-3
    cortes: list = field(default_factory=list)

    def __post_init__(self):
        mdl = self.com.model
        p = mdl.problem
        ys = [p.Variable(f"y__{k}", type="binary") for k in self.com.organismos]
        mdl.add_cons_vars(ys)
        Xs, Y = self.com.X, self.y
        cons = []
        for k in self.com.organismos:
            cons.append(p.Constraint(Xs[k] - Y[k], ub=0, name=f"X_le_y__{k}"))
            cons.append(p.Constraint(Xs[k] - self.eps * Y[k], lb=0, name=f"X_ge_eps_y__{k}"))
        mdl.add_cons_vars(cons)
        if self.ancla is not None:
            self.y[self.ancla].lb = 1
        # Tolerancias estrictas: con cotas de flujo de 1000, una binaria "casi cero"
        # (p.ej. y = 1e-5, dentro de la tolerancia de integralidad por defecto) deja
        # X = 1e-5 y flujos de hasta 1000 * 1e-5 = 0.01 > theta: el organismo trabaja
        # sin ser contado. Se ajustan las tolerancias y se verifica cada solución.
        try:
            tol = mdl.solver.configuration.tolerances
            tol.integrality = 1e-9
            tol.feasibility = 1e-9
        except Exception:
            pass

    # referencias por nombre (sobreviven a cambios de solver)
    @property
    def y(self) -> dict:
        return {k: self.com.model.variables[f"y__{k}"] for k in self.com.organismos}

    def ancla_funcional(self, reaccion_base: str, delta: float) -> None:
        """Exige v_{ancla, reaccion} >= delta * X_ancla (p.ej. reaccion_base='3HLYTCL')."""
        p = self.com.model.problem
        r = self.com.model.reactions.get_by_id(f"{reaccion_base}__{self.ancla}")
        c = p.Constraint(r.flux_expression - delta * self.com.X[self.ancla], lb=0,
                         name=f"ancla_funcional__{reaccion_base}")
        self.com.model.add_cons_vars([c])

    def quitar_ancla_funcional(self, reaccion_base: str) -> None:
        nombre = f"ancla_funcional__{reaccion_base}"
        if nombre in self.com.model.constraints:
            self.com.model.remove_cons_vars([self.com.model.constraints[nombre]])

    def resolver(self) -> list[str] | None:
        """Una comunidad mínima (lista de etiquetas) o None si es infactible."""
        acompanantes = [v for k, v in self.y.items() if k != self.ancla]
        val = self.com.optimizar(sum(acompanantes), "min")
        if val is None:
            return None
        return [k for k, v in self.y.items() if v.primal > 0.5]

    def verificar(self, comunidad: list[str]) -> bool:
        """
        Re-resuelve con las binarias FIJAS (y_k = 1 en la comunidad, 0 fuera, sin
        tolerancia posible) para confirmar que la comunidad realmente es factible.
        Detecta soluciones donde un organismo "apagado" aportaba flujo por fuga numérica.
        """
        Y = self.y
        respaldo = {k: (v.lb, v.ub) for k, v in Y.items()}
        try:
            for k, v in Y.items():
                b = 1 if k in comunidad else 0
                v.set_bounds(b, b)
            return self.com.factible()
        finally:
            for k, (lb, ub) in respaldo.items():
                Y[k].set_bounds(lb, ub)

    def cortar(self, comunidad: list[str]) -> None:
        """Corte de no superconjunto."""
        p = self.com.model.problem
        Y = self.y
        miembros = [k for k in comunidad if k != self.ancla]
        if not miembros:
            raise ValueError("La comunidad es solo el ancla: no hay alternativas que cortar")
        c = p.Constraint(sum(1 - Y[k] for k in miembros), lb=1,
                         name=f"corte_{len(self.cortes)}")
        self.com.model.add_cons_vars([c])
        self.cortes.append(c.name)

    def enumerar(self, max_soluciones: int = 50, tamano_max: int | None = None,
                 verbose: bool = True) -> list[list[str]]:
        """Todas las comunidades minimales (por inclusión), en orden de tamaño."""
        sols = []
        while len(sols) < max_soluciones:
            c = self.resolver()
            if c is None:
                break
            if tamano_max is not None and len(c) > tamano_max:
                break
            if not self.verificar(c):
                # Cortar aquí excluiría comunidades reales que contienen a c: se detiene.
                print(f"  AVISO: {c} no es factible con las binarias fijas (fuga numérica). "
                      f"Enumeración detenida; subir THETA o ajustar tolerancias.", flush=True)
                break
            sols.append(c)
            if verbose:
                print(f"  comunidad {len(sols)} (tamaño {len(c)}): {c}", flush=True)
            if len(c) == 1 and self.ancla is not None:   # solo el ancla: no hay más
                break
            self.cortar(c)
        return sols

    def limpiar_cortes(self) -> None:
        mdl = self.com.model
        mdl.remove_cons_vars([mdl.constraints[n] for n in self.cortes])
        self.cortes = []