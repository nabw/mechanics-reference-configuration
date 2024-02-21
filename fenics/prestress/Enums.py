class ControlTypes:
    distributed_only = 1
    surface_load_only = 2
    surface_pressure_only = 3
    distributed_and_surface_load = 4
    distributed_and_surface_pressure = 5


class SolverStrategy:
    dynamic = 1
    fixed = 3


class FormulationType:
    monolithic = 1
    split = 2


class MinSolvers:
    moolaNCG = 1
    moolaBFGS = 2
    moolaNewtonCG = 3
    scipyBFGS = 4


class InverseFormulationType:
    direct = 1
    sellier = 2
