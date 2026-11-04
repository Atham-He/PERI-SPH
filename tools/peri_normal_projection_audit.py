"""Independent exact-rational active-set oracle for bounded normal-only probes.

This does not call the C++ projected iteration or import its implementation.
It enumerates at most four active constraints and solves each equality system
with rational Gaussian elimination. It is not a scalable production solver.
"""
from fractions import Fraction as Q
from itertools import combinations
import math


def dot(a, b): return sum((x*y for x, y in zip(a, b)), Q(0))


def solve(matrix, rhs):
    n = len(rhs)
    a = [list(row) + [value] for row, value in zip(matrix, rhs)]
    for col in range(n):
        pivot = next((row for row in range(col, n) if a[row][col]), None)
        if pivot is None: return None
        a[col], a[pivot] = a[pivot], a[col]
        divisor = a[col][col]
        a[col] = [value / divisor for value in a[col]]
        for row in range(n):
            if row != col:
                multiplier = a[row][col]
                a[row] = [x - multiplier*y for x, y in zip(a[row], a[col])]
    return [row[-1] for row in a]


def exact_projection(case):
    mass, position, free = [[Q(str(x)) for x in case[key]] for key in ('mass', 'position', 'velocity')]
    dt = Q(str(case['dt']))
    if not mass or any(m <= 0 for m in mass) or len(mass) != len(position) or len(mass) != len(free) or dt <= 0:
        raise ValueError('invalid exact probe body')
    C, bound, gap, wall_velocity = [], [], [], []
    for row in case['rows']:
        coefficients = [Q(0)] * len(mass)
        for i, weight in zip(row['targets'], row['weights']): coefficients[i] = Q(str(row['normal'])) * Q(str(weight))
        current = dot(coefficients, position)
        old, new = Q(str(row['old_offset'])), Q(str(row['new_offset']))
        old_gap = current - old
        if old_gap < -dt * Q(str(case.get('velocity_scale', 1))) * Q(str(case.get('tolerance', 1e-12))):
            raise ValueError('initially penetrating exact probe')
        C.append(coefficients); bound.append((new-current)/dt); gap.append(old_gap); wall_velocity.append((new-old)/dt)
    if len(C) > 4 or len(mass) > 4: raise ValueError('exact oracle budget exceeded')
    for count in range(min(len(C), len(mass)) + 1):
        for active in combinations(range(len(C)), count):
            matrix = [[sum(C[i][k]*C[j][k]/mass[k] for k in range(len(mass))) for j in active] for i in active]
            rhs = [bound[i] - dot(C[i], free) for i in active]
            multipliers = solve(matrix, rhs)
            if multipliers is None or any(value < 0 for value in multipliers): continue
            impulse = [Q(0)] * len(C)
            for at, value in zip(active, multipliers): impulse[at] = value
            velocity = [free[k] + sum(C[i][k]*impulse[i] for i in range(len(C)))/mass[k] for k in range(len(mass))]
            if any(dot(row, velocity) < value for row, value in zip(C, bound)): continue
            kinetic_before = sum(m*v*v/2 for m, v in zip(mass, free))
            kinetic_after = sum(m*v*v/2 for m, v in zip(mass, velocity))
            work = dot(impulse, wall_velocity)
            dissipation = sum(l*g/dt for l, g in zip(impulse, gap)) + sum(m*(v-f)**2/2 for m, v, f in zip(mass, velocity, free))
            return dict(velocity=velocity, position=[q+dt*v for q, v in zip(position, velocity)], normal_impulse=impulse,
                        free_kinetic=kinetic_before, accepted_kinetic=kinetic_after, wall_work=work,
                        contact_step_dissipation=dissipation, C=C, bound=bound)
    return None


def check_native(case, native, code):
    expected_code = case.get('expected_code', 0)
    if code != expected_code: raise ValueError(f'Unexpected native return {code}, expected {expected_code}')
    if code == 1:
        try: exact_projection(case)
        except ValueError:
            if not native.get('error'): raise ValueError('Invalid state was not rejected explicitly')
            return dict(status='invalid_initial_geometry_rejected')
        raise ValueError('Native rejected an exact-valid probe')
    exact = exact_projection(case)
    if code == 2:
        if exact is not None or native['converged']: raise ValueError('Infeasible negative control incorrectly passed')
        return dict(status='infeasible_geometry_rejected', sweeps=native['sweeps'])
    if exact is None or not native['converged']: raise ValueError('Feasible probe failed')
    total_mass, scale_v = sum(case['mass']), case.get('velocity_scale', 1)
    scale_p, scale_e = total_mass*scale_v, total_mass*scale_v**2
    scale_q = max([abs(x) for x in case['position']] + [case['dt']*scale_v])
    tolerance = 1e-10
    errors = {}
    for name, scale in [('velocity', scale_v), ('position', scale_q)]:
        if len(native[name]) != len(exact[name]): raise ValueError('Native vector shape differs')
        error = max(abs(a - float(b))/scale for a, b in zip(native[name], exact[name]))
        if not math.isfinite(error) or error > tolerance: raise ValueError('Exact rational disagreement: ' + name)
        errors[name] = error
    for name in ('free_kinetic', 'accepted_kinetic', 'wall_work', 'contact_step_dissipation'):
        error = abs(native[name] - float(exact[name])) / max(scale_e, abs(float(exact[name])))
        if not math.isfinite(error) or error > tolerance: raise ValueError('Exact energy/work disagreement: ' + name)
        errors[name] = error
    impulse = native['normal_impulse']
    if len(impulse) != len(case['rows']) or any(not math.isfinite(value) or value < 0 for value in impulse):
        raise ValueError('Native dual impulse invalid')
    # Reconstruct KKT independently from emitted primal/dual values, without
    # assuming the dual distribution of redundant rows is unique.
    C = [[float(v) for v in row] for row in exact['C']]
    b = [float(v) for v in exact['bound']]
    primal = max([0.] + [value-sum(c*v for c, v in zip(row, native['velocity'])) for row, value in zip(C, b)])
    stationarity = max(abs(m*(v-f)-sum(row[k]*l for row, l in zip(C, impulse)))
                       for k, (m, v, f) in enumerate(zip(case['mass'], native['velocity'], case['velocity'])))
    complementarity = sum(abs(l*(sum(c*v for c, v in zip(row, native['velocity']))-value)) for row, l, value in zip(C, impulse, b))
    momentum = sum(m*(v-f) for m, v, f in zip(case['mass'], native['velocity'], case['velocity'])) - sum(row['normal']*l for row, l in zip(case['rows'], impulse))
    energy = native['wall_work'] - (native['accepted_kinetic']-native['free_kinetic']) - native['contact_step_dissipation']
    if primal > tolerance*scale_v or stationarity > tolerance*scale_p or complementarity > tolerance*scale_e or \
            abs(momentum) > tolerance*scale_p or abs(energy) > tolerance*scale_e or native['contact_step_dissipation'] < -tolerance*scale_e:
        raise ValueError('Independent KKT/force-transfer/energy check failed')
    return dict(status='exact_rational_primal_energy_and_native_dual_KKT_checked', normalized_errors=errors,
                primal_velocity_error=primal, stationarity_impulse_error=stationarity,
                complementarity_energy_error=complementarity, momentum_residual=momentum, energy_residual=energy)
