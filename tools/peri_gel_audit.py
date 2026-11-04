#!/usr/bin/env python3
"""Independent stored-frame gel energy/stress/contact reconstruction (stdlib).

This checks every available full-precision VTP frame. It does not claim that
sparse frames provide independent raw energy verification of every time step.
"""
import csv
import hashlib
import math
from pathlib import Path
import struct
import xml.etree.ElementTree as ET


def checkpoint_physical_digest(file, output_directory):
    """Independently identify/hash the full physical byte region, excluding the
    header's wall-clock field and output attachments. Validate both FNV footers
    and every attached output prefix. No particle/history fields are omitted."""
    data = file.read_bytes()
    def fnv(values):
        value = 14695981039346656037
        for byte in values: value = ((value ^ byte) * 1099511628211) & ((1 << 64) - 1)
        return value
    if len(data) < 16: raise ValueError('truncated checkpoint')
    checksum, footer = struct.unpack_from('<QQ', data, len(data) - 16)
    if footer != 0x50455249444F4E45 or fnv(data[:-16]) != checksum:
        raise ValueError('independent complete checkpoint checksum failed')
    at = 0
    def number():
        nonlocal at
        value, = struct.unpack_from('<Q', data, at); at += 8
        return value
    def label():
        nonlocal at
        length = number()
        if length > 1 << 20 or at + length > len(data) - 16: raise ValueError('invalid checkpoint label')
        value = data[at:at+length].decode(); at += length
        return value
    if label() != 'PERI-SPH-CPU-CHECKPOINT-1': raise ValueError('unknown checkpoint format')
    if [number() for _ in range(4)] != [0x0102030405060708, 8, 8, 3]:
        raise ValueError('checkpoint ABI differs from the qualified CPU format')
    label() # Exact executable/config identity remains included in the full-file SHA.
    progress = struct.unpack_from('<ddQQQ', data, at); at += 40
    start = at
    names = ['effective_parameters.cfg', 'diagnostics.csv', 'step_trace.csv']
    marker = struct.pack('<QQ', len(names), len(names[0])) + names[0].encode()
    end = data.find(marker, start)
    if end < start or data.find(marker, end + 1) != -1: raise ValueError('ambiguous checkpoint physical-region boundary')
    at = end
    if number() != len(names): raise ValueError('checkpoint output count differs')
    for expected in names:
        if label() != expected: raise ValueError('checkpoint output schema differs')
        size = number()
        if size > 1 << 32 or at + size > len(data) - 16: raise ValueError('checkpoint output prefix invalid')
        with (output_directory / expected).open('rb') as stream:
            prefix = stream.read(size)
        if prefix != data[at:at+size]: raise ValueError('checkpoint output prefix differs: ' + expected)
        at += size
    if at != len(data) - 16: raise ValueError('unexpected checkpoint trailing payload')
    return dict(physical_sha256=hashlib.sha256(data[start:end]).hexdigest(),
                physical_fnv1a64=format(fnv(data[start:end]), 'x'),
                time_s=progress[0], acoustic_steps=progress[2], pd_steps=progress[3], completed_frame=progress[4])


def dot(a, b): return math.fsum(x * y for x, y in zip(a, b))
def transpose(a): return [list(row) for row in zip(*a)]
def multiply(a, b): return [[dot(row, col) for col in zip(*b)] for row in a]
def matrix(flat): return [[flat[i + 3 * j] for j in range(3)] for i in range(3)]


def determinant(F):
    a, b, c = F
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])


def material_response(F, viscous, config):
    """Different code path from C++: direct Cauchy tensor and scalar energy."""
    young, nu = config['youngs_modulus'], config['poisson_ratio']
    bulk, shear = young / (3 * (1 - 2 * nu)), young / (2 * (1 + nu))
    branch = shear * config['relaxing_shear_ratio']
    J = determinant(F)
    if J < 0.05 or not math.isfinite(J): raise ValueError('raw accepted deformation invalid')
    C = multiply(transpose(F), F)
    B = multiply(F, transpose(F))
    trace = math.fsum(C[d][d] for d in range(3))
    green = [[0.5 * (C[i][j] - (trace / 3 if i == j else 0)) for j in range(3)] for i in range(3)]
    q = [[green[i][j] - viscous[i][j] for j in range(3)] for i in range(3)]
    spatial_q = multiply(multiply(F, q), transpose(F))
    stress = [[shear * J**(-5 / 3) * (B[i][j] - (trace / 3 if i == j else 0)) +
               (0.5 * bulk * (J - 1 / J) if i == j else 0) + 2 * branch * spatial_q[i][j] / J
               for j in range(3)] for i in range(3)]
    energy = shear / 2 * (J**(-2 / 3) * trace - 3) + bulk / 4 * (J * J - 1 - 2 * math.log(J)) + branch * math.fsum(x * x for row in q for x in row)
    return energy, stress, green


def box_response(position, velocity, center, half, wall_velocity, volume, config):
    local = [x - y for x, y in zip(position, center)]
    q = [abs(x) - h for x, h in zip(local, half)]
    outside = [max(value, 0) for value in q]
    length = math.sqrt(dot(outside, outside))
    normal = [0., 0., 0.]
    if length:
        normal = [(1 if x >= 0 else -1) * y / length for x, y in zip(local, outside)]
        distance = length
    else:
        distance = max(q); axis = q.index(distance)
        normal[axis] = 1 if local[axis] >= 0 else -1
    spacing = config['spacing']
    penetration = max(spacing / 2 - distance, 0)
    young, nu = config['youngs_modulus'], config['poisson_ratio']
    penalty = young / (3 * (1 - 2 * nu)) + 4 / 3 * young / (2 * (1 + nu)) * (1 + config['relaxing_shear_ratio'])
    stiffness = penalty * volume / spacing**2
    relative = [v - w for v, w in zip(velocity, wall_velocity)]
    tangent = [v - dot(relative, normal) * n for v, n in zip(relative, normal)]
    speed = math.sqrt(dot(tangent, tangent))
    force = [stiffness * penetration * n for n in normal]
    if speed and penetration:
        friction = config['friction_mu'] * stiffness * penetration * math.tanh(speed / config['friction_velocity']) / speed
        force = [f - friction * t for f, t in zip(force, tangent)]
    return force, stiffness * penetration**2 / 2, penetration


def audit_frames(directory, config):
    with (directory / 'output/diagnostics.csv').open() as stream:
        diagnostics = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    with (directory / 'output/step_trace.csv').open() as stream:
        trace = {int(row['pd_steps']): float(row['dt_s']) for row in csv.DictReader(stream)}
    files = sorted((directory / 'output').glob('GelParticle_ite_*.vtp'))
    if len(files) != len(diagnostics): raise ValueError('full accepted frames are missing or duplicated')
    young, spacing = config['youngs_modulus'], config['spacing']
    duration = config['compression_stroke'] / config['compression_speed']
    half = [config['plate_half_width'], config['plate_thickness'] / 2, config['plate_half_width']]
    error = dict(energy_over_E=0., stress_over_E=0., green=0., contact_force_over_Edp2=0., density_over_rho=0.)
    initial_masses = None
    previous_loss = None
    checked = 0
    for file, row in zip(files, diagnostics):
        fields = {node.get('Name'): [float(value) for value in (node.text or '').split()]
                  for node in ET.parse(file).iter('DataArray') if node.get('Name')}
        count = int(row['solid_count'])
        if fields['TimeValue'] != [row['time_s']]: raise ValueError('VTP time differs from accepted diagnostics')
        if any(not math.isfinite(value) for name, values in fields.items() for value in values):
            raise ValueError('non-finite full raw field')
        masses = fields['Mass']
        if len(masses) != count: raise ValueError('raw particle count changed')
        if initial_masses is None: initial_masses = masses
        if initial_masses != masses: raise ValueError('raw particle masses changed')
        if previous_loss is not None and any(b < a for a, b in zip(previous_loss, fields['GelDissipationDensity'])):
            raise ValueError('a particle lost accumulated viscous dissipation')
        previous_loss = fields['GelDissipationDensity']
        energy_sum, loss_sum, observed_upper, applied_upper, potential = [], [], [], [], []
        for i in range(count):
            F = matrix(fields['DeformationGradient'][9*i:9*i+9])
            error['density_over_rho'] = max(error['density_over_rho'],
                abs(fields['Density'][i] / config['density'] - 1 / determinant(F)))
            if not math.isclose(fields['GelHistoryTime'][i], row['time_s'], rel_tol=1e-12, abs_tol=1e-20):
                raise ValueError('raw accepted history has a different physical time')
            viscous = matrix(fields['GelViscousStrain'][9*i:9*i+9])
            energy, stress, green = material_response(F, viscous, config)
            native_green = matrix(fields['GelGreenStrain'][9*i:9*i+9])
            native_stress = matrix(fields['GelAcceptedStress'][9*i:9*i+9])
            error['energy_over_E'] = max(error['energy_over_E'], abs(energy - fields['GelEnergyDensity'][i]) / young)
            error['stress_over_E'] = max(error['stress_over_E'], max(abs(stress[a][b] - native_stress[a][b]) for a in range(3) for b in range(3)) / young)
            error['green'] = max(error['green'], max(abs(green[a][b] - native_green[a][b]) for a in range(3) for b in range(3)))
            volume = fields['VolumetricMeasure'][i]
            energy_sum.append(energy * volume); loss_sum.append(fields['GelDissipationDensity'][i] * volume)
            if fields['GelDissipationDensity'][i] < 0: raise ValueError('negative raw gel dissipation')
            position, velocity = fields['Position'][3*i:3*i+3], fields['Velocity'][3*i:3*i+3]
            upper = [0, row['upper_face_y_m'] + half[1], 0]
            wall_velocity = [0, -config['compression_speed'] if row['time_s'] < duration else 0, 0]
            force, U, _ = box_response(position, velocity, upper, half, wall_velocity, volume, config)
            observed_upper.append(-force[1]); potential.append(U)
            _, lower_U, _ = box_response(position, velocity, [0, -half[1], 0], half, [0, 0, 0], volume, config)
            potential.append(lower_U)
            if row['pd_steps']:
                dt = trace[int(row['pd_steps'])]
                for d in range(3):
                    predicted = fields['GelContactBaseVelocity'][3*i+d] + 0.5 * dt * fields['PDAccelerationOld'][3*i+d]
                    if abs(predicted - fields['GelContactSampleVelocity'][3*i+d]) > 1e-12 * (1 + abs(predicted)):
                        raise ValueError('raw contact midpoint velocity predictor differs')
                sample_time = row['last_contact_sample_time_s']
                sample_center = [0, 2 * (config['radius'] + config['initial_gap']) + half[1] - min(sample_time / duration, 1) * config['compression_stroke'], 0]
                sample_force, _, _ = box_response(fields['GelContactSamplePosition'][3*i:3*i+3],
                    fields['GelContactSampleVelocity'][3*i:3*i+3], sample_center, half,
                    [0, -config['compression_speed'] if sample_time < duration else 0, 0], volume, config)
                reaction = fields['UpperPlateReaction'][3*i:3*i+3]
                error['contact_force_over_Edp2'] = max(error['contact_force_over_Edp2'],
                    max(abs(reaction[d] + sample_force[d]) for d in range(3)) / (young * spacing**2))
                applied_upper.append(reaction[1])
        energy_scale = young * math.fsum(fields['VolumetricMeasure'])
        force_scale = young * spacing**2 * count
        checks = [('gel_stored_energy_J', math.fsum(energy_sum), energy_scale),
                  ('gel_dissipation_J', math.fsum(loss_sum), energy_scale),
                  ('observed_contact_potential_J', math.fsum(potential), energy_scale),
                  ('upper_observed_Fy_N', math.fsum(observed_upper), force_scale)]
        if row['pd_steps']: checks.append(('upper_last_applied_Fy_N', math.fsum(applied_upper), force_scale))
        for name, independent, scale in checks:
            if abs(row[name] - independent) > 1e-10 * scale:
                raise ValueError(f'{file.name}: independent {name} reconstruction failed')
        checked += count
    if any(value > 1e-10 for value in error.values()): raise ValueError('independent full-raw material/contact arithmetic failed: ' + str(error))
    return dict(stored_frames=len(files), particle_states_reconstructed=checked, maximum_normalized_errors=error,
                scope='every stored full-precision frame; sparse raw is not independent full-step energy')
