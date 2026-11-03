#!/usr/bin/env python3
"""Read capsule experiments without sorting away loading/unloading history.

The workbook is read-only. Exported row order is source order, not physical
time. No missing diameter, loading rate, contact zero or time is inferred.
"""
import argparse
import csv
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
import re
import struct


UNIT_SCALE = {'um': 1e-6, 'nm': 1e-9, 'mN': 1e-3, 'uN': 1e-6}
SHEET_UNITS = {'Micro-H': ('um', 'mN'), 'Micro-D': ('um', 'mN'),
               'Nano-H': ('nm', 'uN'), 'Nano-D': ('nm', 'uN')}


def normalized_unit(value):
    return str(value).strip().replace('\u03bc', 'u').replace('\u00b5', 'u')


def sample_metadata(label, sheet):
    match = re.fullmatch(r'(\d+(?:\.\d+)?)TO-(\d+(?:\.\d+)?)SA-([HD])', label)
    if not match or not sheet.endswith('-' + match[3]):
        raise ValueError(f'{sheet}: invalid or inconsistent sample label {label!r}')
    oil, alginate = float(match[1]), float(match[2])
    return {
        'sample_label': label,
        'nominal_tung_oil_loading_percent': oil,
        'loading_fraction_scope': 'synthesis label; not final oil volume fraction',
        'nominal_alginate_solution_concentration_percent': alginate,
        'hydration_state': 'hydrated' if match[3] == 'H' else 'dehydrated',
        'contains_tung_oil_by_formulation': oil > 0,
        'specimen_diameter_m': None,
        'shell_thickness_m': None,
        'loading_rate_m_s': None,
        'physical_time_available': False,
        'contact_zero_verified': False,
        'replicate_identity': None,
    }


@dataclass
class Curve:
    sheet: str
    label: str
    displacement_unit: str
    force_unit: str
    points: list = field(default_factory=list)  # (source row, raw displacement, raw force)

    def point_digest(self):
        digest = hashlib.sha256()
        for _, displacement, force in self.points:
            digest.update(struct.pack('<dd', displacement * UNIT_SCALE[self.displacement_unit],
                                      force * UNIT_SCALE[self.force_unit]))
        return digest.hexdigest()


def read_sheet_rows(sheet, rows):
    """Consume a worksheet once, including all labelled curve pairs."""
    if sheet not in SHEET_UNITS:
        raise ValueError(f'Unsupported experimental sheet {sheet!r}')
    rows = iter(rows)
    try:
        labels, headings, units = next(rows), next(rows), next(rows)
    except StopIteration as error:
        raise ValueError(f'{sheet}: three header rows are required') from error
    curves = {}
    labels_seen = set()
    for column, label in enumerate(labels):
        if label is None or not str(label).strip():
            continue
        label = str(label).strip()
        sample_metadata(label, sheet)
        if label in labels_seen or column + 1 >= len(units) or column + 1 >= len(headings):
            raise ValueError(f'{sheet}: duplicate label or missing force column: {label}')
        labels_seen.add(label)
        if str(headings[column]).strip().lower() != 'displacement' or str(headings[column + 1]).strip().lower() != 'force':
            raise ValueError(f'{sheet}: unexpected column headings for {label}')
        dunit, funit = normalized_unit(units[column]), normalized_unit(units[column + 1])
        if (dunit, funit) != SHEET_UNITS[sheet]:
            raise ValueError(f'{sheet}: unexpected units {(dunit, funit)} for {label}')
        curves[column] = Curve(sheet, label, dunit, funit)
    if not curves:
        raise ValueError(f'{sheet}: no labelled curves found')
    for row_number, row in enumerate(rows, 4):
        for column, curve in curves.items():
            displacement = row[column] if column < len(row) else None
            force = row[column + 1] if column + 1 < len(row) else None
            if displacement is None and force is None:
                continue
            try:
                if isinstance(displacement, bool) or isinstance(force, bool):
                    raise ValueError('boolean is not a measurement')
                displacement, force = float(displacement), float(force)
                if not math.isfinite(displacement) or not math.isfinite(force):
                    raise ValueError('non-finite measurement')
            except (TypeError, ValueError) as error:
                raise ValueError(f'{sheet} row {row_number}, {curve.label}: invalid measurement pair') from error
            curve.points.append((row_number, displacement, force))
    for curve in curves.values():
        if not curve.points:
            raise ValueError(f'{sheet}: empty curve {curve.label}')
    return list(curves.values())


def duplicate_groups(curves):
    groups = {}
    for curve in curves:
        # Equal SI-valued curves are flagged across labels; they are never deleted.
        groups.setdefault(curve.point_digest(), []).append(f'{curve.sheet}_{curve.label}')
    return [labels for labels in groups.values() if len(labels) > 1]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def export_workbook(workbook, destination):
    # Imported only by this read-only workbook command. Core/unit checks use stdlib.
    import openpyxl
    workbook = Path(workbook).resolve(strict=True)
    destination = Path(destination).resolve()
    if destination.exists():
        raise FileExistsError(f'Refusing to replace an existing experiment export: {destination}')
    source_sha = sha(workbook)
    # Measurements must be literal numeric cells; do not trust stale formula caches.
    book = openpyxl.load_workbook(workbook, read_only=True, data_only=False)
    try:
        missing = set(SHEET_UNITS) - set(book.sheetnames)
        if missing:
            raise ValueError(f'Missing experimental sheets: {sorted(missing)}')
        curves = [curve for sheet in SHEET_UNITS for curve in read_sheet_rows(sheet, book[sheet].values)]
    finally:
        book.close()
    destination.mkdir(parents=True)
    records = []
    for curve in curves:
        file = destination / f'{curve.sheet}_{curve.label}.csv'
        with file.open('w', newline='', encoding='utf-8') as stream:
            writer = csv.writer(stream)
            writer.writerow(['source_row', 'displacement_source', 'force_source',
                             'displacement_m', 'force_N', 'time_s'])
            for row, displacement, force in curve.points:
                writer.writerow([row, format(displacement, '.17g'), format(force, '.17g'),
                                 format(displacement * UNIT_SCALE[curve.displacement_unit], '.17g'),
                                 format(force * UNIT_SCALE[curve.force_unit], '.17g'), ''])
        records.append({
            **sample_metadata(curve.label, curve.sheet),
            'sheet': curve.sheet, 'experiment_kind': 'compression' if curve.sheet.startswith('Micro') else 'nanoindentation',
            'source_units': [curve.displacement_unit, curve.force_unit],
            'source_order_preserved': True, 'rows': len(curve.points),
            'decreasing_displacement_transitions': sum(b[1] < a[1] for a, b in zip(curve.points, curve.points[1:])),
            'sequence_scope': 'worksheet row order; not an inferred physical time axis',
            'SI_point_sha256_in_source_order': curve.point_digest(),
            'csv': file.name, 'csv_sha256': sha(file),
        })
    if sha(workbook) != source_sha:
        raise RuntimeError('Source workbook changed during extraction')
    manifest = {'schema_version': 1, 'source_workbook': str(workbook), 'source_sha256': source_sha,
                'source_modified': False, 'curve_count': len(curves),
                'total_rows': sum(len(c.points) for c in curves),
                'duplicate_groups': duplicate_groups(curves), 'curves': records,
                'material_calibration_performed': False,
                'limits': ['No time axis or loading rate inferred from row number.',
                           'No sorting, averaging, zero shift or force scaling applied.',
                           'No specimen geometry inferred from a group particle-size distribution.',
                           'Equal curves with different labels are not independent validation data.']}
    (destination / 'manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest = export_workbook(args.workbook, args.output)
    print(json.dumps({key: manifest[key] for key in ('curve_count', 'total_rows', 'duplicate_groups', 'material_calibration_performed')}, indent=2))


if __name__ == '__main__':
    main()
