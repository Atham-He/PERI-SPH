"""Keep one validated step, never deepcopy a complete physical trajectory."""
import copy,json
from decimal import localcontext,Decimal as D
from pathlib import Path
from peri_physical_time_friction_audit import PhysicalTimeFrictionInspector as FrictionInspector
from peri_elastic_normal import record


class OneStepCapture:
    def __init__(self):self.selected=None;self.score=-1.
    def __call__(self,audit,old,z,loss,history,row):
        score=max(abs(a-b) for a,b in zip(row['dense_minimum_plate_gap'],row['dense_endpoint_minimum_plate_gap']))
        if score<=self.score:return
        self.score=score
        self.selected=dict(case=audit.case,header=audit.h,old=old,z=z,loss=loss,history=history,row=row,
            material_cache=audit.cache.copy())


def negatives(capture,directory,source_raw):
    assert capture.selected is not None and capture.score>1e-10,'selected full-interval geometry must differ from endpoints'
    snapshot=capture.selected;out=[];directory=Path(directory);directory.mkdir()
    with localcontext() as ctx:
        ctx.prec=95
        for label in ['wrong_physical_mass','altered_tangent_dual_history','kinetic_using_wrong_mass',
            'swapped_material_algorithmic_losses','wrong_friction_loss','wrong_wall_velocity_cursor',
            'corrupt_compensation_history','endpoint_substituted_for_whole_interval','wrong_physical_tick_cursor','wrong_protocol_denominator']:
            h=snapshot['header'];row=copy.deepcopy(snapshot['row'])
            if label=='wrong_physical_tick_cursor':row['interval']['old_tick']+=1
            elif label=='wrong_protocol_denominator':row['interval']['time_denominator']+=1
            elif label=='wrong_physical_mass':
                h=copy.deepcopy(h);h['physical_mass'][0]*=2
            elif label=='altered_tangent_dual_history':
                row['friction_history']['last_tangent_impulse'][0][0]+=.001
            elif label=='kinetic_using_wrong_mass':row['kinetic_new']+=.001
            elif label=='swapped_material_algorithmic_losses':
                assert row['algorithmic_loss']!=row['evaluation']['physical_material_loss']
                row['algorithmic_loss'],row['evaluation']['physical_material_loss']=row['evaluation']['physical_material_loss'],row['algorithmic_loss']
            elif label=='wrong_friction_loss':row['friction_loss']+=.001
            elif label=='wrong_wall_velocity_cursor':row['interval']['wall_tangent_velocity'][1][0]+=.001
            elif label=='corrupt_compensation_history':row['friction_history']['relative_tangent_travel_correction'][0][0]+=.001
            else:row['dense_minimum_plate_gap']=row['dense_endpoint_minimum_plate_gap'].copy()
            path=directory/(label+'.json')
            artifact=dict(type='corrupted_single_step_context',label=label,source_raw=source_raw,
                selected_accepted_step=snapshot['row']['state']['step'],header=h,old_state=snapshot['old'],
                old_friction_history=snapshot['history'],old_reference_viscous=[[str(x) for x in z] for z in snapshot['z']],
                old_reference_loss=list(map(str,snapshot['loss'])),corrupted_record=row)
            with path.open('x') as stream:json.dump(artifact,stream,allow_nan=False)
            try:
                audit=FrictionInspector(snapshot['case'],h,True)
                audit.cache=snapshot['material_cache'].copy()
                audit.step(row,snapshot['old'],snapshot['z'],snapshot['loss'],snapshot['history'])
            except AssertionError:out.append(dict(name=label,rejected=True,artifact=record(path)))
            else:raise AssertionError(('negative accepted',label))
            del row,artifact
    return out
