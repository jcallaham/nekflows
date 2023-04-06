import os

from .flows import(
    NekFlowConfig,
    LidDrivenCavity,
    ShearDrivenCavity,
    MixingLayer,
    CylinderWake
)

from .nek import(
    NekVector,
    NekHandle,
    ComplexNekHandle,
    mean,
    project
)

import os
default_case_dir = os.path.abspath(f'{__file__}/../../cases')

def cyl(**kwargs):
    return CylinderWake(field_path=f'{default_case_dir}/cyl/pod_modes')

def mix(**kwargs):
    return MixingLayer(field_path=f'{default_case_dir}/mix/long_domain/pod_modes')

def short_mix(**kwargs):
    return MixingLayer(field_path=f'{default_case_dir}/mix/short_domain/pod_modes', Lx=60)

def lid2d(**kwargs):
    return LidDrivenCavity(field_path=f'{default_case_dir}/lid2d/pod_modes',
        base_path=f'{default_case_dir}/lid2d/base_flows')

def shear2d(**kwargs):
    return ShearDrivenCavity(field_path=f'{default_case_dir}/shear2d/global/pod_modes',
        base_path=f'{default_case_dir}/shear2d/global/base_flows')

def cavity7500(**kwargs):
    return ShearDrivenCavity(field_path=f'{default_case_dir}/shear2d/7500/pod_modes')

def cavity10k(**kwargs):
    return ShearDrivenCavity(field_path=f'{default_case_dir}/shear2d/10k/pod_modes')