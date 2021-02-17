"""
The :mod:`varannot` module defines functions for annotating variants using a combination of manual, vep and
vcf2maf annotations.
"""

from ._main import run_annotator, Vcf2mafConfig, VepConfig
from ._manual import run_manual_annotator
from ._vcf2maf import run_vcf2maf_annotator
from ._vep import run_vep_annotator

__all__ = [
    'run_annotator',
    'Vcf2mafConfig',
    'VepConfig',
    'run_manual_annotator',
    'run_vcf2maf_annotator',
    'run_vep_annotator'
]

