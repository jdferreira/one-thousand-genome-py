import random
import sys
from typing import Iterable, List, Sequence, TextIO, Tuple, TypeVar
from functools import lru_cache

from vcf.predictor import PredictorFactory, PredictorJoaoFactory
from vcf.population import Population, parse_population
from vcf.stream import VCFStream

__all__ = [
    'process_variant',
    'mean',
    'stdev',
    'split'
]

T = TypeVar('T')

@lru_cache()
def process_variant(variant: str) -> float:
    # Get the genotype of this individual
    genotype = variant.split(':')[0]
    
    if '|' in genotype:
        alleles = genotype.split('|')
    elif '/' in genotype:
        alleles = genotype.split('/')
    else:
        alleles = [genotype]
    
    alternatives = sum(1 for i in alleles if i != '0')
    return alternatives / len(alleles)


def mean(xs: Sequence[float]) -> float:
    return sum(xs) / len(xs)

def stdev(xs: Sequence[float]) -> float:
    m = mean(xs)
    return (sum((x - m) ** 2 for x in xs) / (len(xs) - 1)) ** 0.5


def split(to_split: Iterable[T], split_index: int) -> Tuple[List[T], List[T]]:
    clone = list(to_split)
    random.shuffle(clone)
    return clone[split_index:], clone[:split_index]
