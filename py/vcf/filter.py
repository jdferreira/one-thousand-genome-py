from abc import ABCMeta, abstractmethod
from enum import Enum, IntEnum, auto
from typing import Dict, Iterable, Iterator, List, TextIO, Tuple, Union

from vcf.stream import VCFData, VCFDataAction, VCFStream


class VCFFilter(metaclass=ABCMeta):
    def pipe(self, stream: VCFStream) -> VCFStream:
        return _StreamFromFilter(stream, self)
    
    def filter_metadata(self, metadata: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
        return metadata
    
    def filter_individuals(self, individuals: List[str]) -> List[str]:
        return individuals
    
    def filter_item(self, item: VCFData) -> Union[VCFDataAction, VCFData]:
        return item


class _StreamFromFilter(VCFStream):
    def __init__(self, inner: VCFStream, vcf_filter: VCFFilter) -> None:
        self.inner = inner
        self.filter = vcf_filter
        super().__init__()
    
    def start(self) -> None:
        self.metadata = self.filter.filter_metadata(self.inner.metadata)
        self.individuals = self.filter.filter_individuals(self.inner.individuals)
    
    def __next__(self) -> VCFData:
        while True:
            item = next(self.inner)
            action = self.filter.filter_item(item)
            
            if isinstance(action, VCFData):
                return action
            
            if action == VCFDataAction.STOP:
                raise StopIteration


class IndividualsFilter(VCFFilter):
    def __init__(self, individuals: Iterable[str], strict: bool = False) -> None:
        self.individuals = set(individuals)
        self.indices: List[int] = []

        self.strict = strict
    
    def filter_individuals(self, individuals: List[str]) -> List[str]:
        # For each individual in the VCF file, determine whether we are
        # interested in it and keep it only if so
        filtered_individuals = []

        for idx, individual in enumerate(individuals):
            if individual in self.individuals:
                self.indices.append(idx)
                filtered_individuals.append(individual)

        # Also make sure that all individuals of interest have been found
        if self.strict and set(filtered_individuals) != self.individuals:
            bad_individuals = [i for i in self.individuals if i not in filtered_individuals]
            raise Exception(
                f'Cannot find the following individuals in this file: '
                f'{bad_individuals!r}'
            )

        return filtered_individuals

    def filter_item(self, item: VCFData) -> VCFData:
        filtered_fields = [item.data_fields[i] for i in self.indices]
        return VCFData(item.info_fields, filtered_fields)


class PolymorphismFilter(VCFFilter):
    def __init__(self, polymorphisms: Iterable[str]) -> None:
        self.polymorphisms = set(polymorphisms)
        self.count = 0

    def filter_item(self, item: VCFData) -> Union[VCFDataAction, VCFData]:
        if self.count == len(self.polymorphisms):
            return VCFDataAction.STOP
        elif item.get_identifier() not in self.polymorphisms:
            return VCFDataAction.IGNORE
        else:
            self.count += 1
            return item
