from abc import ABCMeta, abstractmethod
from enum import Enum, IntEnum, auto
from typing import Dict, Iterable, Iterator, List, TextIO, Tuple, Union


class VCFColumn(IntEnum):
    CHROMOSOME = 0
    POS = 1
    ID = 2
    REF = 3
    ALT = 4
    QUAL = 5
    FILTER = 6
    INFO = 7
    FORMAT = 8


class VCFData:
    def __init__(self, info_fields: List[str], data_fields: List[str]) -> None:
        self.info_fields = info_fields
        self.data_fields = data_fields
                                                                                     
    def get_chromosome(self) -> str:
        return self.info_fields[VCFColumn.CHROMOSOME]
                                                                                     
    def get_position(self) -> str:
        return self.info_fields[VCFColumn.POS]
                                                                                     
    def get_identifier(self) -> str:
        return self.info_fields[VCFColumn.ID]
                                                                                     
    def get_reference(self) -> str:
        return self.info_fields[VCFColumn.REF]
                                                                                     
    def get_alternatives(self) -> List[str]:
        return self.info_fields[VCFColumn.ALT].split(',')
                                                                                     
    def get_quality(self) -> str:
        return self.info_fields[VCFColumn.QUAL]
    
    def get_filter(self) -> str:
        return self.info_fields[VCFColumn.FILTER]
    
    def get_info(self) -> str:
        return self.info_fields[VCFColumn.INFO]
    
    def get_format(self) -> str:
        return self.info_fields[VCFColumn.FORMAT]
    
    def get_genotype(self, index: int) -> str:
        return self.data_fields[index]


class VCFDataAction(Enum):
    IGNORE = auto()
    STOP = auto()


class VCFStream(Iterable[VCFData], metaclass=ABCMeta):
    def __init__(self) -> None:
        self.metadata: List[Tuple[str, str]] = []
        self.individuals: List[str] = []
    
    def __iter__(self) -> Iterator[VCFData]:
        return self
    
    @abstractmethod
    def __next__(self) -> VCFData:
        pass
    
    @classmethod
    def from_text_stream(cls, text_stream: TextIO) -> 'VCFStream':
        return _StreamFromFile(text_stream)


class _StreamFromFile(VCFStream):
    
    def __init__(self, text_stream: TextIO) -> None:
        super().__init__()
        
        self.text_stream = text_stream
        self.genotype: bool
        
        self.start()
    
    
    def start(self) -> None:
        # Read all metadata lines and the header line, setting the text stream
        # to the point where only data follows. The metadata and individuals
        # list is stored within this VCFStream object

        while True:
            try:
                line = next(self.text_stream)
            except StopIteration:
                if len(self.metadata) == 0:
                    msg = 'VCF file is empty'
                else:
                    msg = 'VCF file contains only metadata'
                raise Exception(msg)

            line = line.rstrip('\n')

            if line.startswith('##'):
                key, value = line[2:].split('=', 1)
                self.metadata.append((key, value))

            elif line.startswith('#'):
                if len(self.metadata) == 0:
                    raise Exception('VCF file does not contain metadata')

                fields = line.split()

                if 'FORMAT' in fields:
                    self.genotype = True
                    start_index = VCFColumn.FORMAT + 1
                else:
                    start_index = VCFColumn.FORMAT
                self.individuals = fields[start_index:]
                break

            else:
                raise Exception('VCF file does not contain a header line')
    
    
    def __next__(self) -> VCFData:
        line = next(self.text_stream).rstrip('\n')
        fields = line.split()

        if self.genotype:
            split_index = VCFColumn.FORMAT + 1
        else:
            split_index = VCFColumn.FORMAT
        info, data = fields[:split_index], fields[split_index:]

        # If no FORMAT column is given, append an empty string as if one had
        # been given
        if not self.genotype:
            info.append('')

        return VCFData(info, data)
