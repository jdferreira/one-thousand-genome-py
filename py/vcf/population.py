from collections import defaultdict
from typing import DefaultDict, Dict, Iterable, List, TextIO, Tuple


class Population:
    """
    This class represents a population, which is a list of individuals along
    with the group they belong to
    """
    
    def __init__(self) -> None:
        self.individual_to_group: Dict[str, str] = {}
        self.group_to_individuals: DefaultDict[str, List[str]] = defaultdict(list)
    
    
    def add_individual(self, identifier: str, group: str) -> None:
        """
        Add a new individual associated with the given group. Individual
        identifiers must be unique and can only belong to one group.
        """
        
        self.individual_to_group[identifier] = group
        self.group_to_individuals[group].append(identifier)
    
    
    def remove_individual(self, identifier: str) -> None:
        if identifier in self.individual_to_group:
            group = self.individual_to_group[identifier]
            del self.individual_to_group[identifier]
            self.group_to_individuals[group].remove(identifier)
    
    
    def individuals(self) -> Iterable[str]:
        """
        Returns an interator of the individuals in this population
        """
        
        return iter(self.individual_to_group)
    
    
    def groups(self) -> Iterable[str]:
        """
        Returns an interator of the groups in this population
        """
        
        return iter(self.group_to_individuals)
    
    
    def items(self) -> Iterable[Tuple[str, List[str]]]:
        """
        Return an item-like iterator that yields the various groups in this
        population. Each item is of the form
            (group_name, list_of_individuals)
        """
        
        return self.group_to_individuals.items()
    
    
    def has_individual(self, individual: str) -> bool:
        """
        Determines whether this population contains a specific individual.
        """
        
        return individual in self.individual_to_group


def parse_population(stream: TextIO) -> Population:
    """
    Read a text stream, where each non-empty line has two space-separated
    fields: an individual identifier and the group they belong to.
    """
    
    result = Population()
    
    for line in stream:
        assert isinstance(line, str)
        
        line = line.rstrip('\n')
        
        # Remove possible comments
        comment_start = line.find('#')
        if comment_start >= 0:
            line = line[:comment_start]
        line = line.strip()
        
        if not line:
            # Ignore empty lines
            continue
        
        fields = line.split()
        if len(fields) != 2:
            raise Exception(
                f'Line {line!r} is invalid: needs 2 fields, found {len(fields)}'
            )
        
        identifier, group = fields
        result.add_individual(identifier, group)
    
    return result
