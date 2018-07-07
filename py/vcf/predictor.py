from abc import ABCMeta, abstractmethod
from collections import defaultdict, Counter
from typing import Dict, List, TypeVar, Tuple, cast
from functools import lru_cache
from itertools import combinations

T = TypeVar('T')
V = TypeVar('V')

    
class Predictor(metaclass=ABCMeta):
    @abstractmethod
    def predict(self, genotype: float) -> Dict[str, float]:
        pass


class FrequencyPredictor(Predictor, metaclass=ABCMeta):
    
    def __init__(self, genotypes: List[float], labels: List[str]) -> None:
        counts: Dict[str, float] = defaultdict(float)
        totals: Dict[str, int] = defaultdict(int)
        
        for label, frequency in zip(labels, genotypes):
            counts[label] += frequency
            totals[label] += 1
        
        self.frequencies = {
            label: count / totals[label]
            for label, count in counts.items()
        }


class PredictorJoao(FrequencyPredictor):
    
    def __init__(self, genotypes: List[float], labels: List[str], threshold: float) -> None:
        super().__init__(genotypes, labels)
        self.threshold = threshold
    
    @lru_cache(maxsize=16)
    def predict(self, genotype: float) -> Dict[str, float]:
        distances = {
            label: abs(genotype - freq)
            for label, freq in self.frequencies.items()
        }
        ordered = sorted(distances, key=lambda x: distances[x])
        if distances[ordered[1]] - distances[ordered[0]] >= self.threshold:
            return {ordered[0]: 1.0}
        else:
            return {}


class PredictorMariana(FrequencyPredictor):
    
    def __init__(self, genotypes: List[float], labels: List[str], threshold: float) -> None:
        super().__init__(genotypes, labels)
        self.threshold = threshold
    
    @lru_cache(maxsize=16)
    def predict(self, genotype: float) -> Dict[str, float]:
        predictions = []
        pairs = combinations(self.frequencies, 2)
        for p1, p2 in pairs:
            d = abs(genotype - self.frequencies[p1]) - abs(genotype - self.frequencies[p2])
            if abs(d) >= self.threshold:
                if d < 0:
                    predictions.append(p1)
                else:
                    predictions.append(p2)
        
        counts = Counter(predictions)
        size = len(predictions)
        return {key: count / size for key, count in counts.items()}


class PredictorFactory(metaclass=ABCMeta):
    
    @abstractmethod
    def build(self, genotypes: List[float], labels: List[str]) -> Predictor:
        pass


class PredictorJoaoFactory(PredictorFactory):
    
    def __init__(self, threshold: float) -> None:
        super().__init__()
        self.threshold = threshold
    
    def build(self, genotypes: List[float], labels: List[str]) -> Predictor:
        return PredictorJoao(genotypes, labels, self.threshold)


class PredictorMarianaFactory(PredictorFactory):
    
    def __init__(self, threshold: float) -> None:
        super().__init__()
        self.threshold = threshold
    
    def build(self, genotypes: List[float], labels: List[str]) -> Predictor:
        return PredictorMariana(genotypes, labels, self.threshold)
