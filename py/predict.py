import sys
from collections import defaultdict
from typing import Dict, Iterable, List, Sequence, TextIO, Tuple, TypeVar

from utils import process_variant, split
from vcf.population import Population, parse_population
from vcf.predictor import PredictorFactory, PredictorJoaoFactory
from vcf.stream import VCFStream


def predict(stream: VCFStream, population: Population, ratio: float, predictor_factory: PredictorFactory, min_dist: float):
    individuals = stream.individuals
    labels = population.get_labels(individuals)
    
    to_split = range(len(individuals))
    split_index = int(len(individuals) * ratio)
    train_idxs, target_idxs = split(to_split, split_index)
    
    cumulative_labels: List[Dict[str, float]] = [defaultdict(float) for _ in individuals]
    
    for polymorphism in stream:
        genotypes = [process_variant(i) for i in polymorphism.data_fields]
        
        predictor = predictor_factory.build(
            [genotypes[i] for i in train_idxs],
            [labels[i] for i in train_idxs]
        )
        
        for target_idx in target_idxs:
            genotype = genotypes[target_idx]
            predictions = predictor.predict(genotype)
            for prediction, confidence in predictions.items():
                cumulative_labels[target_idx][prediction] += confidence
        
    relative = 0 < min_dist < 1
    
    for target_idx in target_idxs:
        target_predictions = list(cumulative_labels[target_idx].items())
        target_predictions.sort(key=lambda x: x[1], reverse=True)
        
        diff = target_predictions[0][1] - target_predictions[1][1]
        if relative:
            diff /= target_predictions[0][1]
        
        if diff >= min_dist:
            final_prediction = target_predictions[0][0]
        else:
            final_prediction = '---'
        
        individual = individuals[target_idx]
        label = labels[target_idx]
        print(f'{individual}\t{label}\t{final_prediction}')
        

if __name__ == '__main__':
    predictor_factory = PredictorJoaoFactory(0)
    
    with open('populations/superpopulations.txt') as f:
        population = parse_population(f)
    
    stream = VCFStream.from_text_stream(sys.stdin)
    predict(stream, population, 0.1, predictor_factory, 0)
