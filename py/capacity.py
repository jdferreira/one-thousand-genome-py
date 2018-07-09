import sys
from typing import Iterable, List, Sequence, TextIO, Tuple, TypeVar

from utils import process_variant, mean, stdev, split
from vcf.population import Population, parse_population
from vcf.predictor import PredictorFactory, PredictorJoaoFactory
from vcf.stream import VCFStream


def compute_capacity(stream: VCFStream, population: Population, repeats: int, ratio: float, predictor_factory: PredictorFactory):
    individuals = stream.individuals
    split_index = int(len(individuals) * ratio)
    
    labels = population.get_labels(individuals)
    
    # Split once, since the call to `split` takes a lot of time to do for every
    # polymorphism and repeat. We call for the number of repeats and apply the
    # corresponding one, copying the same approach for all polymorphisms.
    if repeats > 1:
        train_target_indices_list = [
            split(range(len(individuals)), split_index)
            for _ in range(repeats)
        ]
        train_target_indices = [
            (set(i), set(j))
            for i, j in train_target_indices_list
        ]
    else:
        train_target_indices = [
            (set(range(i)), set(range(i)))
            for i, _ in enumerate(individuals)
        ]
    
    for polymorphism in stream:
        accuracies: List[float] = []
        
        for repeat_idx in range(repeats):
            genotypes = (process_variant(i) for i in polymorphism.data_fields)
            to_split = zip(individuals, genotypes, labels)
            
            train = []
            target = []
            for i, val in enumerate(to_split):
                if i in train_target_indices[repeat_idx][0]:
                    train.append(val)
                if i in train_target_indices[repeat_idx][1]:
                    target.append(val)
            
            _train_individuals, train_genotypes, train_labels = zip(*train)
            
            predictor = predictor_factory.build(train_genotypes, train_labels)
            
            count = 0.0
            for _target_individual, genotype, label in target:
                predictions = predictor.predict(genotype)
                if label in predictions:
                    count += predictions[label]
            
            accuracies.append(count / len(target))
        
        m = mean(accuracies)
        if repeats == 1:
            print(f'{polymorphism.get_identifier()}\t{m:.5f}')
        else:
            s = stdev(accuracies)
            print(f'{polymorphism.get_identifier()}\t{m:.5f}\t{s:.5f}')
        

if __name__ == '__main__':
    predictor_factory = PredictorJoaoFactory(0)
    
    with open('populations/superpopulations.txt') as f:
        population = parse_population(f)
    
    stream = VCFStream.from_text_stream(sys.stdin)
    compute_capacity(stream, population, 20, 0.1, predictor_factory)
