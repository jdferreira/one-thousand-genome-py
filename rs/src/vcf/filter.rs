use std::collections::HashSet;
use vcf::stream::{DataStream, IndividualsReader, MetadataReader, VCFData};

pub enum DataAction {
    Ignore,
    Stop,
    Data(VCFData),
}

pub trait Filter {
    fn filter_metadata(&mut self, metadata: Vec<(String, String)>) -> Vec<(String, String)> {
        metadata
    }

    fn filter_individuals(&mut self, individuals: Vec<String>) -> Vec<String> {
        individuals
    }

    fn filter_item(&mut self, item: VCFData) -> DataAction {
        DataAction::Data(item)
    }
}

pub struct Pipe<S: MetadataReader>(S);

impl<S: MetadataReader> Pipe<S> {
    pub fn new(inner: S) -> Pipe<S> {
        Pipe(inner)
    }

    pub fn pipe<F: Filter>(self, filter: F) -> Pipe<impl MetadataReader> {
        Pipe(MetadataReaderFromFilter {
            inner: self.0,
            filter,
        })
    }
    
    pub fn stream(self) -> S {
        self.0
    }
}

struct MetadataReaderFromFilter<S: MetadataReader, F: Filter> {
    inner: S,
    filter: F,
}

impl<S: MetadataReader, F: Filter> MetadataReader for MetadataReaderFromFilter<S, F> {
    type Next = IndividualsReaderFromFilter<S::Next, F>;

    fn read_metadata(
        mut self,
        map: &mut Vec<(String, String)>,
    ) -> Result<Self::Next, &'static str> {
        let mut inner_map = Vec::new();
        let next = self.inner.read_metadata(&mut inner_map)?;
        map.extend(self.filter.filter_metadata(inner_map));

        Ok(IndividualsReaderFromFilter {
            inner: next,
            filter: self.filter,
        })
    }
}

struct IndividualsReaderFromFilter<S: IndividualsReader, F: Filter> {
    inner: S,
    filter: F,
}

impl<S: IndividualsReader, F: Filter> IndividualsReader for IndividualsReaderFromFilter<S, F> {
    type Next = DataStreamFromFilter<S::Next, F>;

    fn read_individuals(mut self, list: &mut Vec<String>) -> Result<Self::Next, &'static str> {
        let mut inner_list = Vec::new();
        let next = self.inner.read_individuals(&mut inner_list)?;
        list.extend(self.filter.filter_individuals(inner_list));

        Ok(DataStreamFromFilter {
            inner: next,
            filter: self.filter,
        })
    }
}

struct DataStreamFromFilter<S: DataStream, F: Filter> {
    inner: S,
    filter: F,
}

impl<S: DataStream, F: Filter> Iterator for DataStreamFromFilter<S, F> {
    type Item = VCFData;

    fn next(&mut self) -> Option<VCFData> {
        loop {
            match self.inner.next().map(|data| self.filter.filter_item(data))? {
                DataAction::Ignore => continue,
                DataAction::Stop => return None,
                DataAction::Data(data) => return Some(data),
            }
        }
    }
}

impl<S: DataStream, F: Filter> DataStream for DataStreamFromFilter<S, F> {}

pub struct IndividualsFilter {
    individuals: HashSet<String>,
    indices: Vec<usize>,
}

impl IndividualsFilter {
    pub fn new<I, Istr>(individuals: I) -> Self
    where
        I: IntoIterator<Item = Istr>,
        Istr: Into<String>,
    {
        IndividualsFilter {
            individuals: individuals.into_iter().map(|s| s.into()).collect(),
            indices: Vec::new(),
        }
    }
}

impl Filter for IndividualsFilter {
    fn filter_individuals(&mut self, individuals: Vec<String>) -> Vec<String> {
        // For each individual in the VCF file, determine whether we are
        // interested in it and keep it only if so
        let mut filtered_individuals = Vec::new();

        for (idx, individual) in individuals.iter().enumerate() {
            if self.individuals.contains(individual) {
                self.indices.push(idx);
                filtered_individuals.push(individual.to_string())
            }
        }

        filtered_individuals
    }

    fn filter_item(&mut self, item: VCFData) -> DataAction {
        DataAction::Data(item.keep(self.indices.iter()))
        // let (mut line, mut tabs) = create_info_line(&item);

        // let (data_line, more_tabs) = create_data_line(&item, line.len(), &self.indices);
        // line.push_str(&data_line);
        // tabs.extend(more_tabs);

        // let new_item = VCFData::new(line, tabs);
        // VCFDataAction::Data(new_item)
    }
}

pub struct PolymorphismFilter {
    polymorphisms: HashSet<String>,
    count: usize,
}

impl PolymorphismFilter {
    pub fn new<I, Pstr>(polymorphisms: I) -> Self
    where
        I: IntoIterator<Item = Pstr>,
        Pstr: Into<String>,
    {
        PolymorphismFilter {
            polymorphisms: polymorphisms.into_iter().map(|s| s.into()).collect(),
            count: 0,
        }
    }
}

impl Filter for PolymorphismFilter {
    fn filter_item(&mut self, item: VCFData) -> DataAction {
        if self.count == self.polymorphisms.len() {
            DataAction::Stop
        } else if !self.polymorphisms.contains(item.identifier()) {
            DataAction::Ignore
        } else {
            self.count += 1;
            DataAction::Data(item)
        }
    }
}
