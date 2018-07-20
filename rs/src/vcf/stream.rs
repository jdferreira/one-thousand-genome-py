use std::borrow::Borrow;
use std::io::BufRead;

/// The indices of each column in the VCF data lines
enum VCFColumn {
    CHROMOSOME = 0,
    POS = 1,
    ID = 2,
    REF = 3,
    ALT = 4,
    QUAL = 5,
    FILTER = 6,
    INFO = 7,
    FORMAT = 8,
}

static FIRST_DATA_COLUMN: usize = VCFColumn::FORMAT as usize + 1;

/// A VCF data line
pub struct VCFData {
    line: String,
    ranges: Vec<(usize, usize)>,
}

impl VCFData {
    fn new(line: String, ranges: Vec<(usize, usize)>) -> VCFData {
        VCFData { line, ranges }
    }

    fn field(&self, index: usize) -> &str {
        let (start, end) = self.ranges[index];
        &self.line[start..end]
    }

    pub fn chromosome(&self) -> &str {
        self.field(VCFColumn::CHROMOSOME as usize)
    }

    pub fn position(&self) -> &str {
        self.field(VCFColumn::POS as usize)
    }

    pub fn identifier(&self) -> &str {
        self.field(VCFColumn::ID as usize)
    }

    pub fn reference(&self) -> &str {
        self.field(VCFColumn::REF as usize)
    }

    pub fn alternatives(&self) -> Vec<&str> {
        self.field(VCFColumn::ALT as usize).split(',').collect()
    }

    pub fn quality(&self) -> &str {
        self.field(VCFColumn::QUAL as usize)
    }

    pub fn filter(&self) -> &str {
        self.field(VCFColumn::FILTER as usize)
    }

    pub fn info(&self) -> &str {
        self.field(VCFColumn::INFO as usize)
    }

    pub fn format(&self) -> &str {
        self.field(VCFColumn::FORMAT as usize)
    }

    pub fn genotype(&self, index: usize) -> &str {
        self.field(FIRST_DATA_COLUMN + index)
    }

    pub fn genotypes(&self) -> impl Iterator<Item = &str> {
        let range = FIRST_DATA_COLUMN..self.ranges.len();
        range.map(move |idx| self.field(idx))
    }

    pub(crate) fn keep<I, T>(self, indices: I) -> VCFData
    where
        T: Borrow<usize>,
        I: Iterator<Item = T>,
    {
        let mut ranges: Vec<_> = self.ranges[..FIRST_DATA_COLUMN].into();
        ranges.extend(
            indices
                .into_iter()
                .map(|i| self.ranges[i.borrow() + FIRST_DATA_COLUMN]),
        );

        VCFData {
            line: self.line,
            ranges,
        }
    }

    pub(crate) fn debug(&self) {
        for (idx, (start, end)) in self.ranges.iter().enumerate() {
            let slice = self.line.get(*start..*end);
            println!("{}: {}-{}\t{:?}", idx, start, end, slice);
        }
    }
}

pub struct StreamUnfolded<S: DataStream> {
    pub metadata: Vec<(String, String)>,
    pub individuals: Vec<String>,
    pub stream: S,
}

pub trait MetadataReader {
    type Next: IndividualsReader;

    fn read_metadata(
        self,
        metadata: &mut Vec<(String, String)>,
    ) -> Result<Self::Next, &'static str>;

    fn unfold(self) -> Result<StreamUnfolded<<Self::Next as IndividualsReader>::Next>, &'static str>
    where
        Self: Sized,
    {
        let mut metadata = Vec::new();
        let mut individuals = Vec::new();

        let stream = self.read_metadata(&mut metadata)?;
        let stream = stream.read_individuals(&mut individuals)?;

        Ok(StreamUnfolded {
            metadata,
            individuals,
            stream,
        })
    }
}

pub trait IndividualsReader {
    type Next: DataStream;

    fn read_individuals(self, list: &mut Vec<String>) -> Result<Self::Next, &'static str>;
}

pub trait DataStream: Iterator<Item = VCFData> {}

pub fn from_text_stream<R: BufRead>(text_stream: R) -> impl MetadataReader {
    MetadataReaderFromTextStream::new(text_stream)
}

struct MetadataReaderFromTextStream<R: BufRead> {
    reader: R,
}

impl<R: BufRead> MetadataReaderFromTextStream<R> {
    fn new(text_stream: R) -> MetadataReaderFromTextStream<R> {
        MetadataReaderFromTextStream {
            reader: text_stream,
        }
    }
}

impl<R: BufRead> MetadataReader for MetadataReaderFromTextStream<R> {
    type Next = IndividualsReaderFromTextStream<R>;

    fn read_metadata(
        mut self,
        map: &mut Vec<(String, String)>,
    ) -> Result<Self::Next, &'static str> {
        let mut line = String::new();
        let mut count = 0;

        loop {
            line.clear();

            let line = match self.reader.read_line(&mut line) {
                Err(_) => return Err("Error reading the text stream"),
                Ok(len) if len == 0 && count == 0 => return Err("VCF file is empty"),
                Ok(len) if len == 0 => return Err("VCF file contains only metadata"),
                Ok(_) => &line,
            };

            count += 1;

            let line = line.trim_right_matches("\n");
            if line.starts_with("##") {
                let fields = line[2..].splitn(2, "=").collect::<Vec<_>>();
                map.push((fields[0].to_string(), fields[1].to_string()));
            } else {
                return Ok(IndividualsReaderFromTextStream {
                    reader: self.reader,
                    line: line.to_string(),
                });
            }
        }
    }
}

struct IndividualsReaderFromTextStream<R: BufRead> {
    reader: R,
    line: String,
}
impl<R: BufRead> IndividualsReader for IndividualsReaderFromTextStream<R> {
    type Next = DataStreamFromTextStream<R>;

    fn read_individuals(self, list: &mut Vec<String>) -> Result<Self::Next, &'static str> {
        let start_len = list.len();

        let fields = self.line.split("\t").collect::<Vec<_>>();
        let genotype = fields[VCFColumn::FORMAT as usize] == "FORMAT";
        let start_index = if genotype {
            VCFColumn::FORMAT as usize + 1
        } else {
            VCFColumn::FORMAT as usize
        };
        list.extend(fields[start_index..].iter().map(|&s| s.to_string()));

        if list.len() == start_len {
            Err("Line does not contain any individuals")
        } else {
            Ok(DataStreamFromTextStream {
                reader: self.reader,
                genotype,
            })
        }
    }
}

struct DataStreamFromTextStream<R: BufRead> {
    reader: R,
    genotype: bool,
}

impl<R: BufRead> Iterator for DataStreamFromTextStream<R> {
    type Item = VCFData;

    fn next(&mut self) -> Option<VCFData> {
        let mut line = String::new();
        if self.reader.read_line(&mut line).is_err() || line.len() == 0 {
            return None;
        }

        if line.ends_with('\n') {
            line.pop();
        }

        let mut tabs = line
            .chars()
            .enumerate()
            .filter(|(_, c)| *c == '\t')
            .map(|(idx, _)| idx)
            .collect::<Vec<_>>();

        if !self.genotype {
            let idx = VCFColumn::FORMAT as usize - 1;
            let last_tab = tabs[idx];
            tabs.insert(idx, last_tab);
        }

        let starts = vec![0]
            .into_iter()
            .chain(tabs.iter().cloned().map(|idx| idx + 1));
        let ends = tabs.iter().cloned().chain(vec![line.len()].into_iter());
        let ranges = starts.zip(ends).collect();

        Some(VCFData::new(line, ranges))
    }
}
impl<R: BufRead> DataStream for DataStreamFromTextStream<R> {}

mod test {
    use std::io;

    use super::*;

    fn main() -> Result<(), &'static str> {
        let mut metadata = Vec::new();
        let mut individuals = Vec::new();

        let stdin = io::stdin();
        let handle = stdin.lock();
        let s = from_text_stream(handle);
        let s = s.read_metadata(&mut metadata)?;
        let s = s.read_individuals(&mut individuals)?;

        for _polymorphism in s {
            break;
        }

        Ok(())
    }
}
