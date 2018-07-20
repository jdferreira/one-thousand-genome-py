use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::io::{BufRead, Result, Error, ErrorKind};

pub struct Population {
    individual_to_group: HashMap<String, String>,
    group_to_individuals: HashMap<String, Vec<String>>,
}

impl Population {
    pub fn new() -> Self {
        Population {
            individual_to_group: HashMap::new(),
            group_to_individuals: HashMap::new(),
        }
    }

    pub fn add_individual(&mut self, individual: String, group: String) {
        self.individual_to_group.insert(individual, group);
    }

    pub fn group<Istr>(&self, individual: Istr) -> &str
    where
        Istr: AsRef<str>,
    {
        self.individual_to_group
            .get(individual.as_ref())
            .map(|s| &**s)
            .unwrap_or("???")
    }

    pub fn remove_individual<Istr>(&mut self, individual: Istr)
    where
        Istr: Into<String>,
    {
        if let Entry::Occupied(entry) = self.individual_to_group.entry(individual.into()) {
            let group = entry.remove();
            self.group_to_individuals.remove(&group);
        }
    }

    pub fn individuals(&mut self) -> impl Iterator<Item = &str> {
        self.individual_to_group.keys().map(|s| &**s)
    }

    pub fn groups(&mut self) -> impl Iterator<Item = &str> {
        self.group_to_individuals.keys().map(|s| &**s)
    }

    pub fn items(&mut self) -> impl Iterator<Item = (&str, &str)> {
        self.individual_to_group.iter().map(|(s1, s2)| (&**s1, &**s2))
    }

    pub fn has_individual<Istr>(&self, individual: Istr) -> bool
    where
        Istr: AsRef<str>,
    {
        self.individual_to_group.contains_key(individual.as_ref())
    }

    pub fn labels<Istr, I>(&self, individuals: I) -> impl Iterator<Item = &str>
    where
        I: IntoIterator<Item = Istr>,
        Istr: AsRef<str>,
    {
        individuals.into_iter().map(move |ind| self.group(ind))
    }
    
    pub fn parse<R: BufRead>(reader: R) -> Result<Population> {
        let mut population = Population::new();
        
        for line in reader.lines() {
            let line = line?;
            let line = line.trim_right_matches('\n');
            
            // Remove possible comments
            let line = if let Some(comment_start) = line.find('#') {
                &line[comment_start..]
            } else {
                line
            };
            
            // Remove extra heading and trailing spaces
            let line = line.trim();
            
            // Ignore empty lines
            if line.len() == 0 {
                continue;
            }
            
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if fields.len() != 2 {
                let msg = format!("Line {:?} is invalid: needs 2 fields, found {}", line, fields.len());
                return Err(Error::new(ErrorKind::InvalidData, msg))
            }
            
            let identifier = fields[0];
            let group = fields[1];
            population.add_individual(identifier.to_string(), group.to_string());
        }
        
        Ok(population)
    }
}

