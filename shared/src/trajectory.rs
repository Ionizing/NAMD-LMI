use super::Structure;

pub struct Trajectory {
    pub dat: Vec<Structure>,
}


impl From<Vec<Structure>> for Trajectory {
    fn from(dat: Vec<Structure>) -> Self {
        Self { dat }
    }
}
