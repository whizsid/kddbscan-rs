use kddbscan::{cluster, IntoPoint};

pub struct Coordinate {
    pub id: u32,
    pub x: f64,
    pub y: f64
}

impl IntoPoint for Coordinate {
    fn get_distance(&self, neighbor: &Coordinate) -> f64 {
        0.0
    }
}

fn main() {
    println!("Hello, world!");

    let mut coordinates: Vec<Coordinate> = vec!();
    coordinates.push(Coordinate{ id: 0, x: 12.0, y:11.0});
    coordinates.push(Coordinate{ id: 1, x: 0.0, y:0.0});
    coordinates.push(Coordinate{ id: 2, x: 11.0, y:11.0});
    coordinates.push(Coordinate{ id: 3, x: 1.0, y:2.0});

    cluster(coordinates, 5, None, None);
}
