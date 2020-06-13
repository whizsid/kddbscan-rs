use kddbscan::{cluster, IntoPoint};

pub struct Coordinate {
    pub x: f64,
    pub y: f64,
}

impl IntoPoint for Coordinate {
    fn get_distance(&self, neighbor: &Coordinate) -> f64 {
        ((self.x - neighbor.x).powi(2) + (self.y - neighbor.y).powi(2)).powf(0.5)
    }
}

fn main() {
    let mut coordinates: Vec<Coordinate> = vec![];
    coordinates.push(Coordinate { x: 11.0, y: 12.0 });
    coordinates.push(Coordinate { x: 0.0, y: 0.0 });
    coordinates.push(Coordinate { x: 12.0, y: 11.0 });
    coordinates.push(Coordinate { x: 11.0, y: 11.0 });
    coordinates.push(Coordinate { x: 1.0, y: 2.0 });
    coordinates.push(Coordinate { x: 3.0, y: 1.0 });

    let clustered =  cluster(coordinates, 2, None, None);

    let point_wrapper = clustered.get(0).unwrap();
    let _point = point_wrapper.into_inner();

}
