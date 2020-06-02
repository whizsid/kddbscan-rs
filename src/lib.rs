extern crate ordered_float;

use ordered_float::OrderedFloat;
/// Implement this trait to every point
pub trait IntoPoint: Sized {
    /// Calculating the distance to another point
    /// * `neighbor` - Other point
    fn get_distance(&self, neighbor: &Self) -> f64;

    /// Returning the unique identity of the point
    fn get_id(&self) -> u32;
}

pub struct Cluster<F> {
    points: Vec<F>,
    id: u32,
}

impl<F> Cluster<F> {
    pub fn new(id: u32) -> Self {
        Cluster { id, points: vec![] }
    }

    pub fn add_point(&mut self, point: F) {
        self.points.push(point);
    }

    pub fn get_points(&self) -> &Vec<F> {
        self.points.as_ref()
    }

    pub fn get_id(&self) -> u32 {
        self.id
    }
}

/// k-Deviation Density Based Clustering Algorithm
///
/// [Read More](https://www.researchgate.net/publication/323424266_A_k_-Deviation_Density_Based_Clustering_Algorithm)
pub struct Kddbscan<F: IntoPoint> {
    points: Vec<F>,
    k: u32,
    n: u32,
    deviation_factor: u32,
    clusters: Vec<Cluster<F>>,
}

impl<F: IntoPoint> Kddbscan<F> {
    pub fn new<T: IntoPoint>(points: Vec<T>, k: u32, n: u32, deviation_factor: u32) -> Kddbscan<T> {
        Kddbscan {
            points,
            k,
            n,
            deviation_factor,
            clusters: vec![],
        }
    }

    /// Clustering all points
    pub fn cluster(&self) {
        let mut points_iter = self.points.iter();

        let mut c = 0;
        while let Some(point) = points_iter.next() {
            c += 0;
        }
    }

    fn calculate_deviation_factor(&self, point: &F) -> Result<f64, &'static str> {
        let neighbors = self.get_mutual_neighbors(point);

        let distances: Vec<OrderedFloat<f64>> = neighbors
            .iter()
            .map(|v| OrderedFloat::from(point.get_distance(v)))
            .collect();

        let max = distances.iter().max();
        match max {
            Some(max) => {
                let max = max.into_inner();
                let filtered_distances: Vec<f64> = distances
                    .iter()
                    .filter(|distance| distance.into_inner() != max)
                    .map(|distance|{distance.into_inner()})
                    .collect();

                let avg:f64 = filtered_distances.iter().sum::<f64>() / filtered_distances.len() as f64;

                Ok(max/ avg)
            }
            None => Err("Can not find a neighbor"),
        }
    }

    fn get_mutual_neighbors<'a, 'b: 'a>(&'a self, point: &'b F) -> Vec<&F> {
        let mut neighbors = vec![];

        fn fill_mutual_neighbor<'a, 'b: 'a, T: IntoPoint>(
            inner_neighbors: &mut Vec<&'a T>,
            kddbscan: &'b Kddbscan<T>,
            point: &'b T,
            n: u32,
        ) {
            let mut points_iter = kddbscan.points.iter();
            if n >= kddbscan.n {
                inner_neighbors.push(point);
                while let Some(in_point) = points_iter.next() {
                    if in_point.get_id() != point.get_id() {
                        fill_mutual_neighbor(inner_neighbors, kddbscan, in_point, n + 1);
                    }
                }
            }
        }

        fill_mutual_neighbor(&mut neighbors, self, point, 0);

        neighbors
    }

    fn expand_cluster(&mut self, point: F, n: u32) {
        let exist_cluster_result = self.clusters.get_mut(n as usize);

        match exist_cluster_result {
            Some(exist_cluster) => {
                exist_cluster.add_point(point);
            }
            None => {
                let mut cluster = Cluster::new(n);
                cluster.add_point(point);
                self.clusters.push(cluster);
            }
        }
    }

    /// Returning the clustered points
    pub fn get_clusters(self) -> Vec<Cluster<F>> {
        self.clusters
    }
}

/// Clustering a vec of structs
/// * `points` - Vector of points. All points should implement the `IntoPoint` trait
/// * `k` - Value of k constant
/// * `n` - The minimal number of the mutual neighborhood. Default value is 1
/// * `deviation_factor` - Deviation Factor. Default value is 999999
pub fn cluster<T: IntoPoint>(
    points: Vec<T>,
    k: u32,
    n: Option<u32>,
    deviation_factor: Option<u32>,
) -> Vec<Cluster<T>> {
    let kddbscan = Kddbscan::<T>::new::<T>(
        points,
        k,
        n.unwrap_or(1),
        deviation_factor.unwrap_or(999999),
    );

    kddbscan.cluster();

    let clusters = kddbscan.get_clusters();

    clusters
}
