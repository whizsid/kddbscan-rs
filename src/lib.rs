/// Implement this trait to every point
pub trait IntoPoint: Sized {
    fn get_distance(&self, neighbor: Self) -> f64;
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
pub struct Kddbscan<F> {
    points: Vec<F>,
    k: u32,
    n: u32,
    deviation_factor: u32,
    clusters: Vec<Cluster<F>>,
}

impl<F> Kddbscan<F> {
    pub fn new<T: IntoPoint>(points: Vec<T>, k: u32, n: u32, deviation_factor: u32) -> Kddbscan<T> {
        Kddbscan {
            points,
            k,
            n,
            deviation_factor,
            clusters: vec![],
        }
    }

    /// Clustergin all points
    pub fn cluster(&self) {}

    fn calculate_deviation_factor<T: IntoPoint>(point: T) -> f32 {
        0.0
    }

    fn expand_cluster<T: IntoPoint>(&mut self, point: T, n: u32) {}

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
