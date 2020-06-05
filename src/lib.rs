extern crate ordered_float;

use ordered_float::OrderedFloat;
use std::cell::RefCell;

/// Implement this trait to every point
pub trait IntoPoint: Sized {
    /// Calculating the distance to another point
    /// * `neighbor` - Other point
    fn get_distance(&self, neighbor: &Self) -> f64;
}

/// Cluster id types
/// See `ExpandCluster` procedure in [this](https://www.researchgate.net/publication/323424266_A_k_-Deviation_Density_Based_Clustering_Algorithm) research.
pub enum ClusterId {
    Outline,
    Unclassified,
    Classified(usize),
}

/// This struct is using to store temporary values for points
/// such as cluster_id and index of the point
pub struct PointWrapper<F: IntoPoint> {
    point: F,
    cluster_id: ClusterId,
    index: usize,
}

impl<F: IntoPoint> PointWrapper<F> {
    pub fn new(index: usize, point: F) -> PointWrapper<F> {
        PointWrapper {
            point,
            index,
            cluster_id: ClusterId::Unclassified,
        }
    }

    pub fn set_cluster_id(&mut self, cluster_id: ClusterId) {
        self.cluster_id = cluster_id;
    }

    pub fn get_cluster_id(&self) -> &ClusterId {
        &self.cluster_id
    }

    pub fn get_distance(&self, wrapper: &PointWrapper<F>) -> f64 {
        self.point.get_distance(&wrapper.point)
    }

    pub fn get_id(&self) -> usize {
        self.index
    }
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
    points: Vec<PointWrapper<F>>,
    k: u32,
    n: u32,
    deviation_factor: u32,
    clusters: Vec<Cluster<F>>,
}

impl<F: IntoPoint> Kddbscan<F> {
    pub fn new<T: IntoPoint>(points: Vec<T>, k: u32, n: u32, deviation_factor: u32) -> Kddbscan<T> {
        let mut i = 0;
        let mut wrappers: Vec<PointWrapper<T>> = vec![];

        for point in points {
            wrappers.push(PointWrapper::new(i, point));
            i += 1;
        }

        Kddbscan {
            points: wrappers,
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

    fn calculate_deviation_factor(&self, point: &PointWrapper<F>) -> Result<f64, &'static str> {
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
                    .map(|distance| distance.into_inner())
                    .collect();

                let avg: f64 =
                    filtered_distances.iter().sum::<f64>() / filtered_distances.len() as f64;

                Ok(max / avg)
            }
            None => Err("Can not find a neighbor"),
        }
    }

    fn get_mutual_neighbors<'a, 'b: 'a>(
        &'a self,
        point: &'b PointWrapper<F>,
    ) -> Vec<&PointWrapper<F>> {
        let mut neighbors = vec![];

        fn fill_mutual_neighbor<'a, 'b: 'a, T: IntoPoint>(
            inner_neighbors: &mut Vec<&'a PointWrapper<T>>,
            kddbscan: &'b Kddbscan<T>,
            point: &'b PointWrapper<T>,
            n: u32,
        ) {
            // Checking the minimum number of mutual neihborhood
            if n >= kddbscan.n {
                if n != 0 {
                    inner_neighbors.push(point);
                }
            }

            let mut points: Vec<&PointWrapper<T>> =
                kddbscan.points.iter().map(|point| point).collect();

            // Sorting other points by distance to the selected point
            points.sort_by(|a, b| {
                OrderedFloat::from(a.get_distance(point))
                    .cmp(&OrderedFloat::from(b.get_distance(point)))
            });

            // Selecting only k number of values
            let mut i = 0;
            while i < kddbscan.k {
                let in_point_result = points.get(i as usize);
                if let Some(in_point) = in_point_result {
                    if in_point.get_id() != point.get_id() {
                        // Recursively selecting mutual neighbors until original point met.
                        fill_mutual_neighbor(inner_neighbors, kddbscan, in_point, n + 1);
                    }
                }
                i += 1;
            }
        }

        fill_mutual_neighbor(&mut neighbors, self, point, 0);

        // Removing duplicates
        neighbors.sort_by(|a, b| a.get_id().cmp(&b.get_id()));
        neighbors.dedup_by(|a, b| a.get_id() == b.get_id());

        neighbors
    }

    fn expand_cluster(&mut self, point: &mut PointWrapper<F>, cluster_id: usize) {
        let ref_cell = RefCell::from(self);
        let mut core_points: Vec<&mut PointWrapper<F>> = vec![];
        core_points.push(point);

        for cur_point in core_points.iter_mut() {
            let self_borrow = ref_cell.borrow();
            let neighbor_ids: Vec<usize> = self_borrow
                .get_mutual_neighbors(cur_point)
                .iter()
                .map(|wrapper| wrapper.get_id())
                .collect();

            cur_point.set_cluster_id(ClusterId::Classified(cluster_id));

            let mut self_borrow_mut = ref_cell.borrow_mut();
            let neighbors_iter_mut = self_borrow_mut.points.iter_mut();
            let neighbors: Vec<&mut PointWrapper<F>> = neighbors_iter_mut
                .filter(|wrapper| neighbor_ids.iter().any(|id| id == &wrapper.get_id()))
                .collect();

            for neighbor in neighbors {
                match neighbor.get_cluster_id() {
                    ClusterId::Classified(_) => {}
                    _ => {
                        neighbor.set_cluster_id(ClusterId::Classified(cluster_id));
                    }
                }

                let dev_density = self_borrow
                    .calculate_deviation_factor(neighbor)
                    .expect("Can not calculate deviation factor.")
                    as u32;

                if dev_density <= self_borrow.deviation_factor
                    && self_borrow.density_reachable(neighbor, cur_point)
                {
                    core_points.push(neighbor);
                } else {
                    neighbor.set_cluster_id(ClusterId::Outline);
                }
            }
        }
    }

    /// Checking the weather two points are density reachable
    /// * `p` - First point
    /// * `q` - Second point
    fn density_reachable(&self, p: &PointWrapper<F>, q: &PointWrapper<F>) -> bool {
        true
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
