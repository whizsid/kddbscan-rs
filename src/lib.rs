//! A k -Deviation Density Based Clustering Algorithm (kDDBSCAN)
//!
//! [Research Paper](https://www.researchgate.net/publication/323424266_A_k_-Deviation_Density_Based_Clustering_Algorithm)
//!
//! > Due to the adoption of global parameters, DBSCAN fails to
//! > identify clusters with different and varied densities. To
//! > solve the problem, this paper extends DBSCAN by exploiting
//! > a new density definition and proposes a novel algorithm
//! > called k -deviation density based DBSCAN (kDDBSCAN). Various
//! > datasets containing clusters with arbitrary shapes and
//! > different or varied densities are used to demonstrate the
//! > performance and investigate the feasibility and practicality
//! > of kDDBSCAN. The results show that kDDBSCAN performs
//! > better than DBSCAN.
//!
//! # Installation
//!
//! Add the following to your `Cargo.toml` file:
//!
//! ```toml
//! [dependencies]
//! kddbscan = "0.1.0"
//! ```
//!
//! # Usage
//!
//! ```rust
//! use kddbscan::{cluster, IntoPoint};
//!
//! pub struct Coordinate {
//!     pub x: f64,
//!     pub y: f64,
//! }
//!
//! // Implement IntoPoint trait to your data structur
//! impl IntoPoint for Coordinate {
//!     fn get_distance(&self, neighbor: &Coordinate) -> f64 {
//!         ((self.x - neighbor.x).powi(2) + (self.y - neighbor.y).powi(2)).powf(0.5)
//!     }
//! }
//!
//! fn main() {
//!     // Create a vector with your data
//!     let mut coordinates: Vec<Coordinate> = vec![];    
//!     coordinates.push(Coordinate { x: 11.0, y: 12.0 });
//!     coordinates.push(Coordinate { x: 0.0, y: 0.0 });
//!     coordinates.push(Coordinate { x: 12.0, y: 11.0 });
//!     coordinates.push(Coordinate { x: 11.0, y: 11.0 });
//!     coordinates.push(Coordinate { x: 1.0, y: 2.0 });
//!     coordinates.push(Coordinate { x: 3.0, y: 1.0 });
//!
//!     // Call cluster function
//!     let clustered =  cluster(coordinates, 2, None, None);
//!     let first_cluster_id = clustered.get(0).unwrap().get_cluster_id();
//!     let second_cluster_id = clustered.get(1).unwrap().get_cluster_id();
//!     
//!     assert_eq!(first_cluster_id, clustered.get(2).unwrap().get_cluster_id());
//!     assert_eq!(first_cluster_id, clustered.get(3).unwrap().get_cluster_id());
//!     assert_eq!(second_cluster_id, clustered.get(4).unwrap().get_cluster_id());
//!     assert_eq!(second_cluster_id, clustered.get(5).unwrap().get_cluster_id());
//! }
//! ```

extern crate ordered_float;

use ordered_float::OrderedFloat;
use std::collections::HashMap;

/// You should implement `IntoPoint` for all points
pub trait IntoPoint: Sized {
    /// Returns the distance to another point in the dataset
    fn get_distance(&self, neighbor: &Self) -> f64;
}

/// Cluster id types
///
/// > See `ExpandCluster` procedure in [this](https://www.researchgate.net/publication/323424266_A_k_-Deviation_Density_Based_Clustering_Algorithm) research.
#[derive(Debug, PartialEq)]
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
    fn new(index: usize, point: F) -> PointWrapper<F> {
        PointWrapper {
            point,
            index,
            cluster_id: ClusterId::Unclassified,
        }
    }

    fn set_cluster_id(&mut self, cluster_id: ClusterId) {
        self.cluster_id = cluster_id;
    }

    /// Returns the cluster id
    pub fn get_cluster_id(&self) -> &ClusterId {
        &self.cluster_id
    }

    fn get_distance(&self, wrapper: &PointWrapper<F>) -> f64 {
        self.point.get_distance(&wrapper.point)
    }

    /// Returns the index of the point
    pub fn get_id(&self) -> usize {
        self.index
    }

    /// Returns a reference to the original point
    pub fn into_inner(&self) -> &F {
        &self.point
    }
}

/// k-Deviation Density Based Clustering Algorithm
///
/// [Read More](https://www.researchgate.net/publication/323424266_A_k_-Deviation_Density_Based_Clustering_Algorithm)
struct Kddbscan<F: IntoPoint> {
    points: Vec<PointWrapper<F>>,
    k: u32,
    n: u32,
    deviation_factor: u32,
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
        }
    }

    /// Clustering all points
    pub fn cluster(&mut self) {
        let mut points_iter = self.points.iter();

        let mut c = 0;
        let mut cluster_assigns: HashMap<usize, ClusterId> = HashMap::new();
        while let Some(point) = points_iter.next() {
            let cluster_id = cluster_assigns
                .get(&point.get_id())
                .unwrap_or(&ClusterId::Unclassified);
            match cluster_id {
                ClusterId::Unclassified => {
                    let density = self.deviation_density(point).unwrap();

                    if density <= self.deviation_factor as f64 {
                        let tmp_cluster_assigns = self.expand_cluster(&cluster_assigns, point, c);
                        for (k, v) in tmp_cluster_assigns {
                            cluster_assigns.insert(k, v);
                        }
                        c += 1;
                    } else {
                        cluster_assigns.insert(point.get_id(), ClusterId::Outline);
                    }
                }
                _ => {}
            }
        }

        // Applying all stored allocations to the self points
        for (k, v) in cluster_assigns {
            self.points.get_mut(k).unwrap().set_cluster_id(v);
        }
    }
    /// Calculating the deviation density for a point
    ///
    /// Deviation factor formula:-
    /// > ![Deviation Factor Formula](https://render.githubusercontent.com/render/math?math=Dev_%7Bk%7D%28x%2Cn%29%3D%5Cfrac%7Bmax_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28x%2Cx_%7Bi%7D%29%29%7D%7Bavg_%7Bi%20%5Cvarepsilon%20argmax%28d%28x%2Cx_%7Bi%7D%29%29%20%7D%28d%28x%2Cx_%7Bi%7D%29%29%7D)
    fn deviation_density(&self, point: &PointWrapper<F>) -> Result<f64, &'static str> {
        let neighbors = self.get_mutual_neighbors(point);

        let distances: Vec<OrderedFloat<f64>> = neighbors
            .iter()
            .map(|v| OrderedFloat::from(point.get_distance(v)))
            .collect();

        let max = distances.iter().max();
        match max {
            Some(max) => {
                let max = max.into_inner();

                let all_distances: Vec<OrderedFloat<f64>> = self
                    .points
                    .iter()
                    .map(|p| OrderedFloat::from(p.get_distance(point)))
                    .collect();
                let all_max = all_distances.iter().max().unwrap().into_inner();
                let without_max: Vec<f64> = all_distances
                    .iter()
                    .map(|d| d.into_inner())
                    .filter(|d| d != &all_max)
                    .collect();

                let all_avg = without_max.iter().sum::<f64>() / without_max.len() as f64;

                Ok(max / all_avg)
            }
            None => Err("Can not find a neighbor"),
        }
    }

    /// Getting the mutual k nearest neighbors for a given point
    ///
    /// ### How to get k nearest neighbors?
    ///
    /// > 1. Determine parameter K = number of nearest neighbors
    /// > 2. Calculate the distance between the query-instance and all the training samples
    /// > 3. Sort the distance and determine nearest neighbors based on the K-th minimum distance
    /// > 4. Gather the category of the nearest neighbors
    /// > 5. Use simple majority of the category of nearest neighbors as the prediction value of the query instance
    ///
    /// [Read More](https://people.revoledu.com/kardi/tutorial/KNN/KNN_Numerical-example.html)
    ///
    /// ### How to get mutual k nearest neighbor
    ///
    /// > Given two points 𝑥𝑖 and 𝑥𝑗, if 𝑥𝑖 ∈ 𝑁𝑘(𝑥𝑗) and 𝑥𝑗 ∈ 𝑁𝑘(𝑥𝑖), then 𝑥𝑖
    /// > and 𝑥𝑗 are mutual 𝑘-nearest neighbors. The mutual 𝑘-nearest
    /// > neighborhood (mKNN) of a point 𝑥 is denoted by 𝑀𝑘(𝑥)
    ///
    /// See the [Definition 2 of this research paper](https://www.researchgate.net/publication/323424266_A_k_-Deviation_Density_Based_Clustering_Algorithm)
    fn get_mutual_neighbors<'a>(&'a self, point: &'a PointWrapper<F>) -> Vec<&'a PointWrapper<F>> {
        let mut neighbors = vec![];

        fn fill_mutual_neighbor<'a, 'b: 'a, T: IntoPoint>(
            inner_neighbors: &mut Vec<&'a PointWrapper<T>>,
            kddbscan: &'a Kddbscan<T>,
            main_point: &'a PointWrapper<T>,
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
            let mut i = 1;
            while i <= kddbscan.k {
                let in_point_result = points.get(i as usize);
                if let Some(in_point) = in_point_result {
                    // Do not stay in same point
                    if in_point.get_id() != point.get_id() {
                        // Do not go again to the start point
                        if main_point.get_id() != in_point.get_id() {
                            // Do not go again on same way
                            if !inner_neighbors
                                .iter()
                                .any(|in_neighbor| in_neighbor.get_id() == in_point.get_id())
                            {
                                // Recursively selecting mutual neighbors until original point met.
                                fill_mutual_neighbor(
                                    inner_neighbors,
                                    kddbscan,
                                    main_point,
                                    in_point,
                                    n + 1,
                                );
                            }
                        }
                    }
                }
                i += 1;
            }
        }

        fill_mutual_neighbor(&mut neighbors, self, point, point, 0);
        neighbors
    }

    /// Expanding the cluster
    fn expand_cluster(
        &self,
        core_cluster_assigns: &HashMap<usize, ClusterId>,
        point: &PointWrapper<F>,
        cluster_id: usize,
    ) -> HashMap<usize, ClusterId> {
        // Storing cluster assign details in separate variable
        // Because rust don't allowing to mutate the vector inside the loop
        let mut cluster_assigns: HashMap<usize, ClusterId> = HashMap::new();
        cluster_assigns.insert(point.get_id(), ClusterId::Classified(cluster_id));

        // We are temporary storing points inside this vector
        // See the ExpandCluster procedure in the research publication
        let mut core_points: Vec<&PointWrapper<F>> = vec![];
        core_points.push(point);

        // Creating a runtime borrow checker for the vector.
        // We want to push items into this vector while iterating over the same vector

        let mut i = 0;
        while i < core_points.len() {
            let point_i = core_points[i];

            // Getting all mutual neighbors
            let neighbors = self.get_mutual_neighbors(point_i);
            cluster_assigns.insert(point_i.get_id(), ClusterId::Classified(cluster_id));

            for point_j in neighbors {
                let point_j_cluster_id = cluster_assigns.get(&point_j.get_id()).unwrap_or(
                    core_cluster_assigns
                        .get(&point_j.get_id())
                        .unwrap_or(&ClusterId::Unclassified),
                );

                match point_j_cluster_id {
                    ClusterId::Classified(_) => {}
                    _ => {
                        cluster_assigns.insert(point_j.get_id(), ClusterId::Classified(cluster_id));
                    }
                }

                let dev_density = {
                    self.deviation_density(point_j)
                        .expect("Can not calculate deviation factor.")
                };

                let density_reachable = self.density_reachable(point_j, point_i);

                if dev_density <= self.deviation_factor as f64 && density_reachable {
                    if !core_points.iter().any(|p| p.get_id() == point_j.get_id()) {
                        core_points.push(point_j);
                    }
                } else {
                    cluster_assigns.insert(point_j.get_id(), ClusterId::Outline);
                }
            }

            i += 1;
        }

        cluster_assigns
    }

    /// Checking the that weather two points are density reachable
    /// * `p` - First point
    /// * `q` - Second point
    ///
    /// > A point 𝑝 is density reachable from a point 𝑞 if the
    /// > following conditions are satisfied:
    ///
    /// > (1) ![Density Reachable 1st Equation](https://render.githubusercontent.com/render/math?math=max%28max_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28p%2Cx_%7Bi%7D%29%29%2Fmax_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28q%2Cx_%7Bi%7D%29%29%20%2C%20max_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28q%2Cx_%7Bi%7D%29%29%2Fmax_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28p%2Cx_%7Bi%7D%29%29%20%29%20%5Cleq%20%5Calpha%20)
    /// > (2) ![Density Reachable 2nd Equation](https://render.githubusercontent.com/render/math?math=max%28avg_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28p%2Cx_%7Bi%7D%29%29%2Favg_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28q%2Cx%29%29%20%2C%20avg_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28q%2Cx_%7Bi%7D%29%29%2Favg_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28p%2Cx%29%29%20%29%20%5Cleq%20%5Calpha%20)
    /// > (3) ![Density Reachable 3rd Equation](https://render.githubusercontent.com/render/math?math=max%28max_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28p%2Cx_%7Bi%7D%29%29%2Favg_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28q%2Cx_%7Bi%7D%29%29%20%2C%20max_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28q%2Cx_%7Bi%7D%29%29%2Favg_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28p%2Cx_%7Bi%7D%29%29%20%29%20%5Cleq%20%5Calpha%20)
    /// > (4) ![Density Reachable 4th Equation](https://render.githubusercontent.com/render/math?math=max%28d%28p%2Cq%29%2Fmax_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28p%2Cx_%7Bi%7D%29%29%20%2C%20max_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28p%2Cx_%7Bi%7D%29%29%2Fd%28p%2Cq%29%20%29%20%5Cleq%20%5Calpha%20)
    /// > (5) ![Density Reachable 5th Equation](https://render.githubusercontent.com/render/math?math=max%28d%28p%2Cq%29%2Fmax_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28q%2Cx_%7Bi%7D%29%29%20%2C%20max_%7Bx_%7Bi%7D%5Cepsilon%20NM_%7Bk%7D%28x%2Cn%29%7D%28d%28q%2Cx_%7Bi%7D%29%29%2Fd%28p%2Cq%29%20%29%20%5Cleq%20%5Calpha%20)
    ///
    /// > Where `d()` is the density function and ![NMk(x,n)](https://render.githubusercontent.com/render/math?math=NM_%7Bk%7D%28x%2Cn%29) is the mutual k nearest neighbors of x point
    fn density_reachable(&self, p: &PointWrapper<F>, q: &PointWrapper<F>) -> bool {
        let p_mutual_neighbors = self.get_mutual_neighbors(p);
        let q_mutual_neighbors = self.get_mutual_neighbors(q);

        let p_distances: Vec<OrderedFloat<f64>> = p_mutual_neighbors
            .iter()
            .map(|point| OrderedFloat::from(point.get_distance(p)))
            .collect();
        let q_distances: Vec<OrderedFloat<f64>> = q_mutual_neighbors
            .iter()
            .map(|point| OrderedFloat::from(point.get_distance(q)))
            .collect();

        let p_max = p_distances.iter().max().unwrap().into_inner();
        let q_max = q_distances.iter().max().unwrap().into_inner();
        let p_len = p_distances.len();
        let q_len = q_distances.len();
        let p_sum = p_distances
            .iter()
            .map(|ord_float| ord_float.into_inner())
            .sum::<f64>();
        let q_sum = q_distances
            .iter()
            .map(|ord_float| ord_float.into_inner())
            .sum::<f64>();
        let p_avg = p_sum / (p_len as f64);
        let q_avg = q_sum / (q_len as f64);
        let p_to_q = p.get_distance(q);

        let first = (p_max / q_max).max(q_max / p_max) <= self.deviation_factor as f64;
        let second = (p_avg / q_avg).max(q_avg / p_avg) <= self.deviation_factor as f64;
        let third = (p_max / q_avg).max(q_max / p_avg) <= self.deviation_factor as f64;
        let fourth = (p_to_q / p_max).max(p_max / p_to_q) <= self.deviation_factor as f64;
        let fifth = (p_to_q / q_max).max(q_max / p_to_q) <= self.deviation_factor as f64;

        first && second && third && fourth && fifth
    }

    /// Returning the clustered points
    pub fn get_clustered(self) -> Vec<PointWrapper<F>> {
        self.points
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
) -> Vec<PointWrapper<T>> {
    let mut kddbscan = Kddbscan::<T>::new::<T>(
        points,
        k,
        n.unwrap_or(1),
        deviation_factor.unwrap_or(999999),
    );

    kddbscan.cluster();

    let clusters = kddbscan.get_clustered();

    clusters
}

#[cfg(test)]
mod tests {
    use super::*;

    pub struct Coordinate {
        pub x: f64,
        pub y: f64,
    }

    impl IntoPoint for Coordinate {
        fn get_distance(&self, neighbor: &Coordinate) -> f64 {
            ((self.x - neighbor.x).powi(2) + (self.y - neighbor.y).powi(2)).powf(0.5)
        }
    }

    fn create_kddbscan(k: u32) -> Kddbscan<Coordinate> {
        let mut coordinates: Vec<Coordinate> = vec![];
        coordinates.push(Coordinate { x: 11.0, y: 12.0 });
        coordinates.push(Coordinate { x: 0.0, y: 0.0 });
        coordinates.push(Coordinate { x: 12.0, y: 11.0 });
        coordinates.push(Coordinate { x: 11.0, y: 11.0 });
        coordinates.push(Coordinate { x: 1.0, y: 2.0 });
        coordinates.push(Coordinate { x: 3.0, y: 1.0 });

        Kddbscan::<Coordinate>::new::<Coordinate>(coordinates, k, 1, 999999)
    }

    #[test]
    fn test_mutual_neighbor() {
        let kddbscan_2 = create_kddbscan(2);
        let point_wrapper = kddbscan_2.points.get(0).unwrap();
        let mutual_neighbors = kddbscan_2.get_mutual_neighbors(point_wrapper);
        assert_eq!(mutual_neighbors.len(), 2);
        assert_eq!(mutual_neighbors.get(0).unwrap().get_id(), 3);
        assert_eq!(mutual_neighbors.get(1).unwrap().get_id(), 2);

        let kddbscan_3 = create_kddbscan(3);
        let point_wrapper = kddbscan_3.points.get(0).unwrap();
        let mutual_neighbors = kddbscan_3.get_mutual_neighbors(point_wrapper);
        assert_eq!(mutual_neighbors.len(), 5);
    }

    #[test]
    fn test_deviation_density() {
        let kddbscan = create_kddbscan(2);
        let point_wrapper = kddbscan.points.get(0).unwrap();
        let deviation_density = kddbscan.deviation_density(point_wrapper).unwrap();
        assert!(deviation_density<0.24 && deviation_density>0.22);


        let point_wrapper = kddbscan.points.get(1).unwrap();
        let deviation_density = kddbscan.deviation_density(point_wrapper).unwrap();

        assert!(deviation_density<0.61 && deviation_density>0.60);
    }

    #[test]
    fn test_density_reachable(){
        let kddbscan = create_kddbscan(2);
        let density_reachable = kddbscan.density_reachable(
            kddbscan.points.get(0).unwrap()
            , 
            kddbscan.points.get(1).unwrap()
        );

        assert!(density_reachable);
    }
}
