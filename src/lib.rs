extern crate ordered_float;

use ordered_float::OrderedFloat;
use std::collections::HashMap;

/// Implement this trait to every point
pub trait IntoPoint: Sized {
    /// Calculating the distance to another point
    /// * `neighbor` - Other point
    fn get_distance(&self, neighbor: &Self) -> f64;
}

/// Cluster id types
/// See `ExpandCluster` procedure in [this](https://www.researchgate.net/publication/323424266_A_k_-Deviation_Density_Based_Clustering_Algorithm) research.
#[derive(Debug)]
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
    pub fn cluster(&mut self) {
        let mut points_iter = self.points.iter();

        let mut c = 0;
        let mut cluster_assigns:HashMap<usize, ClusterId>  = HashMap::new();
        while let Some(point) = points_iter.next() {
            match point.get_cluster_id() {
                ClusterId::Unclassified=>{

                    let density = self.calculate_deviation_factor(point).unwrap();

                    if density <= self.deviation_factor as f64{
                        let tmp_cluster_assigns = self.expand_cluster(point, c);
                        for (k, v) in tmp_cluster_assigns {
                            cluster_assigns.insert(k, v);
                        }
                        c+=1;
                    } else {
                        cluster_assigns.insert(point.get_id(), ClusterId::Outline);
                    }
                }
                _=>{}
            }
        }

        // Applying all stored allocations to the self points
        for (k, v) in cluster_assigns {
            println!("{}:{:?}",k,v);
            self.points.get_mut(k).unwrap().set_cluster_id(v);
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
                
                let all_distances: Vec<OrderedFloat<f64>> = self.points.iter().map(|p|{OrderedFloat::from(p.get_distance(point))}).collect();
                let all_max = all_distances.iter().max().unwrap().into_inner();
                let without_max: Vec<f64> = all_distances.iter().map(|d|{d.into_inner()}).filter(|d|{d!=&all_max}).collect();

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
                    if in_point.get_id() != point.get_id() {
                        if main_point.get_id() != in_point.get_id() {
                            if !inner_neighbors.iter().any(|in_neighbor|{in_neighbor.get_id()==in_point.get_id()}) {
                                // Recursively selecting mutual neighbors until original point met.
                                fill_mutual_neighbor(inner_neighbors, kddbscan, main_point, in_point, n + 1);
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

    fn expand_cluster(&self, point: &PointWrapper<F>, cluster_id: usize) -> HashMap<usize, ClusterId> {

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

            for point_j in neighbors {
                match point_j.get_cluster_id() {
                    ClusterId::Classified(_) => {}
                    _ => {
                        cluster_assigns.insert(point_j.get_id(), ClusterId::Classified(cluster_id));
                    }
                }

                let dev_density = {
                    self.calculate_deviation_factor(point_j)
                        .expect("Can not calculate deviation factor.") 
                };

                let density_reachable = self.density_reachable(point_j, point_i);

                if dev_density <= self.deviation_factor as f64 && density_reachable
                {
                    core_points.push(point_j);
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

        let p_distances: Vec<OrderedFloat<f64>> = p_mutual_neighbors.iter().map(|point|{ OrderedFloat::from(point.get_distance(p))}).collect();
        let q_distances: Vec<OrderedFloat<f64>> = q_mutual_neighbors.iter().map(|point|{ OrderedFloat::from(point.get_distance(q))}).collect();

        let p_max = p_distances.iter().max().unwrap().into_inner();
        let q_max = q_distances.iter().max().unwrap().into_inner();
        let p_len = p_distances.len();
        let q_len = q_distances.len();
        let p_sum = p_distances.iter().map(|ord_float|{ord_float.into_inner()}).sum::<f64>();
        let q_sum = q_distances.iter().map(|ord_float|{ord_float.into_inner()}).sum::<f64>();
        let p_avg = p_sum/(p_len as f64);
        let q_avg = q_sum/(q_len as f64);
        let p_to_q = p.get_distance(q);

        let first = (p_max/q_max).max(q_max/p_max) <= self.deviation_factor as f64;
        let second = (p_avg/q_avg).max(q_avg/p_avg) <= self.deviation_factor as f64;
        let third = (p_max/q_avg).max(q_max/p_avg) <= self.deviation_factor as f64;
        let fourth = (p_to_q/p_max).max(p_max/p_to_q) <= self.deviation_factor as f64;
        let fifth = (p_to_q/q_max).max(q_max/p_to_q) <= self.deviation_factor as f64;

        first && second && third && fourth && fifth
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
    let mut kddbscan = Kddbscan::<T>::new::<T>(
        points,
        k,
        n.unwrap_or(1),
        deviation_factor.unwrap_or(999999),
    );

    kddbscan.cluster();

    let clusters = kddbscan.get_clusters();

    clusters
}
