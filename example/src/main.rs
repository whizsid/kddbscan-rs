use kddbscan::{cluster, IntoPoint, ClusterId};
use colored::{Colorize};

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
    coordinates.push(Coordinate { x: 11.0, y: 9.0 });
    coordinates.push(Coordinate { x: 10.0, y: 8.0 });
    coordinates.push(Coordinate { x: 1.0, y: 2.0 });
    coordinates.push(Coordinate { x: 3.0, y: 1.0 });
    coordinates.push(Coordinate { x: 4.0, y: 4.0 });
    coordinates.push(Coordinate { x: 9.0, y: 0.0 });

    let clustered =  cluster(coordinates, 2, None, None);

    println!("  _______________        _______________");
    println!("  |          111|        |          111|");
    println!("  |0123456789012|        |0123456789012|");
    println!("  |-------------|        |-------------|");
    for y in 0..13 {
        
        for x in 0..13 {
            let item = clustered.iter().find(|w|{w.into_inner().x==(x as f64)&&w.into_inner().y==(y as f64)});

            let point = match item {
                Some(_)=>{
                    "X".black()
                }
                None=>{
                    " ".white()
                }
            };

            if x==0 {
                print!("{: <2}|",y);
            } 
            print!("{}", point);
        }

        print!("|{: <6}","");

        for x in 0..13 {
            let item = clustered.iter().find(|w|{w.into_inner().x==(x as f64)&&w.into_inner().y==(y as f64)});

            let point = match item {
                Some(item)=>{
                    if let ClusterId::Classified(cid)  = item.get_cluster_id() {
                        if cid==&0 {
                            "X".blue()
                        } else {
                            "X".green()
                        }
                    } else {
                        "X".red()
                    }
                }
                None=>{
                    " ".white()
                }
            };

            if x==12 {
                println!("{}|",point);
            } else {
                if x==0 {
                    print!("{: <2}|",y);
                }
                print!("{}", point);
            }

        }
    }

    println!("________________________________________");
    println!("       DATA                   RESULT    ");
}
