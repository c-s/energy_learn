extern crate energy_learn;

use energy_learn::*;

fn main() {
    let m = read_mnist_images("/Users/changsoonpark/data/mnist/train-labels-idx1-ubyte"); 
    match m {
        Ok(x) => println!("okay {}", x),
        Err(e) => println!("error {}", e)
    }
    println!("hi");
}
