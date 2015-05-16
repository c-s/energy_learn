extern crate energy_learn;

use energy_learn::*;

fn main() {
    let m = read_mnist_images("/Users/chang-soon/data/mnist/train-images-idx3-ubyte",
                              "/Users/chang-soon/data/mnist/train-labels-idx1-ubyte");
    match m {
        Ok(x) => println!("okay {:?}", x),
        Err(e) => println!("error {}", e)
    }
    println!("hi");
}
