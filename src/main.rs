extern crate energy_learn;
extern crate rand;

//use energy_learn::*;
//use energy_learn::mnist;
use std::thread;
use std::sync::mpsc::channel;

fn main() {
    let root_dir = "/Users/changsoonpark/".to_string();
    let train = energy_learn::mnist::read_mnist_images(&(root_dir.clone() + "data/mnist/train-images-idx3-ubyte"),
                              &(root_dir.clone() + "data/mnist/train-labels-idx1-ubyte")).unwrap();
    let test = energy_learn::mnist::read_mnist_images(&(root_dir.clone() + "data/mnist/t10k-images-idx3-ubyte"),
                              &(root_dir.clone() + "data/mnist/t10k-labels-idx1-ubyte")).unwrap();
    //match m {
    //    Ok(x) => println!("okay {:?}", x),
    //    Err(e) => println!("error {}", e)
    //}
    println!("train and test sets are read.");
    let mut rng = rand::thread_rng();
    println!("creating a cube...");
    let num_trains = 60000;
    let num_tests = 10000;
    let cube = energy_learn::lattice::Cube::new(&mut rng, num_trains + num_tests, (28, 28, 28));
    println!("creating membranes...");
    let membrane1 = energy_learn::lattice::Membrane::xy(&cube, (28, 28), (0, 0, 0));
    let membrane2 = energy_learn::lattice::Membrane::xy(&cube, (28, 28), (0, 0, 27));
    println!("membrane1: {:?}", membrane1);
    println!("membrane2: {:?}", membrane2);
    let membranes = vec![membrane1, membrane2];
    // create the target images;
    let labels = (0..10).map(|x| {
        let mut base = vec![false; 28*28];
        for i in 0..28 {
            base[x*28 + i] = true;
        }
        base
    }).collect::<Vec<Vec<bool>>>();
    let mut lattice = energy_learn::lattice::Lattice::new(cube, membranes);
    let mut index = 0;
    println!("initializing the lattice using the training set...");
    for image in train.iter().take(num_trains) {
        if index % 10000 == 0 {
            println!("index {}", index);
        }
        let bool_data = image.data.iter().map(|&x| x > 128).collect::<Vec<bool>>();
        lattice.set_membrane(0, index, &bool_data);
        lattice.set_membrane(1, index, &bool_data);
        //lattice.set_membrane(1, index, &labels[image.label as usize]);
        index += 1;
    }
    for image in test.iter().take(num_tests) {
        if index % 5000 == 0 {
            println!("index {}", index);
        }
        let bool_data = image.data.iter().map(|&x| x > 128).collect::<Vec<bool>>();
        lattice.set_membrane(0, index, &bool_data);
        //lattice.set_membrane(1, num_trains + index, &labels[image.label as usize]);
        index += 1;
    }
    println!("update spin products after setting membranes...");
    lattice.update_spin_products();
    println!("creating a shared channel...");
    let (tx, rx) = channel();
    let observe_index = num_trains;
    thread::spawn(move || {
        let mut rng = rand::thread_rng();
        lattice.run(&mut rng, 1_000_000_000_000_000, 1_000_000_000, 1000000, 0, num_trains, observe_index, 1.0, 1.0, tx);
    });
    //lattice.run(&mut rng, 10000, 1000, 10, 1, num_trains, 1.0, 1.0, tx);

    loop {
        //let data = rx.recv();
        if let Ok(data) = rx.recv() {
            println!("expected label: {}", if observe_index < num_trains {
            train[observe_index].label } else {
            test[observe_index - num_trains].label });
            println!("received: {}", data);
        }
    }
}

//                       rng: &mut T,
//                       total_steps: usize,
//                       interm_step_size: usize,
//                       num_spin_runs_per_weight: usize,
//                       free_membrane_index: usize,
//                       free_spinset_start_index: usize,
//                       weight_gaussian_width: f64,
//                       temperature: f64,
//                       sender: mpsc::Sender<Array2<bool>>) 
    //println!("membranes set");
    //loop {}
    // println!("{:?}", lattice);
    //println!("hi");
