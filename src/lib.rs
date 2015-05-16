use std::io::prelude::*;
use std::fs::File;
use std::io;

pub fn read_mnist_images(filename: &str) -> io::Result<&'static [u8]> {
    let mut contents = try!(File::open(filename));
    let result = 
    //return io::Error::new(io::ErrorKind::Other, "nyi");
    Ok(f)
}

