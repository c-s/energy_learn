use std::fmt;
use std::io::prelude::*;
use std::fs::File;
use std::io;

pub fn read_mnist_images(imagefile: &str, labelfile: &str) -> io::Result<u8> { //-> io::Result<&'static [u8]> {
    let mut image_buf = io::BufReader::new(try!(File::open(imagefile)));
    let mut label_buf = io::BufReader::new(try!(File::open(labelfile)));

    // first read the magic words.
    let mut u32buf0: [u8; 4] = [0u8; 4];
    let mut u32buf1: [u8; 4] = [0u8; 4];
    image_buf.read(&mut u32buf0[..]);
    label_buf.read(&mut u32buf1[..]);
    if u32buf0 != [0u8, 0, 8, 1] || u32buf1 != [0u8, 0, 8, 3] {
        return Err(io::Error::new(io::ErrorKind::Other, "read_mnist_images: magic number does not match."))
    }
    file.read(&mut u32buf0);
    // i don't know how to get int out of 4 bytes using an rust library.
    let num_image_items = ((u32buf0[0] as u32) << 24) +
                    ((u32buf0[1] as u32) << 16) +
                    ((u32buf0[2] as u32) << 8) +
                    (u32buf0[3] as u32);
    file.read(&mut u32buf0);
    // i don't know how to get int out of 4 bytes using an rust library.
    let num_label_items = ((u32buf1[0] as u32) << 24) +
                    ((u32buf1[1] as u32) << 16) +
                    ((u32buf1[2] as u32) << 8) +
                    (u32buf1[3] as u32);
    if num_image_items != num_label_ites {
        return Err(io::Error::new(io::ErrorKind::Other, "the number of the images do not match with those of the labels."));
    }
    println!("num_items {}", num_items);
    
    Ok(12u8)
/*
        for x in &magic_buf {
            println!("x: {}", x);
        }
        println!("vec: {}", magic_buf[3]);
        Ok(0u8)
    }*/
    /*
    let mut bytestream = file.bytes();
    let first_elem: u8 = match bytestream.next() {
        Some(m) => try!(m),
        None => panic!("no")
    };
    println!("hello");
    println!("{}", first_elem);
    //return io::Error::new(io::ErrorKind::Other, "nyi");
    //Ok(f)
    Ok(first_elem)*/
}

//impl fmt::Display for Vec<u8> {
//    fn fmt(&self, &mut fmt: fmt::Formatter) -> fmt::Result {
//        Ok(())
//    }
//}

struct MnistImage {
    data: [u8; 28 * 28],
    label: u8
}
