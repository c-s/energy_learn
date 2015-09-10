extern crate itertools;

use std::fmt;
use std::io::prelude::*;
use std::fs::File;
use std::io;
use self::itertools::Itertools;

const IMAGE_SIZE: usize = 28 * 28;
const IMAGE_NUM_ROWS: usize = 28;
const IMAGE_NUM_COLS: usize = 28;

pub fn read_mnist_images(imagefile: &str, labelfile: &str) -> io::Result<Vec<MnistImage>> { //-> io::Result<&'static [u8]> {
    let mut image_buf = io::BufReader::new(try!(File::open(imagefile)));
    let mut label_buf = io::BufReader::new(try!(File::open(labelfile)));

    // first read the magic words.
    let mut u32buf0: [u8; 4] = [0u8; 4];
    let mut u32buf1: [u8; 4] = [0u8; 4];
    let _temp = image_buf.read(&mut u32buf0);
    let _temp = label_buf.read(&mut u32buf1);
    if u32buf0 != [0u8, 0, 8, 3] || u32buf1 != [0u8, 0, 8, 1] {
        return Err(io::Error::new(io::ErrorKind::Other, "read_mnist_images: magic numbers do not match."))
    }
    let _temp = image_buf.read(&mut u32buf0);
    // i don't know how to get int out of 4 bytes using a rust library.
    let num_image_items = bytes_to_usize(&u32buf0);
    let _temp = label_buf.read(&mut u32buf1);
    // i don't know how to get int out of 4 bytes using a rust library.
    let num_label_items = bytes_to_usize(&u32buf1);
    if num_image_items != num_label_items {
        return Err(io::Error::new(io::ErrorKind::Other,
                    format!("the number of the images {} do not match with those of the labels {}.",
                            num_image_items, num_label_items)));
    }
    println!("reading {} image and label pairs from {} and {}...",
                num_image_items, imagefile, labelfile);
    let _temp = image_buf.read(&mut u32buf0);
    let row_dimension = bytes_to_usize(&u32buf0);
    if row_dimension != IMAGE_NUM_ROWS {
        return Err(io::Error::new(io::ErrorKind::Other,
                    format!("the row dimension {} is not {}.", row_dimension, IMAGE_NUM_ROWS)))
    }
    let _temp = image_buf.read(&mut u32buf1);
    let col_dimension = bytes_to_usize(&u32buf1);
    if col_dimension != IMAGE_NUM_COLS {
        return Err(io::Error::new(io::ErrorKind::Other,
                    format!("the col dimension {} is not {}.", col_dimension, IMAGE_NUM_COLS)))
    }
    //let mut images: Vec<MnistImage> = Vec::with_capacity(num_image_items as usize);
    let mut images = Vec::with_capacity(num_image_items);
    let image_buf_bytes = image_buf.bytes();
    let mut label_buf_bytes = label_buf.bytes();

    for image_chunk in &image_buf_bytes.chunks_lazy(IMAGE_SIZE) {
        let one_image = image_chunk.map(|x| { x.unwrap() }).collect::<Vec<u8>>();
        let label = label_buf_bytes.next().unwrap().unwrap();
        images.push(MnistImage { data: one_image, label: label });
    }
    /*let mut images = vec![MnistImage { data: vec![0u8; IMAGE_SIZE], label: 0u8 }; num_image_items];
    for one_image in &mut images { //0..num_image_items {
        //let one_image = &mut images[i]; // IMAGE_SIZE * i..IMAGE_SIZE * (i+1)];
        let mut read_sofar = 0;
        loop {
            let temp = image_buf.read(&mut one_image.data);

            match temp {
                Ok(readnow) => {
                    read_sofar += readnow;
                    if(read_sofar < IMAGE_SIZE) {
                        println!("more to read: read this time = {}, read_sofar = {}", readnow, read_sofar);
                    }
                    else {
                        println!("read enough bytes.");
                        break;
                    };
                }
                _ =>
                    println!("image size {:?} is incompatible.", temp)
            }
                //panic!(format!("image size {:?} is incompatible.", temp))
        }
        let mut one_label = [0u8];
        let _temp = label_buf.read(&mut one_label);
        one_image.label = one_label[0];
    }*/
    Ok(images)
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

pub struct MnistImage {
    pub data: Vec<u8>,
    pub label: u8
}

impl fmt::Debug for MnistImage {
    fn fmt(&self, _ : &mut fmt::Formatter) -> Result<(), fmt::Error> {
        println!("MnistImage: data size {}, label {}", IMAGE_SIZE, self.label);
        Ok(())
    }
}

impl Clone for MnistImage {
    fn clone(&self) -> Self {
        MnistImage {
            data: self.data.clone(),
            label: self.label
        }
    }
    fn clone_from(&mut self, source: &Self) {
        self.data = source.data.clone();
        //self.data = *(&source.data).clone();
        // same self.data = *((&source.data).clone());
        self.label = source.label;
    }
}

fn bytes_to_usize(bytes: &[u8; 4]) -> usize {
    ((bytes[0] as usize) << 24) +
        ((bytes[1] as usize) << 16) +
        ((bytes[2] as usize) << 8) +
        (bytes[3] as usize)
}
