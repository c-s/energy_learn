use std::ops::{Index, IndexMut};
use std::fmt;

#[derive(Clone, Debug)]
pub struct Array2<T: Clone> {
    data: Vec<T>,
    stride: (usize, usize),
    size: (usize, usize),
}

#[derive(Clone, Debug)]
struct Array3<T: Clone> {
    data: Vec<T>,
    stride: (usize, usize, usize),
    size: (usize, usize, usize),
}

impl<T: Clone> IndexMut<(usize, usize)> for Array2<T> {
    fn index_mut<'a>(&'a mut self, (x, y): (usize, usize)) -> &'a mut T {
        let pos = self.stride.0 * x + self.stride.1 * y;
        &mut self.data[pos]
    }
}

impl<T: Clone> Index<(usize, usize)> for Array2<T> {
    type Output = T;

    fn index<'a>(&'a self, (x, y): (usize, usize)) -> &'a T {
        let pos = self.stride.0 * x + self.stride.1 * y;
        &self.data[pos]
    }
}

impl<T> Array2<T> where T: Clone {
    pub fn new(v: Vec<T>, stride: (usize, usize)) -> Self {
        let num_steps = (v.len() / stride.0, v.len() / stride.1);

        Array2 {
            data: v,
            stride: stride,
            size: if num_steps.0 <= num_steps.1 { (num_steps.1 / num_steps.0, stride.0) }
                                            else { (stride.1, num_steps.0 / num_steps.1) },
        }
    }
}

impl fmt::Display for Array2<bool> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Array2<bool>\n");
        for i in 0..self.size.0 {
            for j in 0..self.size.1 {
                write!(f, "{}", if self[(i, j)] { '1' } else { '0' });
            }
            write!(f, "\n");
        }
        Ok(())
    }
}
