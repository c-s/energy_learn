use std::ops::{Index, IndexMut};

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
            size: if num_steps.0 <= num_steps.1 { (num_steps.1 / num_steps.0, num_steps.1) }
                                            else { (num_steps.0, num_steps.0 / num_steps.1) },
        }
    }
}

