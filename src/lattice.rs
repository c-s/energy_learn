extern crate rand;
extern crate num;

use self::rand::Rng;
use self::rand::distributions::{IndependentSample, Range, Normal};
use std::ops::Index;
use std::cmp;
use self::num::integer;
use std::sync::mpsc;
use ndarray::Array2;
use std::ops::Neg;

#[derive(Debug)]
pub struct Lattice {
    cube: Cube,
    membranes: Vec<Membrane>,
}

trait Spin {
    fn to_fval(self) -> f64;
    fn to_ival(self) -> i64;
}

impl Spin for bool {
    fn to_fval(self) -> f64 {
        if self { 1.0 } else { -1.0 }
    }
    fn to_ival(self) -> i64 {
        if self { 1 } else { -1 }
    }
}

fn spin_gap() -> i64 { 2 }

trait Sum {
    fn temp_sum(&self) -> i64;
}

impl Sum for Vec<i64> {
    fn temp_sum(&self) -> i64 {
        let mut sum = 0i64;
        for item in self.iter() {
            sum += *item;
        }
        sum
    }
}

impl Lattice {
    pub fn new(cube: Cube, membranes: Vec<Membrane>) -> Self {
        Lattice {
            cube: cube,
            membranes: membranes,
        }
    }

    fn test_and_update_spin<T: Rng>(&mut self,
                            rng: &mut T,
                            flat_spin_range: &Range<usize>,
                            unit_interval: &Range<f64>,
                            free_spinset_start_index: usize,
                            free_membrane: Membrane,
                            inverse_T: f64) {
        let flat_spin_index = flat_spin_range.ind_sample(rng);
        let spinset_index = self.cube.get_spinset_index(flat_spin_index);
        let cond1 = !self.membranes.iter().any(|memb| memb.is_in_flat(flat_spin_index, &self.cube));
        if cond1 ||
            (spinset_index >= free_spinset_start_index &&
                free_membrane.is_in_flat(flat_spin_index, &self.cube)) {
                //if spinset_index >= free_spinset_start_index &&
                //free_membrane.is_in_flat(flat_spin_index, &self.cube) {
                //    println!("this is on a membrane in a test set.");
                //}
                let mut energy_change = 0.0;
                let is_x_boundary = self.cube.is_x_boundary(flat_spin_index);
                let is_y_boundary = self.cube.is_y_boundary(flat_spin_index);
                let is_z_boundary = self.cube.is_z_boundary(flat_spin_index);
                let pos = flat_spin_index;
                let weight_index = self.cube.get_weight_index_from_spin_index(pos);

                if !is_x_boundary.0 { energy_change += self.cube.get_flat_weightX(weight_index - self.cube.weight_stride.0) *
                                                        self.cube.get_flat_spin(pos - self.cube.spin_stride.0).to_fval(); }
                if !is_x_boundary.1 { energy_change += self.cube.get_flat_weightX(weight_index) *
                                                        self.cube.get_flat_spin(pos + self.cube.spin_stride.0).to_fval(); }
                if !is_y_boundary.0 { energy_change += self.cube.get_flat_weightY(weight_index - self.cube.weight_stride.1) *
                                                        self.cube.get_flat_spin(pos - self.cube.spin_stride.1).to_fval(); }
                if !is_y_boundary.1 { energy_change += self.cube.get_flat_weightY(weight_index) *
                                                        self.cube.get_flat_spin(pos + self.cube.spin_stride.1).to_fval(); }
                if !is_z_boundary.0 { energy_change += self.cube.get_flat_weightZ(weight_index - self.cube.weight_stride.2) *
                                                        self.cube.get_flat_spin(pos - self.cube.spin_stride.2).to_fval(); }
                if !is_z_boundary.1 { energy_change += self.cube.get_flat_weightZ(weight_index) *
                                                        self.cube.get_flat_spin(pos + self.cube.spin_stride.2).to_fval(); }
                let spin = self.cube.get_flat_spin(flat_spin_index);
                energy_change *= if spin { spin_gap().neg() as f64 } else { spin_gap() as f64 };
                if energy_change <= 0.0 {
                    *self.cube.get_mut_flat_spin(flat_spin_index) = !spin;
                    self.cube.spin_updated(flat_spin_index);
                }
                else {
                    let prob = unit_interval.ind_sample(rng);
                    if prob < (inverse_T * energy_change).neg().exp() {
                        *self.cube.get_mut_flat_spin(flat_spin_index) = !spin;
                        self.cube.spin_updated(flat_spin_index);
                    }
                }
            }
    }

    fn test_and_update_weightX<T: Rng, W: IndependentSample<f64>>(&mut self,
                                      rng: &mut T,
                                      flat_weight_range: &Range<usize>,
                                      unit_interval: &Range<f64>,
                                      weight_dist: &W,
                                      inverse_T: f64) {
        let flat_weight_index = flat_weight_range.ind_sample(rng);
        let curr_weight = self.cube.get_flat_weightX(flat_weight_index);
        let new_weight = weight_dist.ind_sample(rng);
        let is_x_weight_boundary = self.cube.is_x_weight_boundary(flat_weight_index);
        let mut energy_change = 0.0;
        if !is_x_weight_boundary.0 { energy_change += self.cube.spin_product_x[
                                        flat_weight_index - self.cube.weight_stride.0] as f64; }
        if !is_x_weight_boundary.1 { energy_change += self.cube.spin_product_x[flat_weight_index] as f64; }
        energy_change *= new_weight - curr_weight;
        if energy_change <= 0.0 {
            *self.cube.get_mut_flat_weightX(flat_weight_index) = new_weight;
        } else {
            let prob = unit_interval.ind_sample(rng);
            if prob < (inverse_T * energy_change).neg().exp() {
                *self.cube.get_mut_flat_weightX(flat_weight_index) = new_weight;
            }
        }
    }


    fn test_and_update_weightY<T: Rng, W: IndependentSample<f64>>(&mut self,
                                      rng: &mut T,
                                      flat_weight_range: &Range<usize>,
                                      unit_interval: &Range<f64>,
                                      weight_dist: &W,
                                      inverse_T: f64) {
        let flat_weight_index = flat_weight_range.ind_sample(rng);
        let curr_weight = self.cube.get_flat_weightY(flat_weight_index);
        let new_weight = weight_dist.ind_sample(rng);
        let is_y_weight_boundary = self.cube.is_y_weight_boundary(flat_weight_index);
        let mut energy_change = 0.0;
        if !is_y_weight_boundary.0 { energy_change += self.cube.spin_product_y[
                                        flat_weight_index - self.cube.weight_stride.1] as f64; }
        if !is_y_weight_boundary.1 { energy_change += self.cube.spin_product_y[flat_weight_index] as f64; }
        energy_change *= new_weight - curr_weight;
        if energy_change <= 0.0 {
            *self.cube.get_mut_flat_weightY(flat_weight_index) = new_weight;
        } else {
            let prob = unit_interval.ind_sample(rng);
            if prob < (inverse_T * energy_change).neg().exp() {
                *self.cube.get_mut_flat_weightY(flat_weight_index) = new_weight;
            }
        }
    }


    fn test_and_update_weightZ<T: Rng, W: IndependentSample<f64>>(&mut self,
                                      rng: &mut T,
                                      flat_weight_range: &Range<usize>,
                                      unit_interval: &Range<f64>,
                                      weight_dist: &W,
                                      inverse_T: f64) {
        let flat_weight_index = flat_weight_range.ind_sample(rng);
        let curr_weight = self.cube.get_flat_weightZ(flat_weight_index);
        let new_weight = weight_dist.ind_sample(rng);
        let is_z_weight_boundary = self.cube.is_z_weight_boundary(flat_weight_index);
        let mut energy_change = 0.0;
        if !is_z_weight_boundary.0 { energy_change += self.cube.spin_product_z[
                                        flat_weight_index - self.cube.weight_stride.2] as f64; }
        if !is_z_weight_boundary.1 { energy_change += self.cube.spin_product_z[flat_weight_index] as f64; }
        energy_change *= new_weight - curr_weight;
        if energy_change <= 0.0 {
            *self.cube.get_mut_flat_weightZ(flat_weight_index) = new_weight;
        } else {
            let prob = unit_interval.ind_sample(rng);
            if prob < (inverse_T * energy_change).neg().exp() {
                *self.cube.get_mut_flat_weightZ(flat_weight_index) = new_weight;
            }
        }
    }

    pub fn run<T: Rng>(&mut self,
                       rng: &mut T,
                       total_steps: usize,
                       interm_step_size: usize,
                       num_spin_runs_per_weight: usize,
                       free_membrane_index: usize,
                       free_spinset_start_index: usize,
                       observe_index: usize,
                       weight_gaussian_width: f64,
                       temperature: f64,
                       sender: mpsc::Sender<Array2<bool>>) {
        let flat_spin_range = Range::new(0, self.cube.flat_spin_size());
        let flat_weight_range = Range::new(0, self.cube.flat_weight_size());
        let unit_interval = Range::new(0f64, 1f64);
        let normal_dist = Normal::new(0.0, weight_gaussian_width * weight_gaussian_width);
        let inverse_T = 1.0 / temperature;
        let free_membrane = self.membranes[free_membrane_index];

        for index in 0..total_steps {
            if index % interm_step_size == 0 {
                println!("index: {}", index);
                println!("energy: {}", self.cube.calculate_energy());
                //println!("are spin products consistent? {}", self.cube.check_spin_product_consistency());
                //sender.send(self.cube.get_spinset(fixed_spinset_index).to_vec());
                sender.send(self.get_membrane_spins(free_membrane_index, observe_index));
            }
            if index % num_spin_runs_per_weight == 0 {
                self.test_and_update_weightX(
                    rng,
                    &flat_weight_range,
                    &unit_interval,
                    &normal_dist, inverse_T);
                self.test_and_update_weightY(
                    rng,
                    &flat_weight_range,
                    &unit_interval,
                    &normal_dist, inverse_T);
                self.test_and_update_weightZ(
                    rng,
                    &flat_weight_range,
                    &unit_interval,
                    &normal_dist, inverse_T);
            }
            else {
                self.test_and_update_spin(
                    rng,
                    &flat_spin_range,
                    &unit_interval,
                    free_spinset_start_index, free_membrane, inverse_T);
            }
        }
    }

    pub fn update_spin_products(&mut self) {
        self.cube.update_spin_products();
    }

    pub fn get_spin(&self, membrane_index: usize, spinset_index: usize, pos: (usize, usize)) -> bool {
        self.membranes[membrane_index].get_spin(&self.cube, spinset_index, pos)
    }

    pub fn get_mut_spin<'a>(&'a mut self, membrane_index: usize, spinset_index: usize, pos: (usize, usize)) -> &'a mut bool {
        self.membranes[membrane_index].get_mut_spin(&mut self.cube, spinset_index, pos)
    }

   pub fn set_membrane(&mut self, membrane_index: usize, spinset_index: usize, spins: &[bool]) {
       let membrane = &mut self.membranes[membrane_index];
       membrane.set_membrane(&mut self.cube, spinset_index, spins);
   }

   pub fn get_membrane_spins(&self, membrane_index: usize, spinset_index: usize) -> Array2<bool> {
       let membrane = &self.membranes[membrane_index];
       membrane.get_membrane(&self.cube, spinset_index)
   }
}

// each surface of the cube is identified with the opposite one.
// i.e. it is actually a T^3.
// weightX (i, j, k) connects spin (i, j, k) and (i+1, j, k) modulo stride.0
// note the strides cannot be changed arbitrarily: Their order is tied to implementation.
#[derive(Debug)]
pub struct Cube {
    spins: Vec<bool>,
    weightX: Vec<f64>,
    weightY: Vec<f64>,
    weightZ: Vec<f64>,
    size: (usize, usize, usize),
    spinset_stride: usize,
    spin_stride: (usize, usize, usize),
    weight_stride: (usize, usize, usize),

    // auxiliar variable that can be calculated using the above,
    // but store them for performance.
    num_spinsets: usize,
    x_boundary_helper: (usize, usize),
    y_boundary_helper: (usize, usize),
    z_boundary_helper: (usize, usize),
    x_weight_boundary_helper: (usize, usize),
    y_weight_boundary_helper: (usize, usize),
    z_weight_boundary_helper: (usize, usize),
    spin_product_x: Vec<i64>,
    spin_product_y: Vec<i64>,
    spin_product_z: Vec<i64>,
}

impl Cube {
    pub fn new<T: Rng>(rng: &mut T, num_spin_configs: usize, size: (usize, usize, usize)) -> Self {
        let total_spin_size = num_spin_configs * size.0 * size.1 * size.2;
        let total_weight_size = size.0 * size.1 * size.2;
        let spins = rng.gen_iter::<bool>().take(total_spin_size).collect::<Vec<bool>>();
        let spin_stride = (num_spin_configs, num_spin_configs * size.0, num_spin_configs * size.0 * size.1);
        //let spin_product_x = spins.iter().zip(spins.iter().skip(spin_stride.0).cycle()).map(|(s, n)| s*n).sum<i64>();
        let spin_product_x = vec![0i64; total_weight_size];
        let spin_product_y = vec![0i64; total_weight_size];
        let spin_product_z = vec![0i64; total_weight_size];
        /*let mut spin_product_x = Vec::with_capacity(total_weight_size);
        let mut spin_product_y = Vec::with_capacity(total_weight_size);
        let mut spin_product_z = Vec::with_capacity(total_weight_size);
        for ii in 0..total_weight_size {
            let i = ii * num_spin_configs;
            let mut x_sum = 0i64;
            let mut y_sum = 0i64;
            let mut z_sum = 0i64;
            for j in 0..num_spin_configs {
                x_sum += spins[i + j] as i64 * spins[(i + j + spin_stride.0) % total_spin_size] as i64;
                y_sum += spins[i + j] as i64 * spins[(i + j + spin_stride.1) % total_spin_size] as i64;
                z_sum += spins[i + j] as i64 * spins[(i + j + spin_stride.2) % total_spin_size] as i64;
            }
            spin_product_x.push(x_sum);
            spin_product_y.push(y_sum);
            spin_product_z.push(z_sum);
        }*/

        //let mut spins = Vec::with_capacity(total_spin_size);
        //for i in 0..total_spin_size {
        //    spins.push(rng.gen::<bool>());
        //}
        let between = Range::new(-1f64, 1f64);
        let mut weightX = Vec::with_capacity(total_weight_size);
        let mut weightY = Vec::with_capacity(total_weight_size);
        let mut weightZ = Vec::with_capacity(total_weight_size);
        for i in 0..total_weight_size {
            weightX.push(between.ind_sample(rng));
        }
        for i in 0..total_weight_size {
            weightY.push(between.ind_sample(rng));
        }
        for i in 0..total_weight_size {
            weightZ.push(between.ind_sample(rng));
        }
        Cube {
           spins: spins,
           weightX: weightX,
           weightY: weightY,
           weightZ: weightZ,
           size: size,
           spinset_stride: 1,
           spin_stride: spin_stride,
           weight_stride: (1, size.0, size.0 * size.1),
           // auxiliar variable that can be calculated using the above,
           // but store them just in case.
           num_spinsets: num_spin_configs,
           x_boundary_helper: (num_spin_configs, num_spin_configs * size.0),
           y_boundary_helper: (num_spin_configs * size.0,
                                num_spin_configs * size.0 * size.1),
           z_boundary_helper: (num_spin_configs * size.0 * size.1,
                                num_spin_configs * size.0 * size.1 * size.2),
           x_weight_boundary_helper: (1, size.0),
           y_weight_boundary_helper: (size.0, size.0 * size.1),
           z_weight_boundary_helper: (size.0 * size.1, size.0 * size.1 * size.2),
           spin_product_x: spin_product_x,
           spin_product_y: spin_product_y,
           spin_product_z: spin_product_z,
        }
    }
    fn flat_spin_size(&self) -> usize {
        self.spins.len()
    }
    fn flat_weight_size(&self) -> usize {
        self.size.0 * self.size.1 * self.size.2
    }
    //fn update_membrane(membrane: Membrane, spins: Vec<bool>) {
    //    assert_eq!(membrane.total_size(), spins.len());
    //}
    fn spin_updated(&mut self, pos: usize) {
        // update spin_product_x/y/z.
        let is_x_boundary = self.is_x_boundary(pos);
        let is_y_boundary = self.is_y_boundary(pos);
        let is_z_boundary = self.is_z_boundary(pos);
        //println!("pos: {:?}", pos);
        let weight_index = self.get_weight_index_from_spin_index(pos);
        //println!("weight_index: {:?}", weight_index);
        let total_spin_size = self.flat_spin_size();
        let spin = self.get_flat_spin(pos);
        //println!("spin stride0: {:?}, stride1: {:?}, stride2: {:?}", self.spin_stride.0, self.spin_stride.1, self.spin_stride.2);
        //println!("weight stride0: {:?}, stride1: {:?}, stride2: {:?}", self.weight_stride.0, self.weight_stride.1, self.weight_stride.2);
        let spin_chg = if spin { spin_gap() } else { spin_gap().neg() };
        if weight_index >= self.weight_stride.0 {
            let prev_spin = self.get_flat_spin(pos - self.spin_stride.0).to_ival();
            self.spin_product_x[weight_index - self.weight_stride.0] += spin_chg * prev_spin;
        }
        let next_spin = self.get_flat_spin((pos + self.spin_stride.0) % total_spin_size).to_ival();
        self.spin_product_x[weight_index] += spin_chg * next_spin;
        if weight_index >= self.weight_stride.1 {
            let prev_spin = self.get_flat_spin(pos - self.spin_stride.1).to_ival();
            self.spin_product_y[weight_index - self.weight_stride.1] += spin_chg * prev_spin;
        }
        let next_spin = self.get_flat_spin((pos + self.spin_stride.1) % total_spin_size).to_ival();
        self.spin_product_y[weight_index] += spin_chg * next_spin;
        if weight_index >= self.weight_stride.2 {
            let prev_spin = self.get_flat_spin(pos - self.spin_stride.2).to_ival();
            self.spin_product_z[weight_index - self.weight_stride.2] += spin_chg * prev_spin;
        }
        let next_spin = self.get_flat_spin((pos + self.spin_stride.2) % total_spin_size).to_ival();
        self.spin_product_z[weight_index] += spin_chg * next_spin;
    }

    fn get_spin<'a>(&'a self, spinset_index: usize, coords: (usize, usize, usize)) -> &'a bool {
        let pos = self.get_flat_spin_index(spinset_index, coords);
        &self.spins[pos]
    }
    fn get_mut_spin<'a>(&'a mut self, spinset_index: usize, coords: (usize, usize, usize)) -> &'a mut bool {
        let pos = self.get_flat_spin_index(spinset_index, coords);
        &mut self.spins[pos]
    }
    fn get_flat_spin(&self, pos: usize) -> bool { self.spins[pos] }
    fn get_mut_flat_spin<'a>(&'a mut self, pos: usize) -> &'a mut bool { &mut self.spins[pos] }
    fn get_flat_weightX(&self, pos: usize) -> f64 { self.weightX[pos] }
    fn get_flat_weightY(&self, pos: usize) -> f64 { self.weightY[pos] }
    fn get_flat_weightZ(&self, pos: usize) -> f64 { self.weightZ[pos] }
    fn get_mut_flat_weightX<'a>(&'a mut self, pos: usize) -> &'a mut f64 { &mut self.weightX[pos] }
    fn get_mut_flat_weightY<'a>(&'a mut self, pos: usize) -> &'a mut f64 { &mut self.weightY[pos] }
    fn get_mut_flat_weightZ<'a>(&'a mut self, pos: usize) -> &'a mut f64 { &mut self.weightZ[pos] }
    fn get_weightX(&self, pos: (usize, usize, usize)) -> f64 {
        self.weightX[self.get_flat_weight_index(pos)]
    }
    fn get_mut_weightX<'a>(&'a mut self, pos: (usize, usize, usize)) -> &'a mut f64 {
        let flat_pos = self.get_flat_weight_index(pos);
        &mut self.weightX[flat_pos]
    }
    fn get_weightY(&self, pos: (usize, usize, usize)) -> f64 {
        self.weightY[self.get_flat_weight_index(pos)]
    }
    fn get_mut_weightY<'a>(&'a mut self, pos: (usize, usize, usize)) -> &'a mut f64 {
        let flat_pos = self.get_flat_weight_index(pos);
        &mut self.weightY[flat_pos]
    }
    fn get_weightZ(&self, pos: (usize, usize, usize)) -> f64 {
        self.weightZ[self.get_flat_weight_index(pos)]
    }
    fn get_mut_weightZ<'a>(&'a mut self, pos: (usize, usize, usize)) -> &'a mut f64 {
        let flat_pos = self.get_flat_weight_index(pos);
        &mut self.weightZ[flat_pos]
    }
    fn get_spin_from_flat_index(&self, pos: usize) -> bool {
        self.spins[pos]
    }
    fn get_mut_spin_from_flat_index<'a>(&'a mut self, pos: usize) -> &'a mut bool {
        &mut self.spins[pos]
    }

    //fn get_mut_spin_from_membrane<'a>(
    //    &'a mut self, 
    //    spin_index: usize, 
    //    membrane: &Membrane, (x, y): (usize, usize)) -> &'a mut bool {
    //    let pos = membrane.position + membrane.stride.0 * x + membrane.stride.1 * y;
    //    self.get_mut_spin_from_flat_index(spin_index, pos)
    //}

    fn get_flat_spin_index(&self, spinset_index: usize, (x, y, z): (usize, usize, usize)) -> usize {
        spinset_index * self.spinset_stride + x * self.spin_stride.0 + y * self.spin_stride.1 + z * self.spin_stride.2
    }
    fn get_flat_weight_index(&self, (x, y, z): (usize, usize, usize)) -> usize {
        x * self.weight_stride.0 + y * self.weight_stride.1 + z * self.weight_stride.2
    }
    fn get_weight_index_from_spin_index(&self, flat_spin_index: usize) -> usize {
        flat_spin_index / self.num_spinsets
    }
    fn get_spinset_index(&self, flat_spin_index: usize) -> usize {
        flat_spin_index % self.num_spinsets
    }

    fn is_x_boundary(&self, flat_spin_index: usize) -> (bool, bool) {
        let mod_x = flat_spin_index % self.x_boundary_helper.1;
        (mod_x < self.x_boundary_helper.0,
            mod_x >= self.x_boundary_helper.1 - self.x_boundary_helper.0)
    }

    fn is_y_boundary(&self, flat_spin_index: usize) -> (bool, bool) {
        let mod_y = flat_spin_index % self.y_boundary_helper.1;
        (mod_y < self.y_boundary_helper.0,
            mod_y >= self.y_boundary_helper.1 - self.y_boundary_helper.0)
    }

    fn is_z_boundary(&self, flat_spin_index: usize) -> (bool, bool) {
        let mod_z = flat_spin_index % self.z_boundary_helper.1;
        (mod_z < self.z_boundary_helper.0,
            mod_z >= self.z_boundary_helper.1 - self.z_boundary_helper.0)
    }

    fn is_x_weight_boundary(&self, flat_weight_index: usize) -> (bool, bool) {
        let mod_x = flat_weight_index % self.x_weight_boundary_helper.1;
        (mod_x < self.x_weight_boundary_helper.0,
            mod_x >= self.x_weight_boundary_helper.1 - self.x_weight_boundary_helper.0)
    }

    fn is_y_weight_boundary(&self, flat_weight_index: usize) -> (bool, bool) {
        let mod_y = flat_weight_index % self.y_weight_boundary_helper.1;
        (mod_y < self.y_weight_boundary_helper.0,
            mod_y >= self.y_weight_boundary_helper.1 - self.y_weight_boundary_helper.0)
    }

    fn is_z_weight_boundary(&self, flat_weight_index: usize) -> (bool, bool) {
        let mod_z = flat_weight_index % self.z_weight_boundary_helper.1;
        (mod_z < self.z_weight_boundary_helper.0,
            mod_z >= self.z_weight_boundary_helper.1 - self.z_weight_boundary_helper.0)
    }
    fn calculate_energy(&self) -> f64 {
        let mut energy = 0.0;
        for x in 0..(self.size.0 - 1) {
            for y in 0..(self.size.1 - 1) {
                for z in 0..(self.size.2 - 1) {
                    let coord = (x, y, z);
                    let weight_pos = self.get_flat_weight_index(coord);
                    //let x_plus_pos = self.get_flat_weight_index((x + 1, y, z));
                    //let y_plus_pos = self.get_flat_weight_index((x, y + 1, z));
                    //let z_plus_pos = self.get_flat_weight_index((x, y, z + 1));
                    let weightX = self.get_flat_weightX(weight_pos);
                    let weightY = self.get_flat_weightY(weight_pos);
                    let weightZ = self.get_flat_weightZ(weight_pos);
                    let sx = self.spin_product_x[weight_pos] as f64;
                    let sy = self.spin_product_y[weight_pos] as f64;
                    let sz = self.spin_product_z[weight_pos] as f64;
                    energy += sx * weightX + sy * weightY + sz * weightZ;
                }
            }
        }
        energy
    }

    fn check_spin_product_consistency(&self) -> bool {
        let total_weight_size = self.flat_weight_size();
        let total_spin_size = self.flat_spin_size();
        let mut spin_product_x = Vec::with_capacity(total_weight_size);
        let mut spin_product_y = Vec::with_capacity(total_weight_size);
        let mut spin_product_z = Vec::with_capacity(total_weight_size);
        //for i in (0..).step_by(num_spin_configs).take(total_weight_size) {
        let spins = &self.spins;
        for ii in 0..total_weight_size {
            let i = ii * self.num_spinsets;
            let mut x_sum = 0i64;
            let mut y_sum = 0i64;
            let mut z_sum = 0i64;
            for j in 0..self.num_spinsets {
                x_sum += spins[i + j].to_ival() * spins[(i + j + self.spin_stride.0) % total_spin_size].to_ival();
                y_sum += spins[i + j].to_ival() * spins[(i + j + self.spin_stride.1) % total_spin_size].to_ival();
                z_sum += spins[i + j].to_ival() * spins[(i + j + self.spin_stride.2) % total_spin_size].to_ival();
            }
            spin_product_x.push(x_sum);
            spin_product_y.push(y_sum);
            spin_product_z.push(z_sum);
        }
        println!("spin_product_x sum: {}, self.spin_product_x sum: {}", spin_product_x.temp_sum(), self.spin_product_x.temp_sum());
        println!("spin_product_y sum: {}, self.spin_product_y sum: {}", spin_product_y.temp_sum(), self.spin_product_y.temp_sum());
        println!("spin_product_z sum: {}, self.spin_product_z sum: {}", spin_product_z.temp_sum(), self.spin_product_z.temp_sum());
        assert_eq!(spin_product_x.len(), self.spin_product_x.len());
        assert_eq!(spin_product_y.len(), self.spin_product_y.len());
        assert_eq!(spin_product_z.len(), self.spin_product_z.len());
        spin_product_x == self.spin_product_x &&
            spin_product_y == self.spin_product_y &&
                spin_product_z == self.spin_product_z
    }

    fn update_spin_products(&mut self) {
        let total_weight_size = self.flat_weight_size();
        let total_spin_size = self.flat_spin_size();
        let spins = &self.spins;
        for ii in 0..total_weight_size {
            let i = ii * self.num_spinsets;
            let mut x_sum = 0i64;
            let mut y_sum = 0i64;
            let mut z_sum = 0i64;
            for j in 0..self.num_spinsets {
                x_sum += spins[i + j].to_ival() * spins[(i + j + self.spin_stride.0) % total_spin_size].to_ival();
                y_sum += spins[i + j].to_ival() * spins[(i + j + self.spin_stride.1) % total_spin_size].to_ival();
                z_sum += spins[i + j].to_ival() * spins[(i + j + self.spin_stride.2) % total_spin_size].to_ival();
            }
            self.spin_product_x[ii] = x_sum;
            self.spin_product_y[ii] = y_sum;
            self.spin_product_z[ii] = z_sum;
        }
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Membrane {
    spinset_stride: usize,
    spin_stride: (usize, usize),
    weight_stride: (usize, usize),
    size: (usize, usize),

    spin_position: usize, // the first coordinate of the surface in the cube.
    weight_position: usize, // the first coordinate of the surface in the cube.
}

impl Membrane {
    pub fn xy(cube: &Cube, size: (usize, usize), p: (usize, usize, usize)) -> Self {
        Membrane {
            weight_stride: (cube.weight_stride.0, cube.weight_stride.1),
            spin_stride: (cube.spin_stride.0, cube.spin_stride.1),
            spinset_stride: cube.spinset_stride,
            size: size,
            spin_position: cube.get_flat_spin_index(0, p),
            weight_position: cube.get_flat_weight_index(p),
        }
    }
    pub fn yz(cube: &Cube, size: (usize, usize), p: (usize, usize, usize)) -> Self {
        Membrane {
            weight_stride: (cube.weight_stride.1, cube.weight_stride.2),
            spin_stride: (cube.spin_stride.1, cube.spin_stride.2),
            spinset_stride: cube.spinset_stride,
            size: size,
            spin_position: cube.get_flat_spin_index(0, p),
            weight_position: cube.get_flat_weight_index(p),
        }
    }
    pub fn xz(cube: &Cube, size: (usize, usize), p: (usize, usize, usize)) -> Self {
        Membrane {
            weight_stride: (cube.weight_stride.0, cube.weight_stride.2),
            spin_stride: (cube.spin_stride.0, cube.spin_stride.2),
            spinset_stride: cube.spinset_stride,
            size: size,
            spin_position: cube.get_flat_spin_index(0, p),
            weight_position: cube.get_flat_weight_index(p),
        }
    }
    fn total_size(&self) -> usize {
        self.size.0 * self.size.1
    }

    fn is_in(&self, pos: (usize, usize, usize), cube: &Cube) -> bool {
        self.is_in_flat(cube.get_flat_spin_index(0, pos), cube)
    }

    fn is_in_flat(&self, pos: usize, cube: &Cube) -> bool {
        let rel_pos = pos as i64 - self.spin_position as i64;
        //let (index_in_oneset, _) = integer::div_mod_floor(rel_pos, cube.num_spinsets);
        let result = if self.spin_stride.0 <= self.spin_stride.1 {
            let (rel_1, rem_1) = integer::div_mod_floor(rel_pos, self.spin_stride.1 as i64);
            if rel_1 < 0 || rel_1 >= self.size.1 as i64 { false }
            else {
                let (rel_0, _) = integer::div_mod_floor(rem_1, self.spin_stride.0 as i64);
                rel_0 >= 0 && rel_0 < self.size.0 as i64
            }
        } else {
            let (rel_0, rem_0) = integer::div_mod_floor(rel_pos, self.spin_stride.0 as i64);
            if rel_0 < 0 || rel_0 >= self.size.0 as i64 { false }
            else {
                let (rel_1, _) = integer::div_mod_floor(rem_0, self.spin_stride.1 as i64);
                rel_1 >= 0 || rel_1 < self.size.1 as i64
            }
        };
        //println!("pos: {}, rel_pos: {}, spin_stride: {:?}, size: {:?}, is_in_flat? {}", pos, rel_pos, self.spin_stride, self.size,result);
        result
        //let small_stride = cmp::min(self.spin_stride.0, self.spin_stride.1);
        //let (_, rel_pos) = integer::div_rem(rel_pos, cube.spinset_stride);
        //let (d, m) = integer::div_rem(index_in_oneset, small_stride);
        //m == 0 && d>=0 && d<self.size.0 * self.size.1
    }

    fn get_spin(&self, cube: &Cube, spinset_index: usize, pos: (usize, usize)) -> bool {
        let flat_pos = self.get_flat_spin_index(spinset_index, pos);
        cube.get_spin_from_flat_index(flat_pos)
    }

    fn get_mut_spin<'a>(&self, cube: &'a mut Cube, spinset_index: usize, pos: (usize, usize)) -> &'a mut bool {
        let flat_pos = self.get_flat_spin_index(spinset_index, pos);
        cube.get_mut_spin_from_flat_index(flat_pos)
    }

    fn get_flat_spin_index(&self, spinset_index: usize, (x, y): (usize, usize)) -> usize {
        let pos = self.spin_position + spinset_index * self.spinset_stride + x * self.spin_stride.0 + y * self.spin_stride.1;
        //println!("spin_position: {}, spinset_index: {}, stride: {}, pos: {}", self.spin_position, spinset_index, self.spinset_stride, pos);
        pos
    }
    fn get_flat_weight_index(&self, (x, y): (usize, usize)) -> usize {
        self.weight_position + x * self.weight_stride.0 + y * self.weight_stride.1
    }
    pub fn set_membrane(&mut self, cube: &mut Cube, spinset_index: usize, spins: &[bool]) {
       let mut spins_iter = spins.iter();
       for i in 0..self.size.0 {
           for j in 0..self.size.1 {
               if let Some(s) = spins_iter.next() {
                   *self.get_mut_spin(cube, spinset_index, (i, j)) = *s;
               }
           }
       }
    }
    pub fn get_membrane<'a>(&self, cube: &'a Cube, spinset_index: usize) -> Array2<bool> {
        let spins = (0..self.size.1).flat_map(|j| {
            (0..self.size.0).map(move |i| {
                //println!("({}, {})", i, j);
                self.get_spin(cube, spinset_index, (i, j))
            })
        }).collect::<Vec<bool>>();
        Array2::new(spins, (1, self.size.0))
    }
}
