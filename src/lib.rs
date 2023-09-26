pub mod bench;

mod utils {
    /// compute the minimum of `v` and then subtract it from all of the elements
    #[inline]
    pub(crate) fn make_rel(v: &mut [f64]) {
        let min = *v.iter().min_by(|x, y| x.total_cmp(y)).unwrap();
        for e in v.iter_mut() {
            *e -= min;
        }
    }
}

#[cfg(test)]
mod tests;
