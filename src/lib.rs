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
mod tests {
    use std::fs::read_to_string;

    use approx::assert_abs_diff_eq;

    use crate::bench::MoleculeStore;

    #[test]
    fn get_dde() {
        let store = MoleculeStore::from_json("testfiles/store.json").unwrap();
        let ff = "openff-2.1.0.offxml";
        let got = store.get_dde(ff);
        let s = read_to_string("testfiles/dde.txt").unwrap();
        let want: Vec<_> = s
            .lines()
            .map(|line| {
                let fields: Vec<_> = line.split_ascii_whitespace().collect();
                (fields[0].to_owned(), fields[1].parse::<f64>().unwrap())
            })
            .collect();
        assert_eq!(got.len(), want.len());
        let (gr, gv): (Vec<_>, Vec<_>) = got.into_iter().unzip();
        let (wr, wv): (Vec<_>, Vec<_>) = want.into_iter().unzip();
        assert_eq!(gr, wr);
        let eps = 0.63;
        let mut disagree = 0;
        for (i, (g, w)) in gv.iter().zip(wv.iter()).enumerate() {
            if (g - w).abs() > eps {
                println!("{i:5}{:>12}{g:12.8}{w:12.8}{:12.8}", gr[i], g - w);
                disagree += 1;
            }
        }
        if disagree > 0 {
            println!("{} / {} disagree", disagree, gv.len());
        }
        assert_abs_diff_eq!(gv.as_slice(), wv.as_slice(), epsilon = eps);
    }
}
