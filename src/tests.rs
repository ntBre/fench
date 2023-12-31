use std::fs::read_to_string;

use approx::assert_abs_diff_eq;

use crate::bench::MoleculeStore;

/// load (String, f64) pairs from filename
fn load_pairs(filename: &'static str) -> Vec<(String, f64)> {
    let s = read_to_string(filename).unwrap();
    let want: Vec<_> = s
        .lines()
        .map(|line| {
            let fields: Vec<_> = line.split_ascii_whitespace().collect();
            (fields[0].to_owned(), fields[1].parse::<f64>().unwrap())
        })
        .collect();
    want
}

fn check(got: Vec<(String, f64)>, want: Vec<(String, f64)>, eps: f64) {
    let (gr, gv): (Vec<_>, Vec<_>) = got.into_iter().unzip();
    let (wr, wv): (Vec<_>, Vec<_>) = want.into_iter().unzip();
    assert_eq!(gr, wr);
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

#[test]
fn get_dde() {
    let store = MoleculeStore::from_json("testfiles/store.json").unwrap();
    let ff = "openff-2.1.0.offxml";
    let got = store.get_dde(ff);
    let want = load_pairs("testfiles/dde.txt");
    assert_eq!(got.len(), want.len());
    check(got, want, 0.63);
}

#[test]
fn get_rmsd() {
    let store = MoleculeStore::from_json("testfiles/store.json").unwrap();
    let ff = "openff-2.1.0.offxml";
    let got = store.get_rmsd(ff);
    let want = load_pairs("testfiles/rmsd.txt");
    assert_eq!(got.len(), want.len());
    check(got, want, 1.5);
}

#[test]
fn get_tfd() {
    let store = MoleculeStore::from_json("testfiles/store.json").unwrap();
    let ff = "openff-2.1.0.offxml";
    let got = store.get_tfd(ff);
    let want = load_pairs("testfiles/tfd.txt");
    assert_eq!(got.len(), want.len());
    check(got, want, 0.48);
}
