pub mod bench;

#[cfg(test)]
mod tests {
    use crate::bench::MoleculeStore;

    #[test]
    fn get_dde() {
        let store = MoleculeStore::from_json("testfiles/store.json").unwrap();
        let ff = "openff-2.1.0.offxml";
        store.get_dde(ff);
    }
}
