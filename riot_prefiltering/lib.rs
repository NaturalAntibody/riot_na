mod model;
mod prefiltering;

#[cfg(test)]
mod tests;

use model::{GeneMatch, GeneSegment, PrefilteringResult, PrefilteringSegmentResult, SegmentMatch};
use prefiltering::Prefiltering;
use pyo3::prelude::*;

#[pymodule]
fn riot_na(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Prefiltering>()?;
    m.add_class::<PrefilteringResult>()?;
    m.add_class::<PrefilteringSegmentResult>()?;
    m.add_class::<GeneMatch>()?;
    m.add_class::<GeneSegment>()?;
    m.add_class::<SegmentMatch>()?;
    Ok(())
}
