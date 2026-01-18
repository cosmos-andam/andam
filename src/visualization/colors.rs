//! Color schemes for cosmological visualizations

use plotters::style::RGBColor;

/// CMB color scheme (blue to red)
pub fn cmb_colormap(value: f64) -> RGBColor {
    // value should be in [0, 1]
    let v = value.clamp(0.0, 1.0);

    if v < 0.5 {
        let t = v * 2.0;
        RGBColor((255.0 * (1.0 - t)) as u8, 0, (255.0 * t) as u8)
    } else {
        let t = (v - 0.5) * 2.0;
        RGBColor((255.0 * t) as u8, 0, (255.0 * (1.0 - t)) as u8)
    }
}

/// Viridis-like colormap
pub fn viridis(value: f64) -> RGBColor {
    let v = value.clamp(0.0, 1.0);

    let r = (253.0 * v.powf(0.8)) as u8;
    let g = (231.0 * (1.0 - (1.0 - v).powf(2.0))) as u8;
    let b = (37.0 + 218.0 * v) as u8;

    RGBColor(r, g, b)
}
