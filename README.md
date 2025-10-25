# 🌾 Precision Agriculture Analysis Toolkit

A comprehensive Python toolkit for precision agriculture analysis using satellite imagery and remote sensing data. This project provides both **analysis libraries** and **educational notebooks** to help farmers, researchers, and agronomists make data-driven decisions.

## 🎯 Project Vision

This toolkit aims to provide accessible, open-source tools for:
- 📊 **Agricultural Monitoring**: Track crop health, growth stages, and field conditions
- 🛰️ **Remote Sensing Analysis**: Process and analyze satellite imagery at scale
- 💧 **Resource Management**: Optimize irrigation, fertilization, and other inputs
- 📈 **Yield Prediction**: Forecast crop yields using multi-temporal analysis
- 🌍 **Environmental Impact**: Monitor soil health, carbon sequestration, and sustainability metrics

## 🚀 Current Capabilities

### NDVI Analysis (Available Now)
Compute vegetation indices from free satellite data (Sentinel-2, Landsat) to monitor crop health and vegetation vigor.

**Key Features:**
- Automatic data retrieval via STAC APIs
- Support for multiple AOI formats (GeoJSON, bounding boxes, coordinates)
- Built-in visualizations and statistics
- Both CLI and Python library interfaces
- Interactive Jupyter tutorial for learning

## 🔮 Coming Soon

- **Multi-spectral Indices**: EVI, SAVI, NDWI, and more
- **Time Series Analysis**: Track changes over growing seasons
- **Yield Modeling**: Predictive analytics for crop production
- **Soil Moisture Estimation**: Remote sensing-based soil monitoring
- **Field Boundary Detection**: Automated field delineation
- **Crop Classification**: Machine learning-based crop type identification

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/chris/precision-ag.git
cd precision-ag

# Create and activate a virtual environment (recommended)
# Option 1: Using uv (fastest)
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Option 2: Using virtualenv
virtualenv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install the package
pip install -e .

# Or with optional dependencies for development and notebooks
pip install -e ".[dev,notebook]"
```

### NDVI Analysis Examples

**Command Line:**
```bash
# Analyze a field using a bounding box
ndvi-compute \
    --aoi "-121.0,37.0,-120.5,37.5" \
    --start-date 2024-06-01 \
    --end-date 2024-06-30

# Or use a GeoJSON file
ndvi-compute --aoi my_field.geojson --start-date 2024-06-01 --end-date 2024-06-30
```

**Python Library:**
```python
from precision_ag.ndvi_from_aoi import compute_ndvi_for_aoi

# Quick NDVI analysis
results = compute_ndvi_for_aoi(
    aoi_input=[-121.0, 37.0, -120.5, 37.5],
    start_date="2024-06-01",
    end_date="2024-06-30"
)
```

### Interactive Tutorials

Launch the Jupyter notebooks to learn interactively:

```bash
jupyter notebook notebooks/NDVI_Tutorial.ipynb
```

For detailed NDVI usage, see the [NDVI module documentation](precision-ag/ndvi_from_aoi.py) or run `ndvi-compute --help`.

## 📁 Project Structure

```
precision-ag/
├── precision-ag/            # Analysis libraries
│   ├── ndvi_from_aoi.py    # NDVI computation module
│   └── ...                  # More analysis modules coming soon
├── notebooks/               # Educational tutorials
│   ├── NDVI_Tutorial.ipynb # NDVI learning notebook
│   └── ...                  # More tutorials coming soon
├── pyproject.toml           # Project configuration
├── README.md                # This file
└── output/                  # Generated outputs (auto-created)
```

## 🔧 Requirements

- Python >= 3.8
- numpy >= 1.20.0
- rasterio >= 1.3.0
- matplotlib >= 3.5.0
- pystac-client >= 0.7.0
- planetary-computer >= 1.0.0
- requests >= 2.28.0
- shapely >= 2.0.0

## 🌍 Data Sources

This tool uses free, public satellite data via STAC APIs:

| Satellite | Resolution | Revisit | Coverage | STAC Catalog |
|-----------|-----------|---------|----------|--------------|
| **Sentinel-2** | 10m | 5 days | Global | [Microsoft Planetary Computer](https://planetarycomputer.microsoft.com/) |
| **Landsat 8/9** | 30m | 16 days | Global | [Earth Search AWS](https://earth-search.aws.element84.com/) |

## 💡 Use Cases

This toolkit supports various precision agriculture applications:

- 🌾 **Crop Health Monitoring**: Detect stress, disease, and nutrient deficiencies early
- 💧 **Irrigation Management**: Optimize water usage based on vegetation and soil data
- 🌱 **Growth Stage Tracking**: Monitor crop development throughout the season
- 📊 **Yield Forecasting**: Predict harvest outcomes using multi-temporal analysis
- 🗺️ **Field Mapping**: Delineate management zones for variable rate applications
- 🌍 **Sustainability Reporting**: Track environmental metrics and carbon footprint
- 🔬 **Research & Development**: Support agricultural research with reproducible analysis

## 🛠️ Development

### Running Tests
```bash
pytest
```

### Code Formatting
```bash
black precision-ag/
```

### Type Checking
```bash
mypy precision-ag/
```

## 📝 Outputs

Analysis tools generate georeferenced raster files (GeoTIFF), visualizations (PNG/PDF), and statistical summaries. All outputs are compatible with standard GIS software (QGIS, ArcGIS) and can be used for further analysis or reporting.

## 📚 Learning & Documentation

This project emphasizes both **practical tools** and **educational resources**:

- **Analysis Libraries**: Production-ready Python modules for data processing
- **Tutorial Notebooks**: Interactive Jupyter notebooks explaining concepts, methods, and best practices
- **Documentation**: Inline code documentation and detailed docstrings
- **Examples**: Real-world use cases demonstrating various agricultural scenarios

Whether you're a researcher, agronomist, or developer, you'll find resources suited to your needs.

## 🤝 Contributing

Contributions are welcome! This project is in active development. Areas where you can help:
- Adding new analysis modules (vegetation indices, soil metrics, yield models)
- Creating educational notebooks and tutorials
- Improving documentation and examples
- Bug fixes and performance improvements

Please feel free to submit a Pull Request or open an issue to discuss new features.

## 📄 License

MIT License - see LICENSE file for details

## 🙏 Acknowledgments

- [Microsoft Planetary Computer](https://planetarycomputer.microsoft.com/) - Free satellite data access
- [STAC](https://stacspec.org/) - Standardized geospatial data catalog
- [ESA Copernicus](https://www.copernicus.eu/) - Sentinel missions
- [NASA/USGS](https://www.usgs.gov/landsat-missions) - Landsat program

## 📧 Contact

For questions or issues, please open an issue on GitHub.

---

**Building the future of precision agriculture, one pixel at a time 🌾🛰️**

