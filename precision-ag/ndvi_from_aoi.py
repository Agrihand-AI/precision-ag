#!/usr/bin/env python3
"""
NDVI Computation from Satellite Data using PySTAC

This script allows you to provide an Area of Interest (AOI) and automatically
pulls satellite data using STAC API to compute NDVI (Normalized Difference Vegetation Index).

NDVI = (NIR - Red) / (NIR + Red)

Usage:
    python ndvi_from_aoi.py --aoi <geojson_file> --start-date YYYY-MM-DD --end-date YYYY-MM-DD
    
    Or use it as a module:
    from ndvi_from_aoi import compute_ndvi_for_aoi
"""

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import rasterio
from rasterio.mask import mask
from rasterio.warp import transform_geom
from rasterio.errors import WindowError
import matplotlib.pyplot as plt
from pystac_client import Client
import planetary_computer
import requests
from shapely.geometry import shape, box


class NDVIComputer:
    """Class to handle NDVI computation from satellite data using STAC."""
    
    # STAC Catalog URLs
    PLANETARY_COMPUTER_CATALOG = "https://planetarycomputer.microsoft.com/api/stac/v1"
    EARTH_SEARCH_CATALOG = "https://earth-search.aws.element84.com/v1"
    
    def __init__(self, catalog_url: str = None, collection: str = "sentinel-2-l2a"):
        """
        Initialize the NDVI Computer.
        
        Args:
            catalog_url: URL of the STAC catalog to use. Defaults to Planetary Computer.
            collection: Satellite collection to use. Options:
                - 'sentinel-2-l2a' (default, 10m resolution)
                - 'landsat-c2-l2' (30m resolution)
        """
        self.catalog_url = catalog_url or self.PLANETARY_COMPUTER_CATALOG
        self.collection = collection
        self.catalog = Client.open(self.catalog_url)
        
    def load_aoi(self, aoi_input: Union[str, Path, Dict, List]) -> Dict:
        """
        Load Area of Interest from various input formats.
        
        Args:
            aoi_input: Can be:
                - Path to GeoJSON file
                - GeoJSON dict
                - Bounding box as [min_lon, min_lat, max_lon, max_lat]
                - List of coordinate tuples [(lon1, lat1), (lon2, lat2), ...]
                
        Returns:
            GeoJSON geometry dict
        """
        if isinstance(aoi_input, (str, Path)):
            # Load from file
            with open(aoi_input, 'r') as f:
                geojson = json.load(f)
            
            # Extract geometry if it's a Feature or FeatureCollection
            if geojson.get('type') == 'Feature':
                return geojson['geometry']
            elif geojson.get('type') == 'FeatureCollection':
                # Use the first feature's geometry
                return geojson['features'][0]['geometry']
            else:
                return geojson
                
        elif isinstance(aoi_input, dict):
            # Already a geometry dict
            return aoi_input
            
        elif isinstance(aoi_input, list):
            # Check if it's a bounding box (4 numbers)
            if len(aoi_input) == 4 and all(isinstance(x, (int, float)) for x in aoi_input):
                # Convert bounding box to GeoJSON polygon
                min_lon, min_lat, max_lon, max_lat = aoi_input
                return {
                    "type": "Polygon",
                    "coordinates": [[
                        [min_lon, min_lat],
                        [max_lon, min_lat],
                        [max_lon, max_lat],
                        [min_lon, max_lat],
                        [min_lon, min_lat]
                    ]]
                }
            # Check if it's a list of tuples/lists (polygon coordinates)
            elif len(aoi_input) >= 3 and all(
                isinstance(coord, (tuple, list)) and len(coord) >= 2 
                for coord in aoi_input
            ):
                # Convert list of tuples to GeoJSON polygon
                coords = [[float(coord[0]), float(coord[1])] for coord in aoi_input]
                
                # Close the polygon if not already closed
                if coords[0] != coords[-1]:
                    coords.append(coords[0])
                
                return {
                    "type": "Polygon",
                    "coordinates": [coords]
                }
            else:
                raise ValueError(
                    "List input must be either a bounding box [min_lon, min_lat, max_lon, max_lat] "
                    "or a list of coordinate tuples [(lon1, lat1), (lon2, lat2), ...]"
                )
        else:
            raise ValueError(
                "AOI input must be a file path, GeoJSON dict, bounding box list, "
                "or list of coordinate tuples"
            )
    
    def search_satellite_data(
        self,
        aoi_geometry: Dict,
        start_date: str,
        end_date: str,
        max_cloud_cover: float = 20
    ) -> List:
        """
        Search for satellite imagery covering the AOI.
        
        Args:
            aoi_geometry: GeoJSON geometry dict
            start_date: Start date in YYYY-MM-DD format
            end_date: End date in YYYY-MM-DD format
            max_cloud_cover: Maximum cloud cover percentage (0-100)
            
        Returns:
            List of STAC items
        """
        # Create search query
        search = self.catalog.search(
            collections=[self.collection],
            intersects=aoi_geometry,
            datetime=f"{start_date}/{end_date}",
            query={"eo:cloud_cover": {"lt": max_cloud_cover}}
        )
        
        items = list(search.items())
        print(f"Found {len(items)} scenes matching criteria")
        
        return items
    
    def get_band_assets(self, item) -> Tuple[str, str]:
        """
        Get the NIR and Red band asset keys for the given item.
        
        Args:
            item: STAC item
            
        Returns:
            Tuple of (nir_key, red_key)
        """
        if 'sentinel-2' in self.collection.lower():
            # Sentinel-2 bands
            nir_key = 'B08'  # NIR band
            red_key = 'B04'  # Red band
        elif 'landsat' in self.collection.lower():
            # Landsat bands
            nir_key = 'nir08'  # NIR band
            red_key = 'red'    # Red band
        else:
            raise ValueError(f"Unknown collection: {self.collection}")
            
        return nir_key, red_key
    
    def compute_ndvi(
        self,
        item,
        aoi_geometry: Dict,
        output_path: Optional[str] = None
    ) -> Tuple[np.ndarray, Dict]:
        """
        Compute NDVI from a STAC item.
        
        Args:
            item: STAC item
            aoi_geometry: GeoJSON geometry to clip to
            output_path: Optional path to save the NDVI raster
            
        Returns:
            Tuple of (ndvi_array, metadata_dict)
        """
        # Sign the assets if using Planetary Computer
        if 'planetarycomputer' in self.catalog_url:
            item = planetary_computer.sign(item)
        
        # Get band asset keys
        nir_key, red_key = self.get_band_assets(item)
        
        # Get the href (URL) for each band
        nir_href = item.assets[nir_key].href
        red_href = item.assets[red_key].href
        
        print(f"Processing scene: {item.id}")
        print(f"Date: {item.datetime}")
        
        # Read and clip the bands
        aoi_shape = [aoi_geometry]
        
        try:
            with rasterio.open(nir_href) as nir_src:
                # Transform AOI to match the raster's CRS if needed
                aoi_transformed = aoi_geometry
                if aoi_geometry.get('crs') != nir_src.crs or 'crs' not in aoi_geometry:
                    # Assume AOI is in EPSG:4326 (WGS84) if no CRS specified
                    aoi_transformed = transform_geom(
                        'EPSG:4326',
                        nir_src.crs,
                        aoi_geometry
                    )
                
                # Check if AOI intersects with raster bounds
                raster_bounds = box(*nir_src.bounds)
                aoi_shapely = shape(aoi_transformed)
                
                if not raster_bounds.intersects(aoi_shapely):
                    raise ValueError(
                        f"AOI does not intersect with this scene's coverage area.\n"
                        f"AOI bounds: {aoi_shapely.bounds}\n"
                        f"Scene bounds: {nir_src.bounds}\n"
                        f"Scene CRS: {nir_src.crs}"
                    )
                
                # Clip to AOI
                aoi_shape_transformed = [aoi_transformed]
                nir_clipped, nir_transform = mask(nir_src, aoi_shape_transformed, crop=True)
                nir_clipped = nir_clipped[0].astype(float)  # Get first band
                
                # Store metadata
                meta = nir_src.meta.copy()
                meta.update({
                    "driver": "GTiff",
                    "height": nir_clipped.shape[0],
                    "width": nir_clipped.shape[1],
                    "transform": nir_transform,
                    "count": 1,
                    "dtype": 'float32'
                })
            
            with rasterio.open(red_href) as red_src:
                # Transform AOI for red band as well
                aoi_transformed_red = aoi_geometry
                if aoi_geometry.get('crs') != red_src.crs or 'crs' not in aoi_geometry:
                    aoi_transformed_red = transform_geom(
                        'EPSG:4326',
                        red_src.crs,
                        aoi_geometry
                    )
                aoi_shape_transformed_red = [aoi_transformed_red]
                red_clipped, _ = mask(red_src, aoi_shape_transformed_red, crop=True)
                red_clipped = red_clipped[0].astype(float)
        
        except ValueError as e:
            if "do not overlap" in str(e) or "does not intersect" in str(e):
                print(f"\n⚠️  Warning: {str(e)}")
                print("   This scene doesn't cover your AOI. Skipping...")
                raise
            else:
                raise
        
        # Compute NDVI
        # Handle division by zero
        denominator = nir_clipped + red_clipped
        ndvi = np.where(
            denominator != 0,
            (nir_clipped - red_clipped) / denominator,
            np.nan
        )
        
        # Clip NDVI to valid range [-1, 1]
        ndvi = np.clip(ndvi, -1, 1)
        
        # Save if output path provided
        if output_path:
            with rasterio.open(output_path, 'w', **meta) as dst:
                dst.write(ndvi.astype('float32'), 1)
            print(f"NDVI saved to: {output_path}")
        
        return ndvi, meta
    
    def visualize_ndvi(
        self,
        ndvi: np.ndarray,
        title: str = "NDVI",
        save_path: Optional[str] = None
    ):
        """
        Visualize NDVI array.
        
        Args:
            ndvi: NDVI array
            title: Plot title
            save_path: Optional path to save the figure
        """
        plt.figure(figsize=(12, 8))
        
        # Use RdYlGn colormap (red-yellow-green)
        im = plt.imshow(ndvi, cmap='RdYlGn', vmin=-1, vmax=1)
        plt.colorbar(im, label='NDVI', shrink=0.8)
        plt.title(title, fontsize=14, fontweight='bold')
        plt.xlabel('Column')
        plt.ylabel('Row')
        
        # Add statistics
        valid_ndvi = ndvi[~np.isnan(ndvi)]
        stats_text = f"Mean: {valid_ndvi.mean():.3f}\nStd: {valid_ndvi.std():.3f}\nMin: {valid_ndvi.min():.3f}\nMax: {valid_ndvi.max():.3f}"
        plt.text(0.02, 0.98, stats_text,
                transform=plt.gca().transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                fontsize=10)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {save_path}")
        
        plt.show()


def compute_ndvi_for_aoi(
    aoi_input: Union[str, Path, Dict, List],
    start_date: str,
    end_date: str,
    output_dir: str = "output",
    max_cloud_cover: float = 20,
    collection: str = "sentinel-2-l2a",
    visualize: bool = True
) -> List[Tuple[str, np.ndarray]]:
    """
    High-level function to compute NDVI for an AOI.
    
    Args:
        aoi_input: AOI as file path, GeoJSON dict, or bounding box
        start_date: Start date (YYYY-MM-DD)
        end_date: End date (YYYY-MM-DD)
        output_dir: Directory to save outputs
        max_cloud_cover: Maximum cloud cover percentage
        collection: Satellite collection to use
        visualize: Whether to create visualizations
        
    Returns:
        List of tuples (scene_id, ndvi_array)
    """
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Initialize computer
    computer = NDVIComputer(collection=collection)
    
    # Load AOI
    print("Loading AOI...")
    aoi_geometry = computer.load_aoi(aoi_input)
    
    # Search for data
    print(f"\nSearching for {collection} data from {start_date} to {end_date}...")
    items = computer.search_satellite_data(
        aoi_geometry,
        start_date,
        end_date,
        max_cloud_cover
    )
    
    if not items:
        print("No scenes found matching criteria!")
        return []
    
    # Process each scene
    results = []
    skipped = 0
    
    for i, item in enumerate(items):
        scene_id = item.id
        output_file = output_path / f"ndvi_{scene_id}.tif"
        
        print(f"\n--- Processing scene {i+1}/{len(items)} ---")
        
        try:
            # Compute NDVI
            ndvi, meta = computer.compute_ndvi(
                item,
                aoi_geometry,
                output_path=str(output_file)
            )
            
            results.append((scene_id, ndvi))
            
            # Visualize
            if visualize:
                fig_path = output_path / f"ndvi_{scene_id}.png"
                computer.visualize_ndvi(
                    ndvi,
                    title=f"NDVI - {scene_id}\n{item.datetime.strftime('%Y-%m-%d')}",
                    save_path=str(fig_path)
                )
        
        except ValueError as e:
            if "does not intersect" in str(e) or "do not overlap" in str(e):
                # Scene doesn't cover AOI, skip it
                skipped += 1
                continue
            else:
                # Other error, re-raise
                raise
        except Exception as e:
            print(f"\n⚠️  Error processing scene: {str(e)}")
            skipped += 1
            continue
    
    print(f"\n✓ Successfully processed {len(results)} scenes")
    if skipped > 0:
        print(f"⚠️  Skipped {skipped} scenes (no AOI coverage or errors)")
    print(f"✓ Outputs saved to: {output_path}")
    
    if len(results) == 0:
        print("\n❌ No scenes successfully processed!")
        print("   Tip: Try a different date range or check your AOI coordinates")
    
    return results


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description="Compute NDVI from satellite data using STAC API"
    )
    parser.add_argument(
        '--aoi',
        required=True,
        help='Area of Interest: path to GeoJSON file or bounding box as "min_lon,min_lat,max_lon,max_lat"'
    )
    parser.add_argument(
        '--start-date',
        required=True,
        help='Start date in YYYY-MM-DD format'
    )
    parser.add_argument(
        '--end-date',
        required=True,
        help='End date in YYYY-MM-DD format'
    )
    parser.add_argument(
        '--output-dir',
        default='output',
        help='Output directory for NDVI rasters and figures (default: output)'
    )
    parser.add_argument(
        '--max-cloud-cover',
        type=float,
        default=20,
        help='Maximum cloud cover percentage (default: 20)'
    )
    parser.add_argument(
        '--collection',
        default='sentinel-2-l2a',
        choices=['sentinel-2-l2a', 'landsat-c2-l2'],
        help='Satellite collection to use (default: sentinel-2-l2a)'
    )
    parser.add_argument(
        '--no-visualize',
        action='store_true',
        help='Skip visualization'
    )
    
    args = parser.parse_args()
    
    # Parse AOI
    if args.aoi.endswith('.geojson') or args.aoi.endswith('.json'):
        aoi_input = args.aoi
    else:
        # Try to parse as bounding box
        try:
            bbox = [float(x.strip()) for x in args.aoi.split(',')]
            if len(bbox) != 4:
                raise ValueError
            aoi_input = bbox
        except:
            parser.error("AOI must be a GeoJSON file path or bounding box as 'min_lon,min_lat,max_lon,max_lat'")
    
    # Run NDVI computation
    compute_ndvi_for_aoi(
        aoi_input=aoi_input,
        start_date=args.start_date,
        end_date=args.end_date,
        output_dir=args.output_dir,
        max_cloud_cover=args.max_cloud_cover,
        collection=args.collection,
        visualize=not args.no_visualize
    )


if __name__ == "__main__":
    main()

