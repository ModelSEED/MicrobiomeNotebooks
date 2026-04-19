#!/usr/bin/env python3
"""
Create an Escher map from ASV community flux data.

This script reads a CSV file where:
- Rows are ASVs (Amplicon Sequence Variants)
- Columns are metabolites
- Values are consumption/production fluxes (negative = consumed, positive = produced)

It generates an Escher-compatible JSON file with reactions for each ASV.
"""

import csv
import json
import math
import uuid
import argparse
from pathlib import Path


def generate_map_id():
    """Generate a random map ID similar to Escher's format."""
    return uuid.uuid4().hex[:12]


def read_flux_csv(csv_path):
    """
    Read the flux CSV file and return ASV data.

    Returns:
        tuple: (metabolite_names, asv_data)
            - metabolite_names: list of metabolite column names
            - asv_data: list of dicts with 'name' and 'fluxes' keys
    """
    asv_data = []
    metabolite_names = []

    with open(csv_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        metabolite_names = header[1:]  # Skip first column (ASV name)

        for row in reader:
            if not row or not row[0]:  # Skip empty rows
                continue

            asv_name = row[0]
            fluxes = {}

            for i, val in enumerate(row[1:]):
                if val and val.strip():
                    try:
                        flux = float(val)
                        if flux != 0:
                            fluxes[metabolite_names[i]] = flux
                    except ValueError:
                        pass  # Skip non-numeric values

            # Only include ASVs with actual flux data
            if fluxes:
                asv_data.append({
                    'name': asv_name,
                    'fluxes': fluxes
                })

    return metabolite_names, asv_data


def categorize_metabolites(asv_data):
    """
    Categorize metabolites based on their flux patterns across all ASVs.

    Returns:
        tuple: (only_consumed, only_produced, mixed)
            - only_consumed: list of metabolites that are only consumed (negative flux)
            - only_produced: list of metabolites that are only produced (positive flux)
            - mixed: list of metabolites that are both consumed and produced
    """
    # Track whether each metabolite is consumed or produced
    consumed = set()
    produced = set()

    for asv in asv_data:
        for met_name, flux in asv['fluxes'].items():
            if flux < 0:
                consumed.add(met_name)
            elif flux > 0:
                produced.add(met_name)

    only_consumed = sorted(consumed - produced)
    only_produced = sorted(produced - consumed)
    mixed = sorted(consumed & produced)

    return only_consumed, only_produced, mixed


def create_escher_map(asv_data, map_name="ASV_Community_Fluxes",
                     reaction_spacing=400, metabolite_spacing=80, prefix=""):
    """
    Create a complete Escher map from ASV flux data.

    This creates a proper Escher map where:
    - Each ASV becomes a reaction
    - Metabolite nodes are shared between reactions
    - Reactions are laid out in a grid
    - Metabolites are positioned by their role:
      - Only consumed: left side
      - Only produced: right side
      - Mixed (both): center near ASV nodes

    Args:
        asv_data: list of dicts with 'name' and 'fluxes' keys
        map_name: name for the map
        reaction_spacing: spacing between reaction centers
        metabolite_spacing: spacing for metabolite fan layout
        prefix: optional prefix for reaction IDs

    Returns:
        list: Escher map JSON structure
    """
    nodes = {}
    reactions = {}
    node_id_counter = 0
    segment_id_counter = 0  # Global segment ID counter

    # Categorize metabolites by their flux patterns
    only_consumed, only_produced, mixed = categorize_metabolites(asv_data)

    print(f"  - Only consumed (reactants): {len(only_consumed)}")
    print(f"  - Only produced (products): {len(only_produced)}")
    print(f"  - Mixed (both): {len(mixed)}")

    # Calculate layout dimensions
    num_reactions = len(asv_data)
    num_cols = max(1, int(math.sqrt(num_reactions * 2)))  # Wider than tall
    num_rows = math.ceil(num_reactions / num_cols)

    # Calculate the center of the ASV reaction grid
    grid_width = num_cols * reaction_spacing
    grid_height = num_rows * reaction_spacing
    grid_center_x = 800 + grid_width / 2
    grid_center_y = 400 + grid_height / 2

    # Position metabolites in three regions
    metabolite_nodes = {}  # metabolite_name -> node_id

    # Left side: only consumed (reactants)
    left_x = 200
    consumed_y_start = grid_center_y - (len(only_consumed) * metabolite_spacing) / 2

    for i, met_name in enumerate(only_consumed):
        node_id = str(node_id_counter)
        bigg_id = met_name.replace(" ", "_").replace("-", "_")
        y = consumed_y_start + i * metabolite_spacing

        nodes[node_id] = {
            "node_type": "metabolite",
            "x": left_x,
            "y": y,
            "bigg_id": bigg_id,
            "name": met_name,
            "label_x": left_x - 100,
            "label_y": y,
            "node_is_primary": True
        }
        metabolite_nodes[met_name] = node_id
        node_id_counter += 1

    # Right side: only produced (products)
    right_x = grid_center_x * 2 - 200
    produced_y_start = grid_center_y - (len(only_produced) * metabolite_spacing) / 2

    for i, met_name in enumerate(only_produced):
        node_id = str(node_id_counter)
        bigg_id = met_name.replace(" ", "_").replace("-", "_")
        y = produced_y_start + i * metabolite_spacing

        nodes[node_id] = {
            "node_type": "metabolite",
            "x": right_x,
            "y": y,
            "bigg_id": bigg_id,
            "name": met_name,
            "label_x": right_x + 20,
            "label_y": y,
            "node_is_primary": True
        }
        metabolite_nodes[met_name] = node_id
        node_id_counter += 1

    # Center: mixed metabolites (arranged in a column near reactions)
    center_x = grid_center_x
    mixed_y_start = grid_center_y - (len(mixed) * metabolite_spacing) / 2

    for i, met_name in enumerate(mixed):
        node_id = str(node_id_counter)
        bigg_id = met_name.replace(" ", "_").replace("-", "_")
        y = mixed_y_start + i * metabolite_spacing

        nodes[node_id] = {
            "node_type": "metabolite",
            "x": center_x,
            "y": y,
            "bigg_id": bigg_id,
            "name": met_name,
            "label_x": center_x + 20,
            "label_y": y,
            "node_is_primary": True
        }
        metabolite_nodes[met_name] = node_id
        node_id_counter += 1

    # Create reactions in a grid layout
    # Position reactions around the mixed metabolites in center
    reaction_start_x = 500
    reaction_start_y = 400

    for rxn_idx, asv in enumerate(asv_data):
        col = rxn_idx % num_cols
        row = rxn_idx // num_cols

        # Reaction center position
        center_x = reaction_start_x + col * reaction_spacing
        center_y = reaction_start_y + row * reaction_spacing

        # Create multimarker (input side - left)
        multimarker_in_id = str(node_id_counter)
        nodes[multimarker_in_id] = {
            "node_type": "multimarker",
            "x": center_x - 20,
            "y": center_y
        }
        node_id_counter += 1

        # Create midmarker (center)
        midmarker_id = str(node_id_counter)
        nodes[midmarker_id] = {
            "node_type": "midmarker",
            "x": center_x,
            "y": center_y
        }
        node_id_counter += 1

        # Create multimarker (output side - right)
        multimarker_out_id = str(node_id_counter)
        nodes[multimarker_out_id] = {
            "node_type": "multimarker",
            "x": center_x + 20,
            "y": center_y
        }
        node_id_counter += 1

        # Build segments (using globally unique segment IDs)
        segments = {}

        # Segment: input multimarker -> midmarker
        segments[str(segment_id_counter)] = {
            "from_node_id": multimarker_in_id,
            "to_node_id": midmarker_id,
            "b1": None,
            "b2": None
        }
        segment_id_counter += 1

        # Segment: output multimarker -> midmarker
        segments[str(segment_id_counter)] = {
            "from_node_id": multimarker_out_id,
            "to_node_id": midmarker_id,
            "b1": None,
            "b2": None
        }
        segment_id_counter += 1

        # Build metabolite list and segments
        metabolites = []

        for met_name, flux in asv['fluxes'].items():
            met_node_id = metabolite_nodes[met_name]
            met_node = nodes[met_node_id]

            bigg_id = met_name.replace(" ", "_").replace("-", "_")
            metabolites.append({
                "bigg_id": bigg_id,
                "coefficient": flux
            })

            # Create segment from reaction marker to metabolite
            if flux < 0:
                # Consumed: from multimarker_in (left side)
                from_node = multimarker_in_id
                from_x = center_x - 20
            else:
                # Produced: from multimarker_out (right side)
                from_node = multimarker_out_id
                from_x = center_x + 20

            # Bezier control points for curved line
            met_x = met_node['x']
            met_y = met_node['y']

            segments[str(segment_id_counter)] = {
                "from_node_id": from_node,
                "to_node_id": met_node_id,
                "b1": {
                    "x": from_x + (met_x - from_x) * 0.3,
                    "y": center_y + (met_y - center_y) * 0.3
                },
                "b2": {
                    "x": met_x + (from_x - met_x) * 0.3,
                    "y": met_y + (center_y - met_y) * 0.3
                }
            }
            segment_id_counter += 1

        # Create reaction
        bigg_id = f"{prefix}{asv['name']}" if prefix else asv['name']
        reactions[str(rxn_idx)] = {
            "name": asv['name'],
            "bigg_id": bigg_id,
            "reversibility": False,
            "label_x": center_x,
            "label_y": center_y - 30,
            "gene_reaction_rule": "",
            "genes": [],
            "metabolites": metabolites,
            "segments": segments
        }

    # Calculate canvas size based on actual node positions
    all_x = [n['x'] for n in nodes.values()]
    all_y = [n['y'] for n in nodes.values()]
    min_x, max_x = min(all_x), max(all_x)
    min_y, max_y = min(all_y), max(all_y)

    canvas_width = max(2000, max_x - min_x + 800)
    canvas_height = max(1500, max_y - min_y + 400)

    # Create the Escher map structure
    escher_map = [
        {
            "map_name": map_name,
            "map_id": generate_map_id(),
            "map_description": f"ASV community flux map with {len(asv_data)} reactions",
            "homepage": "https://escher.github.io",
            "schema": "https://escher.github.io/escher/jsonschema/1-0-0#"
        },
        {
            "reactions": reactions,
            "nodes": nodes,
            "text_labels": {},
            "canvas": {
                "x": min_x - 200,
                "y": min_y - 200,
                "width": canvas_width,
                "height": canvas_height
            }
        }
    ]

    return escher_map


def save_html(escher_map, output_path):
    """Save the Escher map as a standalone HTML file."""
    html = '''<!DOCTYPE html>
<html>
<head>
    <title>ASV Escher Map</title>
    <script src="https://d3js.org/d3.v4.min.js"></script>
    <script src="escher.min.js"></script>
    <style>
        html, body { margin: 0; padding: 0; width: 100%; height: 100%; }
        #map-container { width: 100%; height: 100%; position: absolute; top: 0; left: 0; }
    </style>
</head>
<body>
    <div id="map-container"></div>
    <script>
        var mapData = MAP_DATA_PLACEHOLDER;

        document.addEventListener('DOMContentLoaded', function() {
            var builder = escher.Builder(mapData, null, null,
                d3.select('#map-container'), {
                    menu: 'all',
                    fill_screen: true,
                    scroll_behavior: 'zoom',
                    enable_editing: true,
                    enable_keys: true,
                    enable_search: true
                });

            if (builder.map) {
                setTimeout(function() { builder.map.zoom_extent_canvas(); }, 500);
            }
        });
    </script>
</body>
</html>
'''
    html = html.replace('MAP_DATA_PLACEHOLDER', json.dumps(escher_map))

    with open(output_path, 'w') as f:
        f.write(html)


def main():
    parser = argparse.ArgumentParser(
        description="Create an Escher map from ASV community flux CSV data"
    )
    parser.add_argument(
        "input_csv",
        type=str,
        nargs="?",
        default="RC_-1.5_communityFluxes.csv",
        help="Input CSV file with ASV flux data (default: RC_-1.5_communityFluxes.csv)"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Output JSON file path (default: <input_name>_escher.json)"
    )
    parser.add_argument(
        "-n", "--map-name",
        type=str,
        default=None,
        help="Name for the Escher map (default: derived from input filename)"
    )
    parser.add_argument(
        "-p", "--prefix",
        type=str,
        default="",
        help="Prefix for reaction IDs (default: none)"
    )
    parser.add_argument(
        "--html",
        action="store_true",
        help="Also generate a standalone HTML file"
    )

    args = parser.parse_args()

    # Determine file paths
    input_path = Path(args.input_csv)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}")
        return 1

    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.with_suffix("").with_name(
            f"{input_path.stem}_escher.json"
        )

    if args.map_name:
        map_name = args.map_name
    else:
        map_name = input_path.stem.replace("_", " ").replace("-", " ")

    # Read and process data
    print(f"Reading flux data from: {input_path}")
    metabolite_names, asv_data = read_flux_csv(input_path)

    print(f"Found {len(metabolite_names)} metabolites")
    print(f"Found {len(asv_data)} ASVs with non-zero fluxes")

    if not asv_data:
        print("Error: No ASVs with valid flux data found")
        return 1

    # Create the Escher map
    print("Creating Escher map...")
    escher_map = create_escher_map(
        asv_data,
        map_name=map_name,
        prefix=args.prefix
    )

    # Get stats
    num_nodes = len(escher_map[1]["nodes"])
    num_reactions = len(escher_map[1]["reactions"])
    canvas = escher_map[1]["canvas"]
    print(f"Created map with {num_reactions} reactions and {num_nodes} nodes")
    print(f"Canvas size: {canvas['width']}x{canvas['height']}")

    # Write JSON output
    print(f"Writing Escher map to: {output_path}")
    with open(output_path, 'w') as f:
        json.dump(escher_map, f, indent=2)

    # Write HTML if requested
    if args.html:
        html_path = output_path.with_suffix('.html')
        print(f"Writing HTML to: {html_path}")
        save_html(escher_map, html_path)

    print("Done!")
    return 0


if __name__ == "__main__":
    exit(main())
