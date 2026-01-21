# Desolvation Volume Calculator

Calculates AutoDock4 desolvation volumes for atoms in protein structures using SMARTS-based atom typing.

## Overview

For each atom, the desolvation volume is computed as:

**Σ (Vₙ × exp(-0.5 × (distance / σ)²))**

where the sum is over all atoms n within the cutoff distance.

- **Vₙ**: Volume of neighboring atom
- **distance**: Distance between atoms
- **σ**: Gaussian width parameter (default: 3.6 Å)
- **cutoff**: Maximum interaction distance (default: 20.42 Å)

## Required Input Files

### 1. `atom_types.json`
Defines atom types using SMARTS patterns and associated volumes.

```json
{
  "atom_types": {
    "6_0": {
      "smarts": "[#6;D0]",
      "volume": 26.54,
      "description": "Any carbon with 0 bonds"
    },
    "6_1": {
      "smarts": "[#6;D1]",
      "volume": 42.26,
      "description": "Any carbon with 1 bond"
    },
    "6_2": {
      "smarts": "[#6;D2]",
      "volume": 25.35,
      "description": "Any carbon with 2 bonds"
    }
  }
}
```

**Fields:**
- `smarts`: SMARTS pattern for atom matching
- `volume`: Atomic volume (Å³) used in desolvation calculation

### 2. `pdb_sources.json`
Specifies PDB files to process and optional custom residue templates.

```json
{
  "sources": {
    "protein_complex": [
      "/path/to/structure1.pdb",
      "/path/to/structure2.pdb"
    ],
    "ligand_bound": [
      "/path/to/structure3.pdb"
    ]
  },
  "custom_templates": {
    "sdf_templates": [
      {
        "resname": "LIG",
        "sdf_path": "/path/to/ligand.sdf"
      }
    ],
    "smiles_templates": [
      {
        "resname": "POP",
        "smiles": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC"
      }
    ]
  }
}
```

**Fields:**
- `sources`: Dictionary mapping source names to lists of PDB file paths
- `custom_templates` (optional): Residue templates for non-standard residues
  - `sdf_templates`: sdf templates for meeko
  - `smiles_templates`: smiles to make templates for meeko

## Usage

```python
from desolvation_volume_calculator import DesolvationCalculatorMeeko

calculator = DesolvationCalculatorMeeko(
    atom_types_file='atom_types.json',
    pdb_sources_file='pdb_sources.json',
    sigma=3.6,              # Gaussian width (Å)
    cutoff=20.42,           # Cutoff distance (Å)
    debug_dir='debug_pdbs'  # Optional: save debug files
)

calculator.run(output_dir='desolvation_results')
```

Or run directly:
```bash
python desolvation_volume_calculator.py
```

## Outputs

### Main Results Directory
```
desolvation_results/
├── desolvation_summary.csv          # Summary statistics by atom type and source
├── desolvation_heatmap.png          # Heatmap of maximum desolvation values
└── {atom_type}/
    ├── {atom_type}_detailed_data.csv      # Per-atom desolvation values
    └── {atom_type}_distribution.png       # Distribution plots by source
```

### Debug Directory (if enabled)
```
debug_pdbs/
├── meeko_pdbs/                       # PDBs after Meeko processing
├── rdkit_pdbs/                       # RDKit molecule PDBs
└── pymol_commands/
    ├── {id}_pymol.pml                # PyMOL visualization script
    └── {id}_colorbar.png             # Viridis colorbar for desolvation values
```

## PyMOL Visualization

Debug mode generates PyMOL scripts that create two types of visualizations:

1. **Atom Type Objects**: Separate objects for each SMARTS-defined atom type
2. **Desolvation Coloring**: `desolv_colored` object with atoms colored by desolvation volume using a continuous viridis gradient

To use:
```python
# In PyMOL
load your_structure.pdb
run {id}_pymol.pml
```

The colorbar PNG shows the desolvation value scale.

## Dependencies

```
rdkit
numpy
pandas
matplotlib
seaborn
tqdm
meeko
prody
```


