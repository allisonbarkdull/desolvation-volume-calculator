#!/usr/bin/env python
"""
Desolvation Volume Calculator
"""

import json
import numpy as np
from pathlib import Path
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import pandas as pd
import time
import tempfile
from meeko.polymer import Polymer, ResidueChemTemplates
from meeko import MoleculePreparation
import prody



def sdf_to_template_dict(sdf_path: str, resname: str) -> dict:
    """Convert SDF file to residue template dictionary for meekoo."""
    mol = Chem.SDMolSupplier(str(sdf_path), removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Could not read molecule from {sdf_path}")
    
    mol = Chem.AddHs(mol)
    smiles = Chem.MolToSmiles(mol, allHsExplicit=False)
    atom_names = [f"A{i+1}" for i in range(mol.GetNumAtoms())]
    
    return {
        "ambiguous": {resname: [resname]},
        "residue_templates": {
            resname: {
                "smiles": smiles,
                "atom_name": atom_names,
                "link_labels": {}
            }
        }
    }


def smiles_to_template_dict(smiles: str, resname: str) -> dict:
    """Convert SMILES string to residue template dictionary."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    smiles_with_h = Chem.MolToSmiles(mol, allHsExplicit=False)
    atom_names = [f"A{i+1}" for i in range(mol.GetNumAtoms())]
    
    return {
        "ambiguous": {resname: [resname]},
        "residue_templates": {
            resname: {
                "smiles": smiles_with_h,
                "atom_name": atom_names,
                "link_labels": {}
            }
        }
    }


def strip_dum_residues(pdb_path):
    """Remove DUM residues from PDB file, need this if using outputs from my metadynamics 
    prody will not be happy."""
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    
    cleaned_lines = []
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            resname = line[17:20].strip()
            if resname == 'DUM':
                continue
        cleaned_lines.append(line)
    
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
    temp_file.writelines(cleaned_lines)
    temp_file.close()
    
    return temp_file.name

#plot settings
sns.set_style("whitegrid")
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10


class DesolvationCalculatorMeeko:
    """Calculate desolvation volumes from pdbs"""
    
    def __init__(self, atom_types_file, pdb_sources_file, sigma=3.6, cutoff=20.42, debug_dir=None):

        
        self.sigma = sigma
        self.cutoff = cutoff
        self.cutoff_sq = cutoff ** 2
        self.inv_sigma_sq = 1.0 / (sigma ** 2)
        self.debug_dir = Path(debug_dir) if debug_dir else None
        
        if self.debug_dir:
            self.debug_dir.mkdir(parents=True, exist_ok=True)
            self.meeko_pdb_dir = self.debug_dir / "meeko_pdbs"
            self.rdkit_pdb_dir = self.debug_dir / "rdkit_pdbs"
            self.pymol_cmd_dir = self.debug_dir / "pymol_commands"
            
            self.meeko_pdb_dir.mkdir(exist_ok=True)
            self.rdkit_pdb_dir.mkdir(exist_ok=True)
            self.pymol_cmd_dir.mkdir(exist_ok=True)
        
        # Load configurations
        with open(atom_types_file) as f:
            self.atom_types_config = json.load(f)
        
        with open(pdb_sources_file) as f:
            pdb_config = json.load(f)
            self.pdb_sources = pdb_config['sources']
            self.custom_templates = pdb_config.get('custom_templates', {})
        
        # Compile SMARTS patterns for different atom types
        self.patterns = {}
        for atom_type_name, config in self.atom_types_config['atom_types'].items():
            pattern = Chem.MolFromSmarts(config['smarts'])
            if pattern is None:
                raise ValueError(f"Invalid SMARTS for {atom_type_name}: {config['smarts']}")
            self.patterns[atom_type_name] = {
                'pattern': pattern,
                'volume': config['volume']
            }
        
        self.results = defaultdict(lambda: defaultdict(list))
        self.detailed_results = defaultdict(list)  # detailed data per atom
        self.coverage_stats = []
        self.chem_templates = self._initialize_templates()
        self.pdb_counter = 0  # Counter for unique PDB naming

    
    def _initialize_templates(self):
        """add residue templates to meeko from input json"""
        templates = ResidueChemTemplates.create_from_defaults()
        
        if 'sdf_templates' in self.custom_templates:
            print("\nLoading SDF templates:")
            for template_info in self.custom_templates['sdf_templates']:
                sdf_path = template_info['sdf_path']
                resname = template_info['resname']
                print(f"  - {resname} from {sdf_path}")
                try:
                    template_dict = sdf_to_template_dict(sdf_path, resname)
                    templates.add_dict(template_dict)
                except Exception as e:
                    print(f"    WARNING: Failed to load SDF template: {e}")
        
        if 'smiles_templates' in self.custom_templates:
            print("\nLoading SMILES templates:")
            for template_info in self.custom_templates['smiles_templates']:
                smiles = template_info['smiles']
                resname = template_info['resname']
                print(f"  - {resname} from SMILES: {smiles}")
                try:
                    template_dict = smiles_to_template_dict(smiles, resname)
                    templates.add_dict(template_dict)
                except Exception as e:
                    print(f"    WARNING: Failed to load SMILES template: {e}")
        
        return templates
    
    def save_rdkit_pdb(self, mol, unique_id):
        """Save RDKit molecule as PDB with unique identifier for debugging."""
        if not self.debug_dir:
            return None
        
        output_path = self.rdkit_pdb_dir / f"{unique_id}_rdkit.pdb"
        Chem.MolToPDBFile(mol, str(output_path))
        return output_path
    
    def generate_pymol_commands(self, unique_id, atom_assignments, rdkit_pdb_path, total_atoms, desolvation_values=None):
        """Generate PyMOL commands to visualize atom types and desolvation volumes for debugging."""
        if not self.debug_dir:
            return
        
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib.cm import ScalarMappable
        
        # Group atoms by type
        atoms_by_type = defaultdict(list)
        all_assigned_ids = set()
        
        for atom_idx, (atom_type_name, _) in atom_assignments.items():
            atom_id = atom_idx + 1  # 1-based indexing
            atoms_by_type[atom_type_name].append(atom_id)
            all_assigned_ids.add(atom_id)
        
        # Generate PyMOL script
        script_lines = [
            f"# PyMOL visualization for {unique_id}",
            f"#load {rdkit_pdb_path.name}",
            "hide everything",
            "show sticks",
            ""
        ]
        
        # ===== ATOM TYPE VISUALIZATION (EXISTING) =====
        script_lines.append("# ===== ATOM TYPE VISUALIZATION =====")
        script_lines.append("")
        
        # Create selections and show as spheres for each atom type
        for atom_type_name in sorted(atoms_by_type.keys()):
            atom_indices = atoms_by_type[atom_type_name]
            selection_name = f"type_{atom_type_name}"
            id_string = "+".join(map(str, atom_indices))
            
            script_lines.extend([
                f"# {atom_type_name} ({len(atom_indices)} atoms)",
                f"create {selection_name}, id {id_string}",
                f"show spheres, {selection_name}",
                ""
            ])
        
        # Create object for unassigned atoms
        if len(all_assigned_ids) < total_atoms:
            unassigned_ids = sorted(set(range(1, total_atoms + 1)) - all_assigned_ids)
            if unassigned_ids:
                id_string = "+".join(map(str, unassigned_ids))
                script_lines.extend([
                    f"# Unassigned atoms ({len(unassigned_ids)} atoms)",
                    f"create unassigned, id {id_string}",
                    "show spheres, unassigned",
                    ""
                ])
        
        # ===== DESOLVATION VOLUME VISUALIZATION =====
        if desolvation_values is not None:
            script_lines.extend([
                "",
                "# ===== DESOLVATION VOLUME VISUALIZATION =====",
                ""
            ])
            
            # Get desolvation values for assigned atoms
            desolv_dict = {}  # atom_id (1-based) -> desolvation value
            for atom_idx, desolv_val in desolvation_values.items():
                desolv_dict[atom_idx + 1] = desolv_val  # Convert to 1-based
            
            if desolv_dict:
                # Get min/max for normalization
                min_val = min(desolv_dict.values())
                max_val = max(desolv_dict.values())
                
                # Create a copy for desolvation visualization
                assigned_ids_str = "+".join(map(str, sorted(desolv_dict.keys())))
                script_lines.extend([
                    f"# Create copy colored by desolvation volume",
                    f"create desolv_colored, id {assigned_ids_str}",
                    "show spheres, desolv_colored",
                    f"# Desolvation range: {min_val:.2f} to {max_val:.2f}",
                    ""
                ])
                
                # Use viridis colormap
                cmap = plt.cm.viridis
                norm = mcolors.Normalize(vmin=min_val, vmax=max_val)
                
                # Color each atom individually
                for atom_id, desolv_val in sorted(desolv_dict.items()):
                    # Get RGB from viridis
                    rgba = cmap(norm(desolv_val))
                    r, g, b = rgba[0], rgba[1], rgba[2]
                    
                    script_lines.append(
                        f"set_color desolv_{atom_id}, [{r:.4f}, {g:.4f}, {b:.4f}]"
                    )
                    script_lines.append(
                        f"color desolv_{atom_id}, desolv_colored and id {atom_id}"
                    )
                
                script_lines.append("")
                
                # Create colorbar using matplotlib
                fig, ax = plt.subplots(figsize=(1.5, 6))
                fig.subplots_adjust(right=0.3)
                
                sm = ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                
                cbar = plt.colorbar(sm, cax=ax)
                cbar.set_label('Desolvation Volume', fontsize=12, fontweight='bold')
                
                # Save colorbar
                colorbar_path = self.pymol_cmd_dir / f"{unique_id}_colorbar.png"
                plt.savefig(colorbar_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                script_lines.extend([
                    f"# Colorbar saved to: {colorbar_path.name}",
                    ""
                ])
        
        script_lines.extend([
            "deselect",
            "zoom",
        ])
        
        # Write script
        script_path = self.pymol_cmd_dir / f"{unique_id}_pymol.pml"
        with open(script_path, 'w') as f:
            f.write('\n'.join(script_lines))
        
        print(f"  PyMOL script: {script_path}")

    
    def load_pdb_with_meeko(self, pdb_path, unique_id):
        """Load PDB using Meeko."""
        cleaned_pdb_path = None
        
        try:
            mk_prep = MoleculePreparation()
            cleaned_pdb_path = strip_dum_residues(pdb_path)
            
            input_obj = prody.parsePDB(cleaned_pdb_path, altloc="all")
            
            start = time.time()
            polymer = Polymer.from_prody(
                input_obj,
                self.chem_templates,
                mk_prep,
                set_template={},
                allow_bad_res=True,
                wanted_altloc=None,
                default_altloc=None,
            )
            polymer_time = time.time() - start
            
            valid_monomers = polymer.get_valid_monomers()
            if not valid_monomers:
                print(f"  Warning: No valid monomers found")
                return None
            
            print(f"  Polymer creation: {polymer_time:.3f}s ({len(valid_monomers)} residues)")
            
            start = time.time()
            mol = polymer.stitch() #i rly think this should be named to_rdkit_mol
            stitch_time = time.time() - start
            print(f"  Stitching: {stitch_time:.3f}s")
            
            # Save Meeko PDB with unique identifier
            if self.debug_dir:
                debug_path = self.meeko_pdb_dir / f"{unique_id}_meeko.pdb"
                pdb_string = polymer.to_pdb()
                with open(debug_path, 'w') as f:
                    f.write(pdb_string)
            
            if mol.GetNumAtoms() > 0:
                mol = Chem.RemoveHs(mol) #could keep hydrogens and change my smarts pattern to ignore mayb
                return mol
            else:
                return None
                
        except Exception as e:
            print(f"  Error loading {pdb_path}: {e}")
            return None
        finally:
            if cleaned_pdb_path:
                try:
                    Path(cleaned_pdb_path).unlink()
                except:
                    pass
    
    def assign_atom_types(self, mol):
        """Assign atom types based on SMARTS patterns."""
        Chem.SanitizeMol(mol)
        atom_assignments = {}
        matched_atoms = set()
        n_atoms = mol.GetNumAtoms()

        for atom_type_name, pattern_info in self.patterns.items():
            pattern = pattern_info['pattern']
            volume = pattern_info['volume']
            matches = mol.GetSubstructMatches(pattern, maxMatches=n_atoms)
            
            for match in matches:
                for atom_idx in match:
                    if atom_idx not in matched_atoms:
                        atom_assignments[atom_idx] = (atom_type_name, volume)
                        matched_atoms.add(atom_idx)
        
        return atom_assignments
    
    def calculate_desolvation_volume(self, coords, volumes, atom_idx):
        """Calculate desolvation volume for a single atom."""
        atom_coord = coords[atom_idx]
        diff = coords - atom_coord
        dist_sq = np.sum(diff ** 2, axis=1)
        mask = (dist_sq < self.cutoff_sq) & (dist_sq > 0)
        
        if not np.any(mask):
            return 0.0
        
        dist_sq_masked = dist_sq[mask]
        volumes_masked = volumes[mask]
        gaussian_weights = np.exp(-0.5 * dist_sq_masked * self.inv_sigma_sq)
        
        return np.sum(volumes_masked * gaussian_weights)
    
    def process_pdb(self, pdb_path, source_name):
        """Process a single PDB file - covert to rdkit mol,
          calcucate desolvation for each atom, 
          write out debug pdbs and pymol commnds to visulize results."""
        # Create unique identifier for this PDB
        self.pdb_counter += 1
        pdb_basename = Path(pdb_path).stem
        unique_id = f"{self.pdb_counter:04d}_{source_name}_{pdb_basename}"
        
        mol = self.load_pdb_with_meeko(pdb_path, unique_id)
        
        if mol is None:
            print(f"  Could not load {pdb_path}")
            return None
        
        # Get coordinates
        conf = mol.GetConformer()
        n_atoms = mol.GetNumAtoms()
        coords = np.zeros((n_atoms, 3))
        
        for i in range(n_atoms):
            pos = conf.GetAtomPosition(i)
            coords[i] = [pos.x, pos.y, pos.z]
        
        # Assign atom types
        atom_assignments = self.assign_atom_types(mol)
        
        coverage = len(atom_assignments) / n_atoms * 100
        self.coverage_stats.append({
            'source': source_name,
            'pdb': Path(pdb_path).name,
            'coverage_percent': coverage,
            'matched_atoms': len(atom_assignments),
            'total_atoms': n_atoms
        })
        
        print(f"  Coverage: {coverage:.1f}% ({len(atom_assignments)}/{n_atoms} atoms)")
        
        if not atom_assignments:
            print(f"  Warning: No atoms matched patterns")
            # Still save debug outputs even with no assignments
            if self.debug_dir:
                rdkit_pdb_path = self.save_rdkit_pdb(mol, unique_id)
                self.generate_pymol_commands(unique_id, atom_assignments, rdkit_pdb_path, n_atoms)
            return None
        
        # Calculate desolvation volumes
        matched_indices = sorted(atom_assignments.keys())
        matched_coords = coords[matched_indices]
        matched_volumes = np.array([atom_assignments[i][1] for i in matched_indices])
        
        pdb_results = defaultdict(list)
        desolvation_dict = {}  # Store for PyMOL visualization (0-based indexing)
        
        for local_idx, global_idx in enumerate(tqdm(matched_indices, 
                                                    desc="  Calculating desolvation",
                                                    leave=False)):
            atom_type = atom_assignments[global_idx][0]
            desolvation = self.calculate_desolvation_volume(
                matched_coords, matched_volumes, local_idx
            )
            pdb_results[atom_type].append(desolvation)
            self.results[atom_type][source_name].append(desolvation)
            
            # Store desolvation value for this atom (0-based index)
            desolvation_dict[global_idx] = desolvation
            
            # Store detailed data
            self.detailed_results[atom_type].append({
                'source': source_name,
                'input_pdb': Path(pdb_path).name,
                'pdb_path': pdb_path,
                'atom_id': global_idx,
                'desolvation_volume': desolvation
            })
        
        # Save debug outputs with desolvation values
        if self.debug_dir:
            rdkit_pdb_path = self.save_rdkit_pdb(mol, unique_id)
            self.generate_pymol_commands(unique_id, atom_assignments, rdkit_pdb_path, 
                                        n_atoms, desolvation_values=desolvation_dict)
        
        return pdb_results
    
    def process_all_pdbs(self):
        """Process all PDB files."""
        print("\n" + "="*80)
        print("PROCESSING PDB FILES")
        print("="*80 + "\n")
        
        for source_name, pdb_list in self.pdb_sources.items():
            print(f"\nSource: {source_name}")
            print("-" * 80)
            
            for pdb_path in pdb_list:
                print(f"\nProcessing: {pdb_path}")
                self.process_pdb(pdb_path, source_name)
        
        self.print_coverage_summary()
    
    def print_coverage_summary(self):
        """Print coverage statistics."""
        print("\n" + "="*80)
        print("COVERAGE SUMMARY")
        print("="*80 + "\n")
        
        if not self.coverage_stats:
            print("No files successfully processed")
            return
        
        df = pd.DataFrame(self.coverage_stats)
        print(df.to_string(index=False))
        
        print("\n\nBy Source:")
        summary = df.groupby('source').agg({
            'coverage_percent': 'mean',
            'matched_atoms': 'sum',
            'total_atoms': 'sum'
        })
        summary['overall_coverage'] = (summary['matched_atoms'] / summary['total_atoms'] * 100)
        print(summary)
    
    def create_output_directory(self, base_dir="desolvation_results"):
        """Create output directory."""
        self.output_dir = Path(base_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        for atom_type_name in self.patterns.keys():
            (self.output_dir / atom_type_name).mkdir(exist_ok=True)
    
    def save_detailed_csvs(self):
        """Save detailed CSV files for each atom type."""
        print("\n" + "="*80)
        print("SAVING DETAILED CSV FILES")
        print("="*80 + "\n")
        
        for atom_type_name in tqdm(self.detailed_results.keys(), desc="Saving CSVs"):
            df = pd.DataFrame(self.detailed_results[atom_type_name])
            csv_path = self.output_dir / atom_type_name / f'{atom_type_name}_detailed_data.csv'
            df.to_csv(csv_path, index=False)
            print(f"  {atom_type_name}: {len(df)} records -> {csv_path}")
    
    def plot_distributions(self):
        """Create distribution plots."""
        print("\n" + "="*80)
        print("GENERATING PLOTS")
        print("="*80 + "\n")
        
        for atom_type_name in tqdm(self.results.keys(), desc="Creating plots"):
            fig, ax = plt.subplots(figsize=(10, 6))
            
            all_data = []
            sources = []
            
            for source_name, values in self.results[atom_type_name].items():
                all_data.extend(values)
                sources.extend([source_name] * len(values))
            
            if not all_data:
                continue
            
            df = pd.DataFrame({
                'Desolvation Volume': all_data,
                'Source': sources
            })
            
            sns.histplot(data=df, x='Desolvation Volume', hue='Source',
                        bins=50, alpha=0.6, ax=ax, stat='density')
            
            for source in df['Source'].unique():
                source_data = df[df['Source'] == source]['Desolvation Volume']
                sns.kdeplot(data=source_data, ax=ax, linewidth=2, label=f'{source} KDE')
            
            ax.set_xlabel('Desolvation Volume', fontsize=12, fontweight='bold')
            ax.set_ylabel('Density', fontsize=12, fontweight='bold')
            ax.set_title(f'{atom_type_name} (n={len(all_data):,})', fontsize=14, fontweight='bold')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(self.output_dir / atom_type_name / f'{atom_type_name}_distribution.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    def create_summary_table(self):
        """Create summary table."""
        print("\n" + "="*80)
        print("CREATING SUMMARY")
        print("="*80 + "\n")
        
        summary_data = []
        
        for atom_type_name in sorted(self.results.keys()):
            row = {'Atom Type': atom_type_name}
            
            for source_name in sorted(self.pdb_sources.keys()):
                if source_name in self.results[atom_type_name]:
                    values = self.results[atom_type_name][source_name]
                    if values:
                        row[f'{source_name}_max'] = max(values)
                        row[f'{source_name}_mean'] = np.mean(values)
                        row[f'{source_name}_std'] = np.std(values)
                        row[f'{source_name}_n'] = len(values)
            
            summary_data.append(row)
        
        df = pd.DataFrame(summary_data)
        
        csv_path = self.output_dir / 'desolvation_summary.csv'
        df.to_csv(csv_path, index=False)
        print(f"Summary saved: {csv_path}")
        print("\n" + df.to_string(index=False))
        
        self.create_heatmap(df)
        return df
    
    def create_heatmap(self, df):
        """Create heatmap of maximum values."""
        max_cols = [col for col in df.columns if col.endswith('_max')]
        if not max_cols:
            return
        
        heatmap_data = df[['Atom Type'] + max_cols].set_index('Atom Type')
        heatmap_data.columns = [col.replace('_max', '') for col in heatmap_data.columns]
        
        fig, ax = plt.subplots(figsize=(12, 8))
        sns.heatmap(heatmap_data, annot=True, fmt='.1f', cmap='YlOrRd',
                   cbar_kws={'label': 'Max Desolvation Volume'},
                   ax=ax, linewidths=0.5)
        
        ax.set_xlabel('Source', fontsize=12, fontweight='bold')
        ax.set_ylabel('Atom Type', fontsize=12, fontweight='bold')
        ax.set_title('Maximum Desolvation Volumes', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'desolvation_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def run(self, output_dir="desolvation_results"):
        """Run complete analysis."""
        self.create_output_directory(output_dir)
        self.process_all_pdbs()
        self.save_detailed_csvs()
        self.plot_distributions()
        self.create_summary_table()
        
        print("\n" + "="*80)
        print("COMPLETE!")
        print("="*80)
        print(f"\nResults: {self.output_dir.absolute()}")
        
        if self.debug_dir:
            print(f"\nDebug files:")
            print(f"  Meeko PDFs: {self.meeko_pdb_dir.absolute()}")
            print(f"  RDKit PDFs: {self.rdkit_pdb_dir.absolute()}")
            print(f"  PyMOL scripts: {self.pymol_cmd_dir.absolute()}")


def main():
    calculator = DesolvationCalculatorMeeko(
        atom_types_file='atom_types.json',
        pdb_sources_file='pdb_sources.json',
        sigma=3.6,
        cutoff=20.42,
        debug_dir='debug_pdbs'
    )
    
    calculator.run()


if __name__ == '__main__':
    main()