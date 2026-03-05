import cobrakbase
from cobrakbase.kbaseapi import KBaseAPI
from tqdm import tqdm
import os

# Initialize KBase API
kbase = KBaseAPI()

def download_assembly(kbase, obj_name, workspace_id, output_dir):
    """Download an assembly's fasta file from KBase"""
    # Get the assembly object to find the fasta handle reference
    assembly_data = kbase.get_object(obj_name, workspace_id)

    if assembly_data is None:
        print(f"  Could not get assembly data for {obj_name}")
        return False

    # Assembly objects have a fasta_handle_ref that points to the actual fasta file
    if 'fasta_handle_ref' in assembly_data:
        handle_ref = assembly_data['fasta_handle_ref']
        output_path = os.path.join(output_dir, f"{obj_name}.fa")
        try:
            kbase.download_file_from_kbase2(handle_ref, output_path)
            return True
        except Exception as e:
            print(f"  Error downloading {obj_name}: {e}")
            return False
    else:
        print(f"  No fasta_handle_ref found for {obj_name}")
        return False

# Narrative 147022 -> MUCC assemblies
print("Pulling assemblies from narrative 147022 (MUCC)...")
os.makedirs("mucc_assemblies", exist_ok=True)
mucc_count = 0
# Filter by Assembly type directly - much faster than iterating all objects
assemblies_147022 = kbase.list_objects(147022, object_type='KBaseGenomeAnnotations.Assembly')
print(f"Found {len(assemblies_147022)} Assembly objects in workspace 147022")
for o in tqdm(assemblies_147022):
    if download_assembly(kbase, o[1], 147022, "mucc_assemblies"):
        mucc_count += 1

print(f"Downloaded {mucc_count} assemblies to mucc_assemblies/")

# Narrative 145226 -> GROW assemblies
print("\nPulling assemblies from narrative 145226 (GROW)...")
os.makedirs("grow_assemblies", exist_ok=True)
grow_count = 0
# Filter by Assembly type directly - much faster than iterating all objects
assemblies_145226 = kbase.list_objects(145226, object_type='KBaseGenomeAnnotations.Assembly')
print(f"Found {len(assemblies_145226)} Assembly objects in workspace 145226")
for o in tqdm(assemblies_145226):
    if download_assembly(kbase, o[1], 145226, "grow_assemblies"):
        grow_count += 1

print(f"Downloaded {grow_count} assemblies to grow_assemblies/")

print("\nDone!")
