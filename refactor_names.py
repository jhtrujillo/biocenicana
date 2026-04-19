import os
import subprocess

# Full mapping for class renames
CLASS_MAP = {
    "GeneDosis": "AlleleDosageCalculator",
    "SimilitudGeneticaCCdistTargeted": "TargetedGeneticSimilarity",
    "GenDosisTargeted": "TargetedAlleleDosage",
    "GenerarAlelos": "AlleleGenerator",
    "CompararDosisHuellaVsTargeted": "CompareDosageFingerprintVsTargeted",
    "FiltrosAlelologicosAUS": "AllelicFiltersAus",
    "FiltrarVCF": "VcfFilter",
    "OrdenarVCFporIndividuos": "SortVcfBySample",
    "SeleccionarAbanicoDosis": "SelectDosageRange",
    "VCFgetfilterprint": "VcfFilterPrint",
    "VCFgetHaplotipes": "VcfHaplotypeExtractor",
    "GeneDosisRapidGenomicVCF": "RapidGenomicAlleleDosage",
    "GenDosisTargetedCommand": "TargetedAlleleDosageCommand",
    "SimilitudGeneticaCommand": "GeneticSimilarityCommand",
    "Carolina": "CarolinaPipeline",
    "VcfToTabTargeted": "TargetedVcfToTab",
    "VcfToStructure2": "VcfToStructure",
    "FileUtilsV1": "LegacyFileUtilsV1",
    "FileUtils2": "LegacyFileUtilsV2",
    "Protocol_notes": "ProtocolNotes"
}

# Mapping for specific methods
METHOD_MAP = {
    "genDosisAlelicas": "computeAlleleDosage",
    "generarDosisTranspuesta": "computeDosageTransposed",
    "metodoImputar": "imputationMethod",
    "imputar": "impute"
}

SRC_DIR = "src/main/java/org/cenicana/bio"

# Step 1: Git MV files explicitly
def rename_files():
    for root, dirs, files in os.walk(SRC_DIR):
        for old_file in files:
            if not old_file.endswith(".java"):
                continue
            
            old_class = old_file[:-5]
            if old_class in CLASS_MAP:
                new_class = CLASS_MAP[old_class]
                old_path = os.path.join(root, old_file)
                new_path = os.path.join(root, new_class + ".java")
                print(f"Renaming file: {old_path} -> {new_path}")
                subprocess.run(["git", "mv", old_path, new_path], check=True)

# Step 2: Search & Replace inside all Java files
def update_references():
    for root, dirs, files in os.walk(SRC_DIR):
        for file in files:
            if not file.endswith(".java"):
                continue
                
            file_path = os.path.join(root, file)
            with open(file_path, "r", encoding="utf-8") as f:
                content = f.read()

            original_content = content
            
            # String replacement for Classes
            for old_name, new_name in CLASS_MAP.items():
                content = content.replace(old_name, new_name)
            
            # String replacement for Methods
            for old_meth, new_meth in METHOD_MAP.items():
                content = content.replace(old_meth, new_meth)

            if original_content != content:
                print(f"Updating references inside {file}")
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(content)

if __name__ == "__main__":
    rename_files()
    update_references()
    print("Done refactoring.")
