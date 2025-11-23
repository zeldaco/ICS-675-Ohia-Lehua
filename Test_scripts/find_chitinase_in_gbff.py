#!/usr/bin/env python3
"""
Scan entire .gbff file to see what feature types exist
"""

def analyze_gbff_structure(gbff_file):
    """Analyze what's actually in the GBFF file"""
    
    feature_types = {}
    total_lines = 0
    in_features = False
    in_origin = False
    
    print("Analyzing GBFF file structure...")
    
    with open(gbff_file, 'r') as f:
        for line in f:
            total_lines += 1
            
            # Track sections
            if line.startswith('FEATURES'):
                in_features = True
                print(f"Found FEATURES section at line {total_lines}")
            
            elif line.startswith('ORIGIN'):
                in_features = False
                in_origin = True
                print(f"Found ORIGIN section at line {total_lines}")
            
            # Count feature types in FEATURES section
            elif in_features and line.startswith('     '):
                # Feature lines start with 5 spaces, then feature type
                parts = line.strip().split()
                if parts:
                    feature_type = parts[0]
                    feature_types[feature_type] = feature_types.get(feature_type, 0) + 1
            
            # Look for chitinase anywhere in the file
            if 'chitinase' in line.lower():
                print(f"CHITINASE MENTION at line {total_lines}: {line.strip()[:100]}")
    
    print(f"\nTotal lines in file: {total_lines:,}")
    print(f"\nFeature types found:")
    
    if feature_types:
        for feature, count in sorted(feature_types.items()):
            print(f"   {feature}: {count:,}")
    else:
        print("   NO FEATURES FOUND")
    
    return feature_types

def main():
    gbff_file = r"Genome_Data\Metrosideros\metrosideros.gbff"
    
    try:
        feature_types = analyze_gbff_structure(gbff_file)
        
        print("\n" + "="*60)
        
        # Check for gene annotations
        if 'gene' in feature_types or 'CDS' in feature_types:
            print("This file HAS gene annotations!")
            print("   → Your original script should work, let me fix it...")
        else:
            print("This file DOES NOT have gene annotations")
            print("   → You'll need to use tBLASTn approach instead")
            print("\nThe file contains:")
            print("  • Assembly gaps (scaffold structure)")
            print("  • DNA sequence")
            print("  • But NO protein-coding gene predictions")
            print("\nNext step: Use tBLASTn to search the DNA for chitinase")
        
    except FileNotFoundError:
        print(f"Could not find file: {gbff_file}")
        print("   Make sure the path is correct!")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
    