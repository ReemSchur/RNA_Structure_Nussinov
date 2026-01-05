import sys
from nussinov_algo import NussinovRNA
from rna_utils import calculate_edit_distance

# --- Main Interactive Execution ---
if __name__ == "__main__":
    print("="*60)
    print("   GENERIC RNA STRUCTURE COMPARISON TOOL")
    print("="*60)

    # 1. Get Sequence and Calculate NJ
    seq_input = input("1. Enter RNA Sequence:\n   >> ").strip()
    if not seq_input:
        print("Error: Empty sequence.")
        sys.exit()
        
    print(f"\n... Running Nussinov Algorithm on {len(seq_input)} bases ...")
    solver = NussinovRNA(seq_input)
    nj_structure = solver.predict()
    pair_count = solver.get_pair_count()
    
    print(f"   -> NJ Structure Calculated: {nj_structure}")
    print(f"   -> Number of Base Pairs: {pair_count}")

    # 2. Get Reference Structure
    print("\n2. Enter BIOLOGICAL REFERENCE Structure (Dot-Bracket):")
    ref_structure = input("   >> ").strip()

    # 3. Get mfold Structure
    print("\n3. Enter MFOLD Structure (Dot-Bracket):")
    mfold_structure = input("   >> ").strip()

    # 4. Perform Comparisons
    if ref_structure and mfold_structure:
        print("\n" + "-"*60)
        print("COMPARISON RESULTS")
        print("-" * 60)

        # Comparison A: NJ vs Reference
        dist_nj = calculate_edit_distance(nj_structure, ref_structure)
        
        # Comparison B: mfold vs Reference
        dist_mfold = calculate_edit_distance(mfold_structure, ref_structure)

        print(f"{'Algorithm':<15} | {'Distance to Reference':<25} | {'Conclusion'}")
        print("-" * 60)
        print(f"{'Nussinov (NJ)':<15} | {dist_nj:<25} | {'High Error' if dist_nj > 10 else 'Good Match'}")
        print(f"{'mfold (Energy)':<15} | {dist_mfold:<25} | {'High Error' if dist_mfold > 10 else 'Good Match'}")
        print("-" * 60)
    else:
        print("\nError: Missing structure inputs. Cannot compare.")
