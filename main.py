import sys
from nussinov_algo import NussinovRNA
from rna_utils import calculate_rnadistance

# --- Main Interactive Execution ---
if __name__ == "__main__":
    print("="*60)
    print("   GENERIC RNA STRUCTURE COMPARISON TOOL")
    print("   (Uses ViennaRNA RNAdistance executable)")
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
    print(f"   -> Number of Base Pairs:    {pair_count}")

    # 2. Get Reference Structure
    print("\n2. Enter BIOLOGICAL REFERENCE Structure (Dot-Bracket):")
    ref_structure = input("   >> ").strip()

    # 3. Get mfold Structure
    print("\n3. Enter MFOLD Structure (Dot-Bracket):")
    mfold_structure = input("   >> ").strip()

    # 4. Perform Comparisons
    if ref_structure and mfold_structure:
        print("\n" + "-"*60)
        print("COMPARISON RESULTS (Metric: Tree Edit Distance)")
        print("-" * 60)

        # Comparison A: NJ vs Reference (Using Default Tree Edit Distance)
        dist_nj = calculate_rnadistance(nj_structure, ref_structure, method='default')
        
        # Comparison B: mfold vs Reference (Using Default Tree Edit Distance)
        dist_mfold = calculate_rnadistance(mfold_structure, ref_structure, method='default')

        # Handle cases where RNAdistance might fail or return None
        nj_display = dist_nj if dist_nj is not None else "Error"
        mfold_display = dist_mfold if dist_mfold is not None else "Error"

        print(f"{'Algorithm':<20} | {'Tree Edit Distance':<20} | {'Notes'}")
        print("-" * 60)
        print(f"{'Nussinov (NJ)':<20} | {str(nj_display):<20} | {'Your Implementation'}")
        print(f"{'mfold (Energy)':<20} | {str(mfold_display):<20} | {'Energy Minimization'}")
        print("-" * 60)
        
        print("\nNote: Lower distance indicates better structural similarity.")
    else:
        print("\nError: Missing structure inputs. Cannot compare.")